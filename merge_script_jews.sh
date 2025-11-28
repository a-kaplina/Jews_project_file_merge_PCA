#!/bin/bash
set -e

if [ $# -lt 1 ]; then
    echo "Usage: $0 <vcf_folder> [sample_to_remove]"
    exit 1
fi

VCF_DIR="$1"
OUT_DIR="$VCF_DIR/output"
mkdir -p "$OUT_DIR"
cp "$VCF_DIR/jews_and_new_groups.csv" "$OUT_DIR"
REMOVE_SAMPLE="${2:-nx1575}"
MERGED="$OUT_DIR/merged.vcf.gz"
SORTED_MERGED="$OUT_DIR/merged_sorted.vcf"
FILTERED="$OUT_DIR/filtered.vcf.gz"
stats_file="$OUT_DIR/normalization_stats.tsv"
echo -e "Sample\tTotal\tSplit\tJoined\tRealigned\tSkipped" > "$stats_file"

echo "STEP 1: SORT, INDEX AND NORMALIZE ALL VCFs"
# Find *.vcf.gz
files=( "$VCF_DIR"/*.vcf.gz )

# If none exist - stop
if [ ! -e "${files[0]}" ]; then
    echo "ERROR: No .vcf.gz files found in $VCF_DIR"
    exit 1
fi

for f in "${files[@]}"; do
    base=$(basename "$f" .vcf.gz)
    echo "→ Processing $base"

    # Sort, normalize
    norm_output=$(bcftools sort "$f" | bcftools norm -m-any -Oz -o "$OUT_DIR/${base}.sorted_norm.vcf.gz" 2>&1)

    # indexation
    bcftools index "$OUT_DIR/${base}.sorted_norm.vcf.gz"

    # During normalization outputs  Lines with stats - save it
    line=$(echo "$norm_output" | grep "Lines.*total/split/joined/realigned/skipped" | tail -n1)

    numbers=$(echo "$line" | sed -E 's/.*total\/split\/joined\/realigned\/skipped:[[:space:]]*([0-9/]+).*/\1/')

    if [ -z "$numbers" ]; then
        total=0; split=0; joined=0; realigned=0; skipped=0
    else
        IFS='/' read -r total split joined realigned skipped <<< "$numbers"
    fi

    echo -e "${base}\t${total}\t${split}\t${joined}\t${realigned}\t${skipped}" >> "$stats_file"
done

echo "STEP 2: MERGE of VCFs"
bcftools merge "$OUT_DIR"/*.sorted_norm.vcf.gz -Oz -o "$MERGED"
bcftools index "$MERGED"

echo "STEP 3: SORT MERGED FILE"
bcftools sort "$MERGED" -Ov -o "$SORTED_MERGED"

echo "STEP 4: REMOVE KNOWN RELATIVE SAMPLE nx1575: $REMOVE_SAMPLE"
vcftools --remove-indv "$REMOVE_SAMPLE" \
    --vcf "$SORTED_MERGED" --recode --out $OUT_DIR/no_$REMOVE_SAMPLE

EXCLUDED="$OUT_DIR/no_${REMOVE_SAMPLE}.recode.vcf"

echo "STEP 5: FILTERING OF VCF"
vcftools --vcf "$EXCLUDED" \
    --remove-indels --maf 0.05 --max-missing 0.9 \
    --recode --stdout | gzip -c > "$FILTERED"

echo "STEP 6: LD-PRUNING"
plink --vcf "$FILTERED" \
    --double-id --allow-extra-chr --set-missing-var-ids @:# \
    --indep-pairwise 50 10 0.1 \
    --out "$OUT_DIR/pruned"

echo "STEP 7: PCA IN PLINK"
plink --vcf "$FILTERED" \
    --double-id --allow-extra-chr --set-missing-var-ids @:# \
    --extract $OUT_DIR/pruned.prune.in \
    --make-bed --pca \
    --out "$OUT_DIR/pruned"

echo "=== STEP 8: PYTHON PCA → PDF + ELLIPSES + TABLE ==="

cat << 'EOF' > plot_pca.py
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import sys
import os
from matplotlib.patches import Ellipse
from matplotlib.path import Path
from scipy.stats import chi2

OUT_DIR = sys.argv[1]

# ============================================================
# Functions for confidence ellipses
# ============================================================
def compute_confidence_ellipse(x, y, confidence=0.95):
    cov = np.cov(x, y)
    lambda_, v = np.linalg.eig(cov)
    order = lambda_.argsort()[::-1]
    lambda_ = lambda_[order]
    v = v[:, order]
    chi2_val = chi2.ppf(confidence, df=2)
    width = 2 * np.sqrt(lambda_[0] * chi2_val)
    height = 2 * np.sqrt(lambda_[1] * chi2_val)
    angle = np.degrees(np.arctan2(*v[:, 0][::-1]))
    return width, height, angle

def ellipse_polygon(cx, cy, width, height, angle_deg, n_points=200):
    t = np.linspace(0, 2*np.pi, n_points)
    coords = np.column_stack([width/2 * np.cos(t), height/2 * np.sin(t)])
    theta = np.radians(angle_deg)
    R = np.array([[np.cos(theta), -np.sin(theta)],
                  [np.sin(theta),  np.cos(theta)]])
    coords = coords @ R.T
    coords[:,0] += cx
    coords[:,1] += cy
    return coords

def assign_samples_to_ellipse(x, y, ellipse_poly):
    poly_path = Path(ellipse_poly)
    return poly_path.contains_points(np.column_stack([x, y]))

# ============================================================
# Load PCA results
# ============================================================
cols = ['FID','IID'] + [f'PC{i}' for i in range(1,21)]

eigenvec_file = os.path.join(OUT_DIR, "pruned.eigenvec")

df = pd.read_csv(eigenvec_file, 
                 delim_whitespace=True, 
                 names=cols, 
                 header=None)

# ============================================================
# Update groups with new samples
# ============================================================
groups_file   = os.path.join(OUT_DIR, "jews_and_new_groups.csv")

if os.path.exists(groups_file):
    groups = pd.read_csv(groups_file)
else:
    groups = pd.DataFrame(columns=["IID", "groups"])

known_ids = set(groups["IID"].astype(str))
pca_ids = set(df["IID"].astype(str))
new_ids = pca_ids - known_ids

if len(new_ids) > 0:
    add_df = pd.DataFrame({"IID": list(new_ids), "groups": ["new"]*len(new_ids)})
    groups = pd.concat([groups, add_df], ignore_index=True)
    groups.to_csv(groups_file, index=False)

df = df.merge(groups, on="IID", how="left")
df["groups"] = df["groups"].fillna("new")

# ============================================================
# PCA Plot: PC1 vs PC3 with ellipses
# ============================================================
sns.set(style="whitegrid")
plt.figure(figsize=(10,8))

unique_groups = df["groups"].unique()
ellipse_dict = {}
ellipse_table_rows = []

ellipse_scale = 0.9 

# make ellipses for known Jewish populations
for grp in unique_groups:
    sub = df[df["groups"] == grp]
    plt.scatter(sub["PC1"], sub["PC3"], label=grp, s=120, alpha=0.85)

    if len(sub) > 2 and grp != "new":
        cx, cy = sub["PC1"].mean(), sub["PC3"].mean()
        width, height, angle = compute_confidence_ellipse(sub["PC1"], sub["PC3"])
        
        width  *= ellipse_scale
        height *= ellipse_scale
        
        ell = Ellipse((cx, cy), width, height, angle=angle,
                      edgecolor="black", fill=False, lw=1.5)
        plt.gca().add_patch(ell)
        poly = ellipse_polygon(cx, cy, width, height, angle)
        ellipse_dict[grp] = poly


for grp, poly in ellipse_dict.items():
    inside_mask = assign_samples_to_ellipse(df["PC1"], df["PC3"], poly)
    inside_samples = df[inside_mask]
    for _, row in inside_samples.iterrows():
        ellipse_table_rows.append({
            "ellipse_group": grp,
            "IID": row["IID"],
            "groups": row["groups"],
            "ellipse_built_on": "PC1 & PC3"
        })

# If new samples are outside ellipses → not_defined
new_samples = df[df["groups"] == "new"]
for _, row in new_samples.iterrows():
    assigned = False
    for poly in ellipse_dict.values():
        if assign_samples_to_ellipse([row["PC1"]], [row["PC3"]], poly)[0]:
            assigned = True
            break
    if not assigned:
        ellipse_table_rows.append({
            "ellipse_group": "not_defined",
            "IID": row["IID"],
            "groups": row["groups"],
            "ellipse_built_on": "PC1 & PC3"
        })

# Always annotate ALL 'new' samples on the PCA plot
for _, row in new_samples.iterrows():
    plt.annotate(
        row["IID"],
        (row["PC1"], row["PC3"]),
        textcoords="offset points",
        xytext=(6, 6),
        fontsize=7
    )


# Save table with groups annotations for new samples
ellipse_table = pd.DataFrame(ellipse_table_rows)
ellipse_table.to_csv(os.path.join(OUT_DIR, "samples_in_ellipses.csv"), index=False)

# Save plot with ellipses
plt.xlabel("PC1")
plt.ylabel("PC3")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "PCA_PC1_PC3_with_ellipses.pdf"))

# ============================================================
# PCA Plot 2: PC1 vs PC2 (without ellipses)
# ============================================================
plt.figure(figsize=(10,8))
for grp in unique_groups:
    sub = df[df["groups"] == grp]
    plt.scatter(sub["PC1"], sub["PC2"], label=grp, s=120, alpha=0.85)

plt.xlabel("PC1")
plt.ylabel("PC2")
plt.title("PCA: PC1 vs PC2 (no ellipses)")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "PCA_PC1_PC2.pdf"))

EOF

python3 plot_pca.py "$OUT_DIR"

echo "DONE → Files are in $OUT_DIR: PCA_PC1_PC2.pdf, PCA_PC1_PC3_with_ellipses.pdf, samples_in_ellipses.csv, normalization_stats.tsv"