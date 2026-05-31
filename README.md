# Jews_project_file_merge_PCA

Script for merging individual `vcf.gz` files with a reference `vcf.gz` file containing samples from Jewish populations, generating PCA plots, and assigning new samples to populations based on PCA results. Output files are saved in the `output` folder (created inside the directory containing the input `vcf.gz` files).

- The script normalizes (normalization statistics are saved in `normalization_stats.tsv`), indexes, and merges individual `vcf.gz` files with the reference Jewish population dataset (`vcf_jews_bgz.vcf.gz`).
- One related sample (`nx1575`) is removed from the reference dataset because it is a previously identified related individual (the only sample pair with PI-HAT > 0.25: PI_HAT = 0.5386 between samples `nx1575` and `kl9502`). Relatedness checking is **not** performed within this script.
- The merged VCF is filtered: all insertions and deletions are removed (only SNPs are retained). Variants are kept only if their minor allele frequency (MAF) is at least 5%, and if genotype data are available for at least 90% of samples.
- LD pruning is performed using a sliding window of 50 SNPs, a step size of 10 SNPs, and an r² threshold of 0.1.
- PCA is performed using PLINK.
- A Python script is generated that loads PCA results, assigns samples to groups using the `jews_and_new_groups.csv` table (IDs of new samples are added automatically with the group label `new`).
- PCA plots for PC1–PC2 and PC1–PC3 are generated. Population ellipses are drawn for known populations on the PC1–PC3 plot, and IDs of new samples are displayed on the plots.
- A table (`samples_in_ellipses.csv`) is generated indicating whether samples fall within population ellipses. New samples are labeled as `new`. If a new sample does not fall within any ellipse, it is labeled as `not_defined`.
- PCA plots are saved in PDF format (`PCA_PC1_PC2.pdf`, `PCA_PC1_PC3_with_ellipses.pdf`).

## Software used

- bcftools 1.20
- VCFtools 0.1.17
- PLINK v1.9.0-b.7.7 64-bit (22 Oct 2024) — https://www.cog-genomics.org/plink/1.9/
- Python 3.10.10 (libraries: pandas, numpy, matplotlib, seaborn, scipy)

## Running the script

```bash
merge_script_jews.sh ~/absolute/path
```

where `~/absolute/path` is the absolute path to the folder containing the VCF files to be merged.

## Input directory contents

The directory should contain:

- VCF files to be merged
- Reference population VCF file (`vcf_jews_bgz.vcf.gz`)
- Population annotation table (`jews_and_new_groups.csv`) containing sample IDs and Jewish population group labels corresponding to samples in `vcf_jews_bgz.vcf.gz`
 <img width="639" height="412" alt="520185374-5c0d87b3-b0dd-4fd0-85d5-613f55adfa1c" src="https://github.com/user-attachments/assets/33503e91-5f9f-4940-b2bc-1ebc3519fb01" />


## Output files

After execution, the `output` directory (created inside the sample directory) will contain:

- Merged and indexed VCF file (unfiltered merged dataset): `merged.vcf.gz` and `merged.vcf.gz.csi`
- Filtered VCF used for LD pruning and PCA after removal of the related sample: `filtered_no_nx1575.vcf.gz`
- LD pruning files: `pruned.prune.in`, `pruned.prune.out`, `pruned.log`
- PCA input files: `pruned.bed`, `pruned.bim`, `pruned.fam`, `pruned.nosex`
- PCA result files: `pruned.eigenval`, `pruned.eigenvec`
- VCF normalization statistics table (`Sample`, `Total`, `Split`, `Joined`, `Realigned`, `Skipped`): `normalization_stats.tsv`
- Updated sample-group table with newly added sample IDs (group = `new`): `jews_and_new_groups.csv`
- Population assignment table for new samples based on PCA ellipses: `samples_in_ellipses.csv`
- PCA plots: `PCA_PC1_PC2.pdf`, `PCA_PC1_PC3_with_ellipses.pdf`
 <img width="487" height="580" alt="520214853-f5e3133c-918b-4d49-b575-a0aeeddb8022" src="https://github.com/user-attachments/assets/2ff5068a-c832-4932-9610-ea88a54123b4" />

