plink2 \
  --bfile [GENOTYPE_INPUT_FILE] \
  --geno 0.01 \
  --hwe 0.05 \
  --keep-allele-order \
  --maf 0.05 \
  --make-bed \
  --out [QC_OUTPUT_FILE]

plink \
  --bfile [QC_OUTPUT_FILE] \
  --indep-pairwise 1500 150 0.03 \
  --out [PRUNE_OUTPUT_FILE]

plink \
  --bfile [QC_OUTPUT_FILE] \
  --extract [PRUNE_OUTPUT_FILE].prune.in \
  --keep-allele-order \
  --make-bed \
  --out [PRUNED_BED_OUTPUT_FILE]

plink \
  --bfile [PRUNED_BED_OUTPUT_FILE] \
  --out [PCA_OUTPUT_FILE] \
  --pca
