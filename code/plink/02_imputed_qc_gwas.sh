for chr in $(seq 1 22); do
  bcftools filter -e 'R2 < 0.8' [IMPUTED_VCF_INPUT_DIR]/chr${chr}.imputed.vcf.gz -O z -o [QC_VCF_OUTPUT_DIR]/chr${chr}.imputed.QC.vcf.gz
done

for chr in $(seq 1 22); do
  plink2 --vcf [QC_VCF_OUTPUT_DIR]/chr${chr}.imputed.QC.vcf.gz --make-pgen --keep-allele-order --out [PGEN_OUTPUT_DIR]/chr${chr}.imputed.QC
done

for chr in $(seq 1 22); do
  plink2 --pfile [PGEN_OUTPUT_DIR]/chr${chr}.imputed.QC --make-pgen --extract [TARGET_SNP_LIST] --out [TARGET_SNP_OUTPUT_DIR]/chr${chr}.imputed.QC.targetSNP
done

for chr in $(seq 1 22); do
  plink2 \
    --pfile [TARGET_SNP_OUTPUT_DIR]/chr${chr}.imputed.QC.targetSNP \
    --pheno [PHENOTYPE_FILE] \
    --pheno-name [STANDARDIZED_PHENO_NAME] \
    --covar [COVARIATE_FILE] \
    --covar-name C1 C2 C3 C4 C5 C6 C7 C8 C9 C10 age1 sex \
    --covar-variance-standardize \
    --glm hide-covar \
    --ci 0.95 \
    --out [GWAS_OUTPUT_DIR]/chr${chr}.imputed.QC.gwas
done
