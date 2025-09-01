# Data Dictionary (Urate–eGFR MR, Japanese cohorts)

This document describes the public/sharable, **de-identified** derived datasets that accompany the manuscript.

## Files

### 1) `MR_instruments_34.tsv`
34 genome-wide significant SNP instruments for SUA used as MR exposure instruments.

| Column | Type | Unit | Description |
|---|---|---|---|
| rsID | string | – | dbSNP ID |
| chr | integer | – | Chromosome (1–22) |
| pos | integer | bp | Base-pair position (GRCh37/38; specify in `provenance.yml`) |
| EA | string | allele | Effect allele for SUA (exposure) |
| OA | string | allele | Other allele |
| EAF | float | proportion | Effect allele frequency in exposure GWAS (J-MICC) |
| beta_exposure | float | SD change | Per-allele effect on SUA (z-scored) |
| se_exposure | float | SD change | Standard error of `beta_exposure` |
| p_exposure | float | – | P-value in exposure GWAS |
| N_exposure | integer | persons | Sample size for exposure effect |
| source_exposure | string | – | `"J-MICC Study"` |

---

### 2) `harmonized_JMICC_UA_ToMMo_eGFR.tsv`
Per-variant **harmonized** table used as MR input (SNP-level).

| Column | Type | Unit | Description |
|---|---|---|---|
| rsID | string | – | dbSNP ID |
| EA | string | allele | Harmonized effect allele (common to exposure/outcome) |
| OA | string | allele | Harmonized other allele |
| EAF | float | proportion | Effect allele frequency (prefer exposure source) |
| beta_exposure | float | SD change | SUA effect (J-MICC; z-scored) |
| se_exposure | float | SD change | SE of `beta_exposure` |
| beta_outcome | float | SD change | eGFR effect (ToMMo; rank-based inverse-normalized residuals) |
| se_outcome | float | SD change | SE of `beta_outcome` |
| N_exposure | integer | persons | Exposure sample size |
| N_outcome | integer | persons | Outcome sample size |
| harmonization_flags | string | – | Flags applied (e.g., `None`, `Allele flipped`, `Palindromic-removed`) |

---

### 3) `rs121907892_effects_overall.tsv`
Effect estimates for rs121907892 across cohorts/sex strata **(model-based)**.

| Column | Type | Unit | Description |
|---|---|---|---|
| cohort | string | – | `ToMMo` / `J-MICC` |
| sex | string | – | `overall` / `male` / `female` |
| phenotype | string | – | e.g., `eGFR` (rank-based inverse-normalized) |
| model_covariates | string | – | e.g., `age + sex + PC1–PC10`  |
| EA | string | allele | Effect allele |
| OA | string | allele | Other allele |
| EAF | float | proportion | Effect allele frequency |
| beta | float | per-allele effect | Scale per `phenotype` |
| se | float | – | Standard error |
| p | float | – | P-value |
| N | integer | persons | Sample size for the model |
| scale_info | string | – | e.g., `eGFR rank-based inverse-normalized`; SUA `z-scored` |

---

### 4) `rs121907892_demographics_by_genotype.tsv`
Group-level demographics by genotype (de-identified; **AA suppressed if n<10**).

| Column | Type | Unit | Description |
|---|---|---|---|
| cohort | string | – | `J-MICC` / `ToMMo` |
| sex | string | – | `overall` / `male` / `female` |
| genotype | string | – | `GG`, `GA` (AA omitted if n<10) |
| N | integer | persons | Group size |
| age_mean / age_sd | float | years | Baseline age |
| eGFR_mean / eGFR_sd | float | mL/min/1.73m² | Raw eGFR |
| SUA_mean / SUA_sd | float | mg/dL | Serum urate |
| BMI_mean / BMI_sd | float | kg/m² | BMI |
| SBP_mean / DBP_mean | float | mmHg | Blood pressure |
| TG_mean / TChol_mean | float | mg/dL | Lipids |
| hypertension_prev | float | % | Prevalence |

---

### 5) `loo_results.tsv`
Leave-one-out MR results for each instrument.

| Column | Type | Unit | Description |
|---|---|---|---|
| method | string | – | `IVW` / `Weighted median` / `Weighted mode` |
| rsID_left_out | string | – | SNP removed in this run |
| beta | float | SD change | MR estimate |
| se | float | – | SE |
| p | float | – | P-value |