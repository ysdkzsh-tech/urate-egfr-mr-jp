// Set working directory
global wd1 "[WORKING_DIRECTORY]/"

// Start log
log using "${wd1}analysis_log.smcl", append

********************************************************************************
******	1. Data cleaning
********************************************************************************
** 
**	Import baseline unique IDs (prepared separately)
import excel "${wd1}baseline_file.xlsx", firstrow sheet("Sheet1") clear
rename 匿名化ID uid
drop if uid==""
sort uid
save temp_uid, replace

import excel "${wd1}followup_file.xlsx", firstrow sheet("Sheet1") clear
rename 匿名化ID uid
drop if uid==""
sort uid
merge 1:1 uid using temp_uid 

drop 死亡フラグ 死亡年月日 ICD10

save temp_uid, replace

// Create egfr1 dataset
use temp_uid, clear
drop _merge
rename (ｸﾚｱﾁﾆﾝ_修正後 尿酸 身長 体重 尿蛋白 ﾍﾓｸﾞﾛﾋﾞﾝ HbA1cJDS従来 開始時年齢 性別 収縮期血圧 拡張期血圧 常用薬血圧 常用薬血糖値 糖尿病_修正 高血圧 高脂血症_修正 狭心症心筋梗塞_修正 脳卒中_修正 喫煙有無 検体採取日_修正後 中性脂肪 総ｺﾚｽﾃﾛｰﾙ 尿素窒素) ///
       (cre1 ua1 height1 weight1 up1 hgb1 HbA1c1 age1 sex sbp1 dbp1 m_ht m_dm q_dm q_ht q_dlp q_mi q_stroke q_smoke entrydt TG TChol BUN)
rename (ｸﾚｱﾁﾆﾝ酵素法_2次 同意日の年齢_2次 受診年月日_2次) (cre2 age2 followdt)
drop if cre1==99999 | ua1==99999
format followdt %td

keep uid cre1 cre2 ua1 height1 weight1 up1 hgb1 HbA1c1 age1 age2 sex sbp1 dbp1 m_ht m_dm q_dm q_ht q_dlp q_mi q_stroke q_smoke q_beer entrydt followdt TG TChol BUN
gen bmi = weight1/(height1/100)^2
gen egfr1= 194*cre1^(-1.094)*age1^(-0.287)
replace egfr1 = egfr1*0.739 if sex==2
gen ckd = 1 if egfr1>=90
replace ckd=2 if egfr1>=60 & egfr1<90
replace ckd=3 if egfr1<60 & egfr1>=30
replace ckd=4 if egfr1<30 & egfr1>=15
replace ckd=5 if egfr1<15
gen c_dm = 0
replace c_dm = 1 if HbA1c1>7.0 & HbA1c1!=99999

gen cut_beer = 0
replace cut_beer = 1 if q_beer>=350

sort uid
save user_id_entrydt_egfr1, replace	

// Prepare GWAS trait dataset for eGFR1
use temp_uid, clear
drop _merge
rename (ｸﾚｱﾁﾆﾝ_修正後 尿酸 身長 体重 尿蛋白 ﾍﾓｸﾞﾛﾋﾞﾝ HbA1cJDS従来 開始時年齢 性別 収縮期血圧 拡張期血圧 糖尿病_修正 高血圧 高脂血症_修正 狭心症心筋梗塞_修正 脳卒中_修正 喫煙有無 検体採取日_修正後) ///
       (cre1 ua1 height1 weight1 up1 hgb1 HbA1c1 age1 sex sbp1 dbp1 q_dm q_ht q_dlp q_mi q_stroke q_smoke entrydt)
drop if cre1==99999 | ua1==99999 

gen egfr1= 194*cre1^(-1.094)*age1^(-0.287)
replace egfr1 = egfr1*0.739 if sex==2
rename uid IID
gen FID = IID, after(IID)
sort FID
save jmicc_egfr1, replace

// Merge with PCA eigenvectors
import delimited "${wd1}JMICC_190729.QC.prune.pca.eigenvec", delimiter(space) encoding(UTF-8) clear
rename v1 FID
sort FID
merge 1:1 FID using jmicc_egfr1
keep if _merge==3
drop _merge
save jmicc_egfr1, replace

use jmicc_egfr1, clear

// Step 1: Regress trait on age, sex, and PCs
regress egfr1 age1 sex v3 v4 v5 v6 v7 v8 v9 v10 v11 v12

// Step 2: Get the residuals from regression
predict residuals, residuals
egen rank = rank(residuals), field

// Step 3: Inverse normal transformation
gen norm_residuals = invnorm((_N + 1 - rank) / (_N + 1))
rename norm_residuals egfr1_rd
drop residuals

// For SUA (ua1) as well
regress ua1 age1 sex v3 v4 v5 v6 v7 v8 v9 v10 v11 v12
predict residuals, residuals
egen mean_resid = mean(residuals)
egen sd_resid = sd(residuals)
gen zscore_resid = (residuals - mean_resid) / sd_resid
rename zscore_resid ua1_z

keep IID egfr1 ua1 egfr1_rd ua1_z 
export delimited using "${wd1}pheno2.tsv", delimiter(tab) replace	

// Create covariate file for egfr1
use user_id_entrydt_egfr1, clear

// Replace missing values with mean (example shown for height1)
summarize height1
local tmp = r(mean)-99999*4/10795
replace height1 = `tmp' if height1 == 99999

summarize weight1
local tmp = r(mean)-99999*4/10795
replace weight1 = `tmp' if weight1 == 99999

summarize TChol
local tmp = r(mean) - 99999*184/10795
replace TChol = `tmp' if TChol == 99999

summarize sbp1
local tmp = r(mean)-99999*29/10795
replace sbp1 = `tmp' if sbp1 == 99999

summarize dbp1
local tmp = r(mean)-99999*30/10795
replace dbp1 = `tmp' if dbp1 == 99999

drop bmi
gen bmi = weight1/(height1/100)^2
replace q_smoke=3 if q_smoke==0 | q_smoke==9

rename uid IID
merge 1:1 IID using covar_egfr1_pca
keep if _merge == 3
drop _merge
rename IID uid

gen ht = 0
replace ht = 1 if sbp > 140 | dbp > 90 | m_ht==1

keep uid ua1 cre1 egfr1 TG TChol age1 sex bmi ht m_ht q_ht q_smoke q_beer C1 C2 C3 C4 C5 C6 C7 C8 C9 C10
export delimited using "${wd1}covariates_egfr1.csv",replace

// GWAS results for SUA and eGFR - merge 22 chromosomes
clear
import delimited "${wd1}ua1_mr/chr1.imputed.IDconverted.filtered07292019.dose.recode.QC.linear.cov.ua1_z.glm.linear", varnames(1) clear

forvalues n=1/22{
    cap: import delimited "${wd1}ua1_mr/chr`n'.imputed.IDconverted.filtered07292019.dose.recode.QC.linear.cov.ua1_z.glm.linear", varnames(1) clear
	cap: save _temp`n', replace
	clear
}
use _temp1, clear
forvalues n=2/22{
    cap: append using _temp`n', force
}
sort id
replace beta = beta * (-1) if a1!=alt
save ua_beta, replace

// GWAS results for eGFR_rd
clear
import delimited "${wd1}egfr1_mr/chr1.imputed.IDconverted.filtered07292019.dose.recode.QC.linear.cov.egfr1_rd.glm.linear", varnames(1) clear

forvalues n=1/22{
    cap: import delimited "${wd1}egfr1_mr/chr`n'.imputed.IDconverted.filtered07292019.dose.recode.QC.linear.cov.egfr1_rd.glm.linear", varnames(1) clear
    cap: save _temp`n', replace
    clear
}
use _temp1, clear
forvalues n=2/22{
    cap: append using _temp`n', force
}
sort id
replace beta = beta * (-1) if a1!=alt
save egfr_beta_jmicc, replace

// Prepare SNP set for publication
use ua_beta, clear
merge 1:1 id using bbjsnps_beta
keep if _merge==3
drop _merge
keep rsID CHR BP REF ALT ALT_Freq ClosestGene
export delimited using "${wd1}SNPLIST.tsv", delimiter(tab) replace

// Prepare MR variables for ToMMo eGFR
import delimited "${wd1}TGA7_egfr_summary_statistics.tsv", varnames(1) clear
keep if imputation_info_score > 0.8
gen str id = substr(variant_id, 1, strpos(variant_id, "_")-1), after(variant_id)
sort id
duplicates drop id, force
merge 1:1 id using bbjsnps_beta
keep if _merge==3
drop _merge

save ToMMo_beta_egfr, replace

// Create harmonized file for MR
use ToMMo_beta_egfr, clear
sort chromosome base_pair_location
replace LeadVariant = rsID
keep id LeadVariant beta reference_allele alternative_allele standard_error
rename beta betaToMMo
merge 1:1 id using ua_beta
keep if _merge == 3
drop _merge 
rename (ref alt beta se) (REF ALT Beta SE)
rename (reference_allele alternative_allele betaToMMo standard_error) (ref alt beta se)
gen refchk=1 if REF==ref, after(REF)
gen altchk=1 if ALT==alt, after(ALT)

gen _REF = REF, after(REF)
gen _ALT = ALT, after(ALT)
replace _REF = ALT if Beta<0
replace _ALT = REF if Beta<0
replace beta=beta * (-1) if Beta<0
replace Beta=Beta * (-1) if Beta<0

replace REF = _REF
replace ALT = _ALT

drop _ALT _REF

keep LeadVariant Beta SE beta se
drop if LeadVariant == "rs671" 
export delimited using "${wd1}beta_JMICC_UA_ToMMo_egfr.tsv", delimiter(tab) replace

log close