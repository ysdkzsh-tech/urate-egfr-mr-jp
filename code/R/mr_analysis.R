library(MendelianRandomization)
library(ggplot2)
library(MRPRESSO)
library(dplyr)
library(ggpubr)

# Set working directory (Recommend using relative path or project root)
# setwd("./")  # <- Uncomment and adjust if necessary

perform_mr_analysis <- function(input_file, output_csv, output_plot, output_ivw_csv, output_format = "jpg") {
  library(MendelianRandomization)
  library(ggplot2)
  library(dplyr)
  library(tools)
  
  # Color settings for regression line legend
  color_map <- c(
    "Inverse variance weighted" = "#1b3c41",
    "Weighted median" = "#e97a20",
    "Weighted mode" = "#2992c9",
    "MR Egger" = "#a1271c"
  )
  
  # Load data
  mydata <- read.delim(input_file)
  
  # Create MRInput object
  MRInputObject <- mr_input(
    bx = mydata$Beta, 
    bxse = mydata$SE, 
    by = mydata$beta, 
    byse = mydata$se, 
    snps = mydata$LeadVariant
  )
  
  num_snps <- length(MRInputObject@snps)
  mr_results <- data.frame()
  models <- list()
  
  # Run each MR method
  if (num_snps >= 1) {
    try({
      ivw_random <- mr_ivw(MRInputObject, model = "random")
      models[["Inverse variance weighted"]] <- list(Estimate = ivw_random@Estimate,
                                                    SE = ivw_random@StdError,
                                                    P = ivw_random@Pvalue)
    }, silent = TRUE)
  }
  
  if (num_snps >= 3) {
    try({
      median_wt <- mr_median(MRInputObject, weighting = "weighted")
      models[["Weighted median"]] <- list(Estimate = median_wt@Estimate,
                                          SE = median_wt@StdError,
                                          P = median_wt@Pvalue)
    }, silent = TRUE)
    try({
      mbe_wt <- mr_mbe(MRInputObject, weighting = "weighted")
      models[["Weighted mode"]] <- list(Estimate = mbe_wt@Estimate,
                                        SE = mbe_wt@StdError,
                                        P = mbe_wt@Pvalue)
    }, silent = TRUE)
    try({
      egger <- mr_egger(MRInputObject)
      models[["MR Egger"]] <- list(Estimate = egger@Estimate,
                                   SE = egger@StdError.Est,
                                   P = egger@Pvalue.Est,
                                   Intercept = egger@Intercept,
                                   InterceptSE = egger@StdError.Int,
                                   InterceptP = egger@Pvalue.Int)
    }, silent = TRUE)
  }
  
  for (model_name in names(models)) {
    result <- models[[model_name]]
    ci_lower <- result$Estimate - 1.96 * result$SE
    ci_upper <- result$Estimate + 1.96 * result$SE
    mr_results <- rbind(mr_results, data.frame(
      Model = model_name,
      Beta = result$Estimate,
      SE = result$SE,
      P_value = result$P,
      CI_lower = ci_lower,
      CI_upper = ci_upper
    ))
  }
  
  write.csv(mr_results, file.path("results", output_csv), row.names = FALSE)
  
  if ("MR Egger" %in% names(models)) {
    egger_row <- models[["MR Egger"]]
    egger_intercept_df <- data.frame(
      Intercept = egger_row$Intercept,
      SE = egger_row$InterceptSE,
      P_value = egger_row$InterceptP,
      CI_lower = egger_row$Intercept - 1.96 * egger_row$InterceptSE,
      CI_upper = egger_row$Intercept + 1.96 * egger_row$InterceptSE
    )
    write.csv(egger_intercept_df, file.path("results", "egger_intercept_results.csv"), row.names = FALSE)
  }
  
  # Summarize regression line information for plot
  line_df <- data.frame(
    Method = character(),
    Intercept = numeric(),
    Slope = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (model_name in names(models)) {
    est <- models[[model_name]]$Estimate
    int <- ifelse(is.null(models[[model_name]]$Intercept), 0, models[[model_name]]$Intercept)
    line_df <- rbind(line_df, data.frame(Method = model_name, Intercept = int, Slope = est))
  }
  
  used_colors <- color_map[line_df$Method]
  used_linetypes <- rep("solid", length(line_df$Method))
  names(used_linetypes) <- line_df$Method
  
  # Create scatter plot with regression lines
  plot <- ggplot(mydata, aes(x = Beta, y = beta)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray70") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray70") +
    geom_errorbar(aes(ymin = beta - 1.96 * se, ymax = beta + 1.96 * se), width = 0, color = "black") +
    geom_errorbarh(aes(xmin = Beta - 1.96 * SE, xmax = Beta + 1.96 * SE), height = 0, color = "black") +
    geom_point(color = "black", size = 2) +
    geom_abline(data = line_df,
                aes(intercept = Intercept, slope = Slope, color = Method, linetype = Method),
                size = 1.2, show.legend = TRUE) +
    scale_color_manual(values = used_colors) +
    scale_linetype_manual(values = used_linetypes) +
    labs(x = "SNP effect on exposure: serum urate level", y = "SNP effect on outcome: eGFR") +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = 11),
      plot.margin = margin(10, 10, 10, 10)
    )
  
  # Save the plot
  plot_path <- file.path("results", paste0(tools::file_path_sans_ext(output_plot), ".", output_format))
  ggsave(plot_path, plot = plot, width = 8, height = 6, dpi = 300, device = output_format)
  
  return(list(MRResults = mr_results, Models = models))
}

# Specify input and output files (use relative path)
analysis_configs <- list(
  list(
    input_file = "data/beta_exposure_outcome.tsv",
    output_csv = "MR_results.csv",
    output_plot = "MR_scatter_plot.png",
    output_ivw_csv = "MR_ivw_results.csv"
  ),
  list(
    input_file = "data/beta_exposure_outcome_wo_outliers.tsv",
    output_csv = "MR_results_wo_outliers.csv",
    output_plot = "MR_scatter_plot_wo_outliers.png",
    output_ivw_csv = "MR_ivw_results_wo_outliers.csv"
  )
)

# Run MR analysis for each config
results <- lapply(analysis_configs, function(config) {
  perform_mr_analysis(
    config$input_file,
    config$output_csv,
    config$output_plot,
    config$output_ivw_csv
  )
})

print("Analysis completed, results and heterogeneity assessments saved.")

# (1) Common function for Leave-One-Out (LOO) analysis
perform_loo_analysis <- function(input_file, method, output_prefix) {
  # Load data
  mydata <- read.delim(input_file)
  
  # Create MRInput object
  MRInputObject <- mr_input(bx = mydata$Beta, bxse = mydata$SE, by = mydata$beta, byse = mydata$se, snps = mydata$LeadVariant)
  
  # Data frame to save results
  loo_results <- data.frame()
  
  for (i in seq_along(MRInputObject@snps)) {
    # Remove one SNP at a time
    bx_subset <- MRInputObject@betaX[-i]
    bxse_subset <- MRInputObject@betaXse[-i]
    by_subset <- MRInputObject@betaY[-i]
    byse_subset <- MRInputObject@betaYse[-i]
    snps_subset <- MRInputObject@snps[-i]
    
    # Reconstruct MRInput for the subset
    MRInput_subset <- mr_input(
      bx = bx_subset,
      bxse = bxse_subset,
      by = by_subset,
      byse = byse_subset,
      snps = snps_subset
    )
    
    # Run analysis using the specified method
    MRResult <- switch(
      method,
      "Penalized IVW" = mr_ivw(MRInput_subset, penalized = TRUE),
      "Robust IVW" = mr_ivw(MRInput_subset, robust = TRUE),
      "Penalized robust IVW" = mr_ivw(MRInput_subset, robust = TRUE, penalized = TRUE),
      "Weighted Median" = mr_median(MRInput_subset, weighting = "weighted"),
      "Weighted Mode" = mr_mbe(MRInput_subset, weighting = "weighted"),
      "IVW" = mr_ivw(MRInput_subset, model = "random"),
      "MREgger" = mr_egger(MRInput_subset)
    )
    
    # Save results
    loo_results <- rbind(
      loo_results,
      data.frame(
        Removed_SNP = MRInputObject@snps[i],
        Estimate = MRResult@Estimate,
        StdError = MRResult@StdError,
        PValue = MRResult@Pvalue
      )
    )
  }
  
  # Save results to CSV
  write.csv(loo_results, file.path("results", paste0(output_prefix, "_", method, "_LeaveOneOut_Results.csv")), row.names = FALSE)
  
  # (2) Generate LOO plot
  loo_results$Removed_SNP <- factor(loo_results$Removed_SNP, levels = loo_results$Removed_SNP)
  
  loo_plot <- ggplot(loo_results, aes(x = Removed_SNP, y = Estimate)) +
    geom_point() +
    geom_errorbar(aes(
      ymin = Estimate - 1.96 * StdError,
      ymax = Estimate + 1.96 * StdError
    ), width = 0.2) +
    labs(
      title = paste0("Leave-One-Out Analysis - ", method),
      x = "Removed SNP",
      y = "Causal Estimate"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Save LOO plot as JPG
  ggsave(
    file.path("results", paste0(output_prefix, "_", method, "_LeaveOneOut_Plot.jpg")),
    loo_plot, 
    width = 8, 
    height = 6, 
    dpi = 300,
    device = "jpeg"
  )
}

# (2) Specify input data and output prefix
input_file = "data/beta_exposure_outcome.tsv"
output_prefix <- "MR_results"

# (3) Run LOO analysis for each method
perform_loo_analysis(input_file, "Penalized IVW", output_prefix)
perform_loo_analysis(input_file, "Robust IVW", output_prefix)
perform_loo_analysis(input_file, "Penalized robust IVW", output_prefix)
perform_loo_analysis(input_file, "Weighted Median", output_prefix)
perform_loo_analysis(input_file, "Weighted Mode", output_prefix)
perform_loo_analysis(input_file, "IVW", output_prefix)
perform_loo_analysis(input_file, "MREgger", output_prefix)

# Run MR-PRESSO

# Load data
mydata <- read.delim(input_file)

# Run MR-PRESSO
press_result <- mr_presso(
  BetaOutcome = "beta",     
  BetaExposure = "Beta",    
  SdOutcome = "se",         
  SdExposure = "SE",        
  OUTLIERtest = TRUE,
  DISTORTIONtest = TRUE,
  NbDistribution = 1000,
  SignifThreshold = 0.05,
  data = mydata
)

# Print result
print(press_result)

# 1. Get Outlier Test data
outlier_test <- press_result$`MR-PRESSO results`$`Outlier Test`

# 2. Handle p-values as string
raw_pvals <- outlier_test$Pvalue

# 3. Convert to numeric ("<0.034" -> 0.034)
pvals <- suppressWarnings(as.numeric(raw_pvals))
pvals_na_idx <- which(is.na(pvals) & grepl("^<", raw_pvals))
pvals[pvals_na_idx] <- as.numeric(sub("^<", "", raw_pvals[pvals_na_idx]))

# 4. Identify outliers (e.g., p < 0.05)
outlier_indices <- which(pvals < 0.05)

# 5. Output
if (length(outlier_indices) > 0) {
  outlier_variants <- mydata$LeadVariant[outlier_indices]
  outlier_df <- data.frame(
    Outlier_Index = outlier_indices,
    LeadVariant = outlier_variants,
    RSSobs = outlier_test$RSSobs[outlier_indices],
    Pvalue = raw_pvals[outlier_indices]
  )
  write.csv(outlier_df, file.path("results", "mrpresso_outliers.csv"), row.names = FALSE)
  message("Outlier SNPs saved to: results/mrpresso_outliers.csv")
} else {
  message("No outliers detected based on p-value < 0.05")
}
# Save other statistics (Global test, etc.)
global_test <- press_result$`MR-PRESSO results`$`Global Test`
write.csv(as.data.frame(global_test), file.path("results", "mrpresso_global_test.csv"), row.names = FALSE)

# Optional: distortion test if exists
if (!is.null(press_result$`MR-PRESSO results`$`Distortion Test`)) {
  distortion_test <- press_result$`MR-PRESSO results`$`Distortion Test`
  write.csv(as.data.frame(distortion_test), file.path("results", "mrpresso_distortion_test.csv"), row.names = FALSE)
}

# (Optional) Remove specific SNP ("11:64361219:G:A")
filtered_data <- subset(mydata, LeadVariant != "11:64361219:G:A")

# Run MR-PRESSO for filtered data
mrpresso_result <- mr_presso(BetaOutcome = "beta", 
                             BetaExposure = "Beta", 
                             SdOutcome = "se", 
                             SdExposure = "SE", 
                             OUTLIERtest = TRUE, 
                             DISTORTIONtest = TRUE, 
                             data = filtered_data, 
                             NbDistribution = 1000)

# Print result
print(mrpresso_result)
