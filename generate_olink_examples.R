# Generate Example Plots and Documentation for Olink Functions
# This script creates example visualizations and README files

# Load required libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyr)
})

# Source the functions
source("R/olink-covariate-adjustment.R")
source("R/olink-regression-analysis.R")

# Create output directory for figures
output_dir <- "olink_examples"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

cat("=== Generating Olink Analysis Examples ===\n\n")

# ==============================================================================
# Example 1: Covariate Adjustment
# ==============================================================================
cat("1. Generating covariate adjustment examples...\n")

set.seed(12345)

# Create synthetic Olink data with batch effects
n_proteins <- 30
n_samples <- 60

count_matrix_cov <- matrix(
  rnorm(n_proteins * n_samples, mean = 100, sd = 20),
  nrow = n_proteins,
  ncol = n_samples
)
rownames(count_matrix_cov) <- paste0("Protein_", 1:n_proteins)
colnames(count_matrix_cov) <- paste0("Sample_", 1:n_samples)

# Create metadata with covariates
metadata_cov <- data.frame(
  sample_id = paste0("Sample_", 1:n_samples),
  age = rnorm(n_samples, mean = 50, sd = 15),
  sex = sample(c("M", "F"), n_samples, replace = TRUE),
  batch = sample(c("Batch1", "Batch2", "Batch3"), n_samples, replace = TRUE),
  stringsAsFactors = FALSE
)

# Add batch effects
batch_effects <- ifelse(
  metadata_cov$batch == "Batch1", 20,
  ifelse(metadata_cov$batch == "Batch2", -15, 0)
)
for (i in 1:n_proteins) {
  count_matrix_cov[i, ] <- count_matrix_cov[i, ] + batch_effects
}

# Add age effects to some proteins
for (i in 1:10) {
  age_effect <- scale(metadata_cov$age) * 10
  count_matrix_cov[i, ] <- count_matrix_cov[i, ] + age_effect
}

# Perform covariate adjustment
result_cov <- adjust_olink_covariates(
  count_matrix = count_matrix_cov,
  metadata = metadata_cov,
  covariates = c("age", "sex", "batch"),
  sample_id_col = "sample_id",
  method = "adjusted",
  verbose = TRUE
)

cat("  - Adjustment complete\n")

# Create batch effect visualization
plot_data_batch <- data.frame()
for (protein in paste0("Protein_", 1:4)) {
  temp_df <- data.frame(
    protein = protein,
    batch = metadata_cov$batch,
    original = count_matrix_cov[protein, ],
    adjusted = result_cov$adjusted_matrix[protein, ]
  )
  plot_data_batch <- rbind(plot_data_batch, temp_df)
}

plot_data_batch_long <- tidyr::pivot_longer(
  plot_data_batch,
  cols = c("original", "adjusted"),
  names_to = "type",
  values_to = "expression"
)

p2 <- ggplot(plot_data_batch_long, aes(x = batch, y = expression, fill = type)) +
  geom_boxplot() +
  facet_wrap(~ protein, scales = "free_y", ncol = 2) +
  labs(
    title = "Batch Effect Removal",
    subtitle = "Protein expression before and after covariate adjustment",
    x = "Batch",
    y = "Expression Level",
    fill = "Data Type"
  ) +
  theme_bw() +
  scale_fill_manual(
    values = c("original" = "#E74C3C", "adjusted" = "#3498DB"),
    labels = c("Original", "Adjusted")
  ) +
  theme(legend.position = "top")

ggsave(
  file.path(output_dir, "batch_effect_removal.png"),
  plot = p2, width = 10, height = 8, dpi = 300
)
cat("  - Saved: batch_effect_removal.png\n")

# Create R-squared distribution plot
r_squared_vals <- sapply(result_cov$model_info, function(x) x$r_squared)
p3 <- ggplot(data.frame(r_squared = r_squared_vals), aes(x = r_squared)) +
  geom_histogram(fill = "#3498DB", color = "white", bins = 20) +
  labs(
    title = "Model R² Distribution",
    subtitle = "Variance explained by covariates (age, sex, batch)",
    x = "R² Value",
    y = "Number of Proteins"
  ) +
  theme_bw() +
  geom_vline(
    xintercept = median(r_squared_vals),
    linetype = "dashed",
    color = "#E74C3C", size = 1
  ) +
  annotate(
    "text",
    x = median(r_squared_vals) + 0.1, y = Inf, vjust = 2,
    label = paste("Median =", round(median(r_squared_vals), 3)),
    color = "#E74C3C"
  )

ggsave(
  file.path(output_dir, "covariate_r_squared.png"),
  plot = p3, width = 8, height = 6, dpi = 300
)
cat("  - Saved: covariate_r_squared.png\n\n")

# ==============================================================================
# Example 2: Regression Analysis - Binary Outcome
# ==============================================================================
cat("2. Generating binary outcome regression examples...\n")

set.seed(54321)

n_proteins <- 50
n_samples <- 80

count_matrix_binary <- matrix(
  rnorm(n_proteins * n_samples, mean = 100, sd = 25),
  nrow = n_proteins,
  ncol = n_samples
)
rownames(count_matrix_binary) <- paste0("Protein_", 1:n_proteins)
colnames(count_matrix_binary) <- paste0("Sample_", 1:n_samples)

# Create metadata with binary outcome
metadata_binary <- data.frame(
  sample_id = paste0("Sample_", 1:n_samples),
  disease_status = sample(c(0, 1), n_samples, replace = TRUE),
  age = rnorm(n_samples, mean = 55, sd = 12),
  sex = sample(c("M", "F"), n_samples, replace = TRUE),
  bmi = rnorm(n_samples, mean = 26, sd = 4),
  stringsAsFactors = FALSE
)

# Add disease association to some proteins
disease_associated <- c(1, 2, 3, 5, 8, 13, 21)
for (i in disease_associated) {
  disease_effect <- ifelse(metadata_binary$disease_status == 1, 30, -10)
  count_matrix_binary[i, ] <- count_matrix_binary[i, ] + disease_effect +
    rnorm(n_samples, 0, 15)
}

# Run regression analysis
result_binary <- olink_regression_analysis(
  count_matrix = count_matrix_binary,
  metadata = metadata_binary,
  outcome = "disease_status",
  outcome_type = "binary",
  covariates = c("age", "sex", "bmi"),
  sample_id_col = "sample_id",
  p_adjust_method = "BH",
  verbose = TRUE
)

cat("  - Analysis complete\n")

# Create forest plot for univariate results
p4 <- plot_forest_results(
  results = result_binary,
  analysis_type = "univariate",
  top_n = 20,
  sort_by = "p_value",
  show_pval = TRUE,
  title = "Univariate Analysis: Disease Association (Top 20 Proteins)"
)

ggsave(
  file.path(output_dir, "forest_plot_univariate_binary.png"),
  plot = p4, width = 12, height = 10, dpi = 300
)
cat("  - Saved: forest_plot_univariate_binary.png\n")

# Create forest plot for multivariate results
p5 <- plot_forest_results(
  results = result_binary,
  analysis_type = "multivariate",
  top_n = 20,
  sort_by = "p_value",
  show_pval = TRUE,
  title = "Multivariate Analysis: Disease Association (Adjusted)"
)

ggsave(
  file.path(output_dir, "forest_plot_multivariate_binary.png"),
  plot = p5, width = 12, height = 10, dpi = 300
)
cat("  - Saved: forest_plot_multivariate_binary.png\n")

# Create volcano plot
p6 <- ggplot(
  result_binary$univariate_results,
  aes(x = log_OR, y = -log10(p_value))
) +
  geom_point(aes(color = p_adjusted < 0.05), alpha = 0.6, size = 3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  scale_color_manual(
    values = c("TRUE" = "#E74C3C", "FALSE" = "#95A5A6"),
    labels = c("TRUE" = "FDR < 0.05", "FALSE" = "Not Significant"),
    name = ""
  ) +
  labs(
    title = "Volcano Plot: Protein Association with Disease",
    subtitle = "Univariate analysis results",
    x = "log(Odds Ratio)",
    y = "-log10(p-value)"
  ) +
  theme_bw() +
  theme(legend.position = "top")

ggsave(
  file.path(output_dir, "volcano_plot_binary.png"),
  plot = p6, width = 10, height = 8, dpi = 300
)
cat("  - Saved: volcano_plot_binary.png\n")

# Comparison: Univariate vs Multivariate
comparison_data <- merge(
  result_binary$univariate_results[, c("protein", "OR", "p_value")],
  result_binary$multivariate_results[, c("protein", "OR", "p_value")],
  by = "protein",
  suffixes = c("_uni", "_multi")
)

p7 <- ggplot(
  comparison_data,
  aes(x = log2(OR_uni), y = log2(OR_multi))
) +
  geom_point(alpha = 0.6, size = 3, color = "#3498DB") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed",
              color = "#E74C3C") +
  labs(
    title = "Univariate vs Multivariate Odds Ratios",
    subtitle = "Effect of covariate adjustment on protein associations",
    x = "log2(OR) - Univariate",
    y = "log2(OR) - Multivariate (Adjusted)"
  ) +
  theme_bw() +
  annotate(
    "text", x = Inf, y = -Inf, hjust = 1.1, vjust = -1,
    label = "Red line: y = x\n(no change)",
    color = "#E74C3C", size = 3.5
  )

ggsave(
  file.path(output_dir, "univariate_vs_multivariate.png"),
  plot = p7, width = 8, height = 8, dpi = 300
)
cat("  - Saved: univariate_vs_multivariate.png\n\n")

# ==============================================================================
# Example 3: Regression Analysis - Continuous Outcome
# ==============================================================================
cat("3. Generating continuous outcome regression examples...\n")

set.seed(99999)

count_matrix_cont <- matrix(
  rnorm(n_proteins * n_samples, mean = 80, sd = 20),
  nrow = n_proteins,
  ncol = n_samples
)
rownames(count_matrix_cont) <- paste0("Protein_", 1:n_proteins)
colnames(count_matrix_cont) <- paste0("Sample_", 1:n_samples)

# Create metadata with continuous outcome
metadata_cont <- data.frame(
  sample_id = paste0("Sample_", 1:n_samples),
  severity_score = rnorm(n_samples, mean = 50, sd = 20),
  age = rnorm(n_samples, mean = 55, sd = 12),
  sex = sample(c("M", "F"), n_samples, replace = TRUE),
  stringsAsFactors = FALSE
)

# Add severity association to some proteins
severity_associated <- c(2, 4, 7, 11, 16, 22, 29)
for (i in severity_associated) {
  severity_effect <- scale(metadata_cont$severity_score) * 15
  count_matrix_cont[i, ] <- count_matrix_cont[i, ] + severity_effect +
    rnorm(n_samples, 0, 10)
}

# Run regression analysis
result_cont <- olink_regression_analysis(
  count_matrix = count_matrix_cont,
  metadata = metadata_cont,
  outcome = "severity_score",
  outcome_type = "continuous",
  covariates = c("age", "sex"),
  sample_id_col = "sample_id",
  p_adjust_method = "BH",
  verbose = TRUE
)

cat("  - Analysis complete\n")

# Create forest plot for continuous outcome
p8 <- plot_forest_results(
  results = result_cont,
  analysis_type = "univariate",
  top_n = 20,
  sort_by = "p_value",
  show_pval = TRUE,
  title = "Beta Coefficients: Severity Score Association (Top 20)"
)

ggsave(
  file.path(output_dir, "forest_plot_continuous.png"),
  plot = p8, width = 12, height = 10, dpi = 300
)
cat("  - Saved: forest_plot_continuous.png\n")

# Create scatter plot for top protein
top_protein <- result_cont$univariate_results$protein[1]
top_protein_data <- data.frame(
  expression = count_matrix_cont[top_protein, ],
  severity = metadata_cont$severity_score
)

p9 <- ggplot(top_protein_data, aes(x = expression, y = severity)) +
  geom_point(alpha = 0.6, size = 3, color = "#3498DB") +
  geom_smooth(method = "lm", se = TRUE, color = "#E74C3C") +
  labs(
    title = paste("Top Association:", top_protein),
    subtitle = sprintf(
      "Beta = %.2f, p = %.2e",
      result_cont$univariate_results$beta[1],
      result_cont$univariate_results$p_value[1]
    ),
    x = "Protein Expression Level",
    y = "Severity Score"
  ) +
  theme_bw()

ggsave(
  file.path(output_dir, "top_protein_association.png"),
  plot = p9, width = 8, height = 6, dpi = 300
)
cat("  - Saved: top_protein_association.png\n")

# Create Manhattan-style plot
manhattan_data <- result_cont$univariate_results
manhattan_data$protein_num <- as.numeric(
  gsub("Protein_", "", manhattan_data$protein)
)

p10 <- ggplot(
  manhattan_data,
  aes(x = protein_num, y = -log10(p_value))
) +
  geom_point(aes(color = p_adjusted < 0.05), size = 3, alpha = 0.7) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             color = "gray50") +
  scale_color_manual(
    values = c("TRUE" = "#E74C3C", "FALSE" = "#95A5A6"),
    labels = c("TRUE" = "FDR < 0.05", "FALSE" = "Not Significant"),
    name = ""
  ) +
  labs(
    title = "Manhattan Plot: Protein-Severity Associations",
    x = "Protein Index",
    y = "-log10(p-value)"
  ) +
  theme_bw() +
  theme(legend.position = "top")

ggsave(
  file.path(output_dir, "manhattan_plot_continuous.png"),
  plot = p10, width = 12, height = 6, dpi = 300
)
cat("  - Saved: manhattan_plot_continuous.png\n\n")

# ==============================================================================
# Generate Summary Statistics Tables
# ==============================================================================
cat("4. Generating summary tables...\n")

# Export top results tables
write.csv(
  head(result_binary$univariate_results, 10),
  file.path(output_dir, "top10_univariate_binary.csv"),
  row.names = FALSE
)

write.csv(
  head(result_binary$multivariate_results, 10),
  file.path(output_dir, "top10_multivariate_binary.csv"),
  row.names = FALSE
)

write.csv(
  head(result_cont$univariate_results, 10),
  file.path(output_dir, "top10_univariate_continuous.csv"),
  row.names = FALSE
)

cat("  - Saved result tables\n")

cat("\n=== All examples generated successfully! ===\n")
cat(paste0("Output directory: ", output_dir, "/\n"))
cat(paste0("Total files created: ", length(list.files(output_dir)), "\n"))
