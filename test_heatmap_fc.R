# Test script for plot_heatmap_with_fc function
# This demonstrates how to use the new combined heatmap and log2FC plot

# Install required packages if needed
required_pkgs <- c("pheatmap", "ggplot2", "grid", "gridExtra")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# Load the mestools package
library(mestools)

# ========================================
# Example 1: Basic usage with simulated data
# ========================================
set.seed(123)

# Create simulated gene expression data
# 30 genes x 10 samples
n_genes <- 30
n_samples <- 10

expression <- matrix(
  rnorm(n_genes * n_samples, mean = 10, sd = 3),
  nrow = n_genes,
  ncol = n_samples
)

# Add some differential expression pattern
# First 15 genes are upregulated in treatment group (samples 6-10)
expression[1:15, 6:10] <- expression[1:15, 6:10] + rnorm(15 * 5, mean = 4, sd = 1)

# Last 10 genes are downregulated in treatment group
expression[21:30, 6:10] <- expression[21:30, 6:10] - rnorm(10 * 5, mean = 3, sd = 1)

# Set row and column names
rownames(expression) <- paste0("Gene", 1:n_genes)
colnames(expression) <- paste0("Sample", 1:n_samples)

# Create log2 fold change values
# Simulating comparison of Treatment vs Control
log2fc <- c(
  rnorm(15, mean = 2, sd = 0.5),    # Upregulated genes
  rnorm(5, mean = 0.2, sd = 0.3),   # No change
  rnorm(10, mean = -1.8, sd = 0.4)  # Downregulated genes
)

# Define sample groups
sample_groups <- c(rep("Control", 5), rep("Treatment", 5))

# Create the combined plot
cat("\n========================================\n")
cat("Example 1: Basic heatmap with log2FC\n")
cat("========================================\n")

result <- plot_heatmap_with_fc(
  expression_matrix = expression,
  log2fc = log2fc,
  sample_groups = sample_groups,
  fc_threshold = 1,
  title = "Differential Gene Expression Analysis",
  scale = "row"
)

# ========================================
# Example 2: Larger dataset without clustering
# ========================================
cat("\n========================================\n")
cat("Example 2: No clustering, different colors\n")
cat("========================================\n")

set.seed(456)
n_genes2 <- 50
expression2 <- matrix(
  rnorm(n_genes2 * 8, mean = 8, sd = 2),
  nrow = n_genes2,
  ncol = 8
)
rownames(expression2) <- paste0("Gene", 1:n_genes2)
colnames(expression2) <- paste0("S", 1:8)

log2fc2 <- rnorm(n_genes2, mean = 0, sd = 1.5)
groups2 <- rep(c("WT", "Mutant"), each = 4)

result2 <- plot_heatmap_with_fc(
  expression_matrix = expression2,
  log2fc = log2fc2,
  sample_groups = groups2,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  fc_threshold = 1.5,
  heatmap_colors = colorRampPalette(c("blue", "white", "red"))(100),
  fc_colors = c("#4575B4", "#D73027"),
  title = "Gene Expression - No Clustering",
  width_ratio = c(3, 1)
)

# ========================================
# Example 3: With custom gene labels
# ========================================
cat("\n========================================\n")
cat("Example 3: Custom gene labels\n")
cat("========================================\n")

set.seed(789)
n_genes3 <- 20
expression3 <- matrix(
  rnorm(n_genes3 * 6, mean = 12, sd = 2.5),
  nrow = n_genes3,
  ncol = 6
)

# Create meaningful gene names
gene_names <- c(
  "BRCA1", "TP53", "EGFR", "MYC", "KRAS",
  "PTEN", "PIK3CA", "AKT1", "BRAF", "ERBB2",
  "VEGFA", "CDK4", "CDKN2A", "RB1", "NRAS",
  "FGFR1", "MET", "ALK", "RET", "NOTCH1"
)

colnames(expression3) <- paste0("Patient", 1:6)
groups3 <- rep(c("Normal", "Tumor"), each = 3)

log2fc3 <- rnorm(n_genes3, mean = 0.5, sd = 1.2)

result3 <- plot_heatmap_with_fc(
  expression_matrix = expression3,
  log2fc = log2fc3,
  sample_groups = groups3,
  gene_labels = gene_names,
  show_rownames = TRUE,
  fc_threshold = 0.8,
  title = "Cancer Gene Expression Panel"
)

cat("\n========================================\n")
cat("Testing completed successfully!\n")
cat("========================================\n")

# Show what's returned by the function
cat("\nReturned object structure:\n")
str(result, max.level = 1)
