#!/usr/bin/env Rscript
# Quick demo script for plot_heatmap_with_fc()
# Run this from R console or RStudio

# Load the development version of mestools
devtools::load_all(".")

# ============================================
# SIMPLE DEMO - Basic Usage
# ============================================
cat("\n===========================================\n")
cat("Creating test data...\n")
cat("===========================================\n")

set.seed(123)

# Create gene expression matrix (20 genes x 8 samples)
expression <- matrix(
  rnorm(160, mean = 10, sd = 3),
  nrow = 20,
  ncol = 8
)

# Add differential expression pattern
# Genes 1-10: upregulated in treatment
expression[1:10, 5:8] <- expression[1:10, 5:8] + 4

# Genes 15-20: downregulated in treatment
expression[15:20, 5:8] <- expression[15:20, 5:8] - 3

# Set names
rownames(expression) <- paste0("Gene_", 1:20)
colnames(expression) <- c(paste0("Ctrl_", 1:4), paste0("Treat_", 1:4))

# Create log2 fold change
log2fc <- c(
  rnorm(10, mean = 2.5, sd = 0.4),   # Upregulated
  rnorm(4, mean = 0.1, sd = 0.2),    # No change
  rnorm(6, mean = -2, sd = 0.3)      # Downregulated
)

# Sample groups
groups <- c(rep("Control", 4), rep("Treatment", 4))

cat("\nData summary:\n")
cat("  - Genes: ", nrow(expression), "\n")
cat("  - Samples: ", ncol(expression), "\n")
cat("  - Groups: ", unique(groups), "\n")
cat("  - log2FC range: ", round(range(log2fc), 2), "\n\n")

cat("===========================================\n")
cat("Generating plot...\n")
cat("===========================================\n\n")

# Create the plot
result <- plot_heatmap_with_fc(
  expression_matrix = expression,
  log2fc = log2fc,
  sample_groups = groups,
  fc_threshold = 1.5,
  title = "Example: Differential Gene Expression",
  scale = "row"
)

cat("\n===========================================\n")
cat("SUCCESS! Plot generated.\n")
cat("===========================================\n\n")

cat("What you should see:\n")
cat("  LEFT: Heatmap showing expression patterns\n")
cat("        - Blue = low expression\n")
cat("        - Red = high expression\n")
cat("        - Control vs Treatment groups annotated\n\n")
cat("  RIGHT: Log2FC bars\n")
cat("        - Red bars (right) = upregulated genes\n")
cat("        - Blue bars (left) = downregulated genes\n")
cat("        - Dashed lines at +/- 1.5 threshold\n\n")

cat("Returned data:\n")
cat("  - result$heatmap: pheatmap object\n")
cat("  - result$fc_plot: ggplot2 object\n")
cat("  - result$fc_data: fold change data frame\n")
cat("  - result$row_order: gene clustering order\n\n")

# Show the first few rows of FC data
cat("First 5 genes in plot order:\n")
print(head(result$fc_data, 5))

cat("\n\nTo save this plot:\n")
cat('  pdf("my_plot.pdf", width = 10, height = 8)\n')
cat('  plot_heatmap_with_fc(...)\n')
cat('  dev.off()\n\n')
