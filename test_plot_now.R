# Test script for plot_heatmap_with_fc()
# Run this in RStudio or R console

# Set working directory
setwd("c:/Users/tranh/OneDrive/mestools")

# Load the package in development mode
devtools::load_all(".")

# Install required packages if needed
required_pkgs <- c("pheatmap", "ggplot2", "grid", "gridExtra")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# ============================================
# Create sample data
# ============================================
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

# Set gene names
rownames(expression) <- c(
  "BRCA1", "TP53", "EGFR", "MYC", "KRAS",
  "PTEN", "PIK3CA", "AKT1", "BRAF", "ERBB2",
  "Gene11", "Gene12", "Gene13", "Gene14",
  "VEGFA", "CDK4", "CDKN2A", "RB1", "NRAS", "FGFR1"
)

colnames(expression) <- c(
  paste0("Ctrl_", 1:4),
  paste0("Treat_", 1:4)
)

# Create log2 fold change
log2fc <- c(
  rnorm(10, mean = 2.5, sd = 0.4),   # Upregulated
  rnorm(4, mean = 0.1, sd = 0.2),    # No change
  rnorm(6, mean = -2, sd = 0.3)      # Downregulated
)

# Sample groups
groups <- c(rep("Control", 4), rep("Treatment", 4))

# ============================================
# Generate the plot
# ============================================
cat("\n===========================================\n")
cat("Generating combined heatmap + log2FC plot\n")
cat("===========================================\n\n")

result <- plot_heatmap_with_fc(
  expression_matrix = expression,
  log2fc = log2fc,
  sample_groups = groups,
  fc_threshold = 1.5,
  title = "Differential Gene Expression",
  scale = "row",
  show_rownames = TRUE
)

cat("\n===========================================\n")
cat("Plot generated successfully!\n")
cat("===========================================\n\n")

cat("What you should see:\n")
cat("  LEFT PANEL: Heatmap with gene names\n")
cat("    - Blue = low expression\n")
cat("    - Red = high expression\n")
cat("    - Genes clustered by similarity\n\n")
cat("  RIGHT PANEL: Log2FC bars with gene names\n")
cat("    - Red bars = upregulated genes\n")
cat("    - Blue bars = downregulated genes\n")
cat("    - Gene names match left panel\n\n")

# Show summary of data
cat("\nData Summary:\n")
cat("  Genes:", nrow(expression), "\n")
cat("  Samples:", ncol(expression), "\n")
cat("  Groups:", unique(groups), "\n")
cat("  log2FC range:", round(range(log2fc), 2), "\n")

# Show first few genes in plot order
cat("\nFirst 5 genes (in clustered order):\n")
print(head(result$fc_data, 5))

cat("\n\nTo save this plot:\n")
cat("  pdf('my_heatmap_fc.pdf', width = 10, height = 8)\n")
cat("  plot_heatmap_with_fc(...)\n")
cat("  dev.off()\n\n")
