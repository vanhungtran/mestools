# Enhanced Test Script for plot_heatmap_with_fc function
# Demonstrates new features: statistical annotations, gene filtering,
# volcano plots, color themes, and enhanced annotations

# Load the mestools package
devtools::load_all()

# Check if required packages are available
required_pkgs <- c("pheatmap", "ggplot2", "grid", "gridExtra")
missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]

if (length(missing_pkgs) > 0) {
  cat("Installing required packages:", paste(missing_pkgs, collapse = ", "), "\n")
  install.packages(missing_pkgs)
}

# ========================================
# Example 1: Statistical Annotations
# ========================================
cat("\n========================================\n")
cat("Example 1: Statistical Significance\n")
cat("========================================\n")

set.seed(123)
n_genes <- 30
n_samples <- 10

# Create expression data
expression <- matrix(
  rnorm(n_genes * n_samples, mean = 10, sd = 3),
  nrow = n_genes,
  ncol = n_samples
)
rownames(expression) <- paste0("Gene", 1:n_genes)
colnames(expression) <- paste0("Sample", 1:n_samples)

# Add differential expression
expression[1:15, 6:10] <- expression[1:15, 6:10] + rnorm(15 * 5, mean = 4, sd = 1)
expression[21:30, 6:10] <- expression[21:30, 6:10] - rnorm(10 * 5, mean = 3, sd = 1)

# Create log2FC and p-values
log2fc <- c(
  rnorm(15, mean = 2, sd = 0.5),
  rnorm(5, mean = 0.2, sd = 0.3),
  rnorm(10, mean = -1.8, sd = 0.4)
)

# Simulate p-values (more significant for higher |log2FC|)
pvalues <- exp(-abs(log2fc) * 2) * runif(n_genes, 0.001, 0.1)
padj <- p.adjust(pvalues, method = "fdr")

sample_groups <- c(rep("Control", 5), rep("Treatment", 5))

# Plot with significance stars
result1 <- plot_heatmap_with_fc(
  expression_matrix = expression,
  log2fc = log2fc,
  pvalues = pvalues,
  padj = padj,
  sample_groups = sample_groups,
  fc_threshold = 1,
  title = "Expression with Significance Annotations",
  scale = "row"
)

cat("  - Significance stars displayed on FC bars\n")
cat("  - *** p < 0.001, ** p < 0.01, * p < 0.05\n")

# ========================================
# Example 2: Gene Filtering (Top N)
# ========================================
cat("\n========================================\n")
cat("Example 2: Top 15 Genes by |log2FC|\n")
cat("========================================\n")

result2 <- plot_heatmap_with_fc(
  expression_matrix = expression,
  log2fc = log2fc,
  pvalues = pvalues,
  padj = padj,
  sample_groups = sample_groups,
  top_n = 15,
  fc_threshold = 1,
  title = "Top 15 Differentially Expressed Genes",
  scale = "row"
)

# ========================================
# Example 3: Gene Highlighting
# ========================================
cat("\n========================================\n")
cat("Example 3: Highlighted Genes of Interest\n")
cat("========================================\n")

genes_of_interest <- c("Gene5", "Gene8", "Gene15", "Gene25")

result3 <- plot_heatmap_with_fc(
  expression_matrix = expression,
  log2fc = log2fc,
  pvalues = pvalues,
  sample_groups = sample_groups,
  highlight_genes = genes_of_interest,
  highlight_color = "#FF6B35",
  fc_threshold = 1,
  title = "Highlighted Genes of Interest",
  scale = "row"
)

cat("  - Genes", paste(genes_of_interest, collapse = ", "), "are highlighted\n")

# ========================================
# Example 4: Color Themes
# ========================================
cat("\n========================================\n")
cat("Example 4: Different Color Themes\n")
cat("========================================\n")

# Viridis theme
cat("  4a. Viridis theme...\n")
result4a <- plot_heatmap_with_fc(
  expression_matrix = expression,
  log2fc = log2fc,
  sample_groups = sample_groups,
  color_theme = "viridis",
  title = "Viridis Color Theme",
  scale = "row"
)

# RdBu theme
cat("  4b. RdBu theme...\n")
result4b <- plot_heatmap_with_fc(
  expression_matrix = expression,
  log2fc = log2fc,
  sample_groups = sample_groups,
  color_theme = "RdBu",
  title = "RdBu Color Theme",
  scale = "row"
)

# Colorblind-friendly theme
cat("  4c. Colorblind-friendly theme...\n")
result4c <- plot_heatmap_with_fc(
  expression_matrix = expression,
  log2fc = log2fc,
  sample_groups = sample_groups,
  color_theme = "colorblind",
  title = "Colorblind-Friendly Theme",
  scale = "row"
)

# Publication theme (grayscale)
cat("  4d. Publication theme (grayscale)...\n")
result4d <- plot_heatmap_with_fc(
  expression_matrix = expression,
  log2fc = log2fc,
  sample_groups = sample_groups,
  color_theme = "publication",
  title = "Publication Theme (Grayscale)",
  scale = "row"
)

# ========================================
# Example 5: Enhanced Annotations
# ========================================
cat("\n========================================\n")
cat("Example 5: Enhanced Row and Column Annotations\n")
cat("========================================\n")

# Create gene categories
gene_categories <- c(
  rep("Metabolic", 10),
  rep("Signaling", 10),
  rep("Structural", 10)
)

# Enhanced column annotations with multiple factors
annotation_col_enhanced <- data.frame(
  Group = sample_groups,
  Batch = rep(c("Batch1", "Batch2"), each = 5),
  TimePoint = rep(c("0h", "24h"), 5),
  row.names = colnames(expression)
)

# Custom annotation colors
annotation_colors <- list(
  Group = c(Control = "#4393C3", Treatment = "#D6604D"),
  Batch = c(Batch1 = "#FDB863", Batch2 = "#B2ABD2"),
  TimePoint = c("0h" = "#E0E0E0", "24h" = "#636363"),
  Category = c(Metabolic = "#66C2A5", Signaling = "#FC8D62", Structural = "#8DA0CB")
)

result5 <- plot_heatmap_with_fc(
  expression_matrix = expression,
  log2fc = log2fc,
  pvalues = pvalues,
  annotation_col_enhanced = annotation_col_enhanced,
  gene_categories = gene_categories,
  annotation_colors = annotation_colors,
  fc_threshold = 1,
  title = "Enhanced Annotations",
  scale = "row"
)

cat("  - Column annotations: Group, Batch, TimePoint\n")
cat("  - Row annotations: Gene categories\n")

# ========================================
# Example 6: Volcano Plot Layout
# ========================================
cat("\n========================================\n")
cat("Example 6: Three-Panel Layout with Volcano Plot\n")
cat("========================================\n")

result6 <- plot_heatmap_with_fc(
  expression_matrix = expression,
  log2fc = log2fc,
  pvalues = pvalues,
  padj = padj,
  sample_groups = sample_groups,
  highlight_genes = genes_of_interest,
  plot_layout = "heatmap_fc_volcano",
  fc_threshold = 1,
  title = "Heatmap + FC + Volcano",
  scale = "row"
)

cat("  - Three panel layout: Heatmap | Fold Change | Volcano\n")
cat("  - Highlighted genes shown in volcano plot\n")

# ========================================
# Example 7: Combined Features
# ========================================
cat("\n========================================\n")
cat("Example 7: All Features Combined\n")
cat("========================================\n")

result7 <- plot_heatmap_with_fc(
  expression_matrix = expression,
  log2fc = log2fc,
  pvalues = pvalues,
  padj = padj,
  sample_groups = sample_groups,
  highlight_genes = genes_of_interest,
  highlight_color = "#FF6B35",
  gene_categories = gene_categories,
  annotation_colors = annotation_colors,
  color_theme = "RdBu",
  plot_layout = "heatmap_fc_volcano",
  top_n = 20,
  fc_threshold = 1,
  sig_thresholds = c(0.05, 0.01, 0.001),
  title = "Comprehensive Analysis (All Features)",
  scale = "row"
)

cat("  - Top 20 genes by |log2FC|\n")
cat("  - Significance stars on FC plot\n")
cat("  - Gene highlighting and categories\n")
cat("  - Three-panel layout with volcano plot\n")
cat("  - RdBu color theme\n")

# ========================================
# Summary
# ========================================
cat("\n========================================\n")
cat("Testing completed successfully!\n")
cat("========================================\n")
cat("\nNew features demonstrated:\n")
cat("1. Statistical significance annotations (*, **, ***)\n")
cat("2. Gene filtering (top_n, filter_by_fc)\n")
cat("3. Gene highlighting with custom colors\n")
cat("4. Multiple color themes (viridis, RdBu, colorblind, publication)\n")
cat("5. Enhanced row/column annotations\n")
cat("6. Volcano plot integration\n")
cat("7. Three-panel layout options\n")
cat("\nReturned object structure:\n")
str(result7, max.level = 1)
