# Combined Heatmap and Log2 Fold Change Plot

## Overview

The `plot_heatmap_with_fc()` function creates a publication-ready visualization that combines:
- **Left panel**: Heatmap of gene expression data (with hierarchical clustering)
- **Right panel**: Bar plot of log2 fold change values

This is particularly useful for visualizing differential gene expression analysis results.

## Installation

The function requires the following packages:
```r
install.packages(c("pheatmap", "ggplot2", "grid", "gridExtra"))
```

## Basic Usage

```r
library(mestools)

# Create sample data
expression <- matrix(rnorm(200), nrow = 20, ncol = 10)
rownames(expression) <- paste0("Gene", 1:20)
colnames(expression) <- paste0("Sample", 1:10)

log2fc <- rnorm(20, mean = 0, sd = 2)
groups <- rep(c("Control", "Treatment"), each = 5)

# Create the plot
plot_heatmap_with_fc(
  expression_matrix = expression,
  log2fc = log2fc,
  sample_groups = groups,
  fc_threshold = 1.5,
  title = "Differential Gene Expression"
)
```

## Key Parameters

### Required Parameters

- `expression_matrix`: Numeric matrix with genes as rows and samples as columns
- `log2fc`: Numeric vector of log2 fold change values (must match row count)

### Optional Parameters

- `sample_groups`: Character vector for sample group annotations (e.g., "Control", "Treatment")
- `cluster_rows`: Cluster genes by expression pattern (default: TRUE)
- `cluster_cols`: Cluster samples by expression pattern (default: TRUE)
- `scale`: Scale data by "row", "column", or "none" (default: "row")
- `show_rownames`: Display gene names (default: TRUE, auto-hidden if >50 genes)
- `show_colnames`: Display sample names (default: TRUE)
- `fc_threshold`: Threshold for highlighting significant fold changes (default: 1)
- `heatmap_colors`: Color palette for heatmap (default: blue-white-red)
- `fc_colors`: Colors for negative and positive fold changes (default: blue and red)
- `title`: Main plot title
- `gene_labels`: Custom gene labels (uses rownames if NULL)
- `width_ratio`: Width ratio of heatmap to FC plot (default: c(4, 1))

## Advanced Examples

### Example 1: Custom Colors and Thresholds

```r
plot_heatmap_with_fc(
  expression_matrix = expression,
  log2fc = log2fc,
  sample_groups = groups,
  fc_threshold = 2,  # More stringent threshold
  heatmap_colors = colorRampPalette(c("navy", "white", "firebrick"))(100),
  fc_colors = c("#4575B4", "#D73027"),  # Custom colors
  title = "High-Confidence Differential Expression",
  width_ratio = c(3, 1)  # Narrower FC panel
)
```

### Example 2: No Clustering, Show All Gene Names

```r
plot_heatmap_with_fc(
  expression_matrix = expression,
  log2fc = log2fc,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  scale = "none",  # No scaling
  title = "Raw Expression Values"
)
```

### Example 3: Real Gene Names

```r
# With real gene symbols
gene_names <- c("BRCA1", "TP53", "EGFR", "MYC", "KRAS", ...)

plot_heatmap_with_fc(
  expression_matrix = expression,
  log2fc = log2fc,
  sample_groups = groups,
  gene_labels = gene_names,
  show_rownames = TRUE,
  fc_threshold = 1,
  title = "Cancer Gene Expression Panel"
)
```

### Example 4: Wide Format for Presentation

```r
# Adjust width ratio for different emphasis
plot_heatmap_with_fc(
  expression_matrix = expression,
  log2fc = log2fc,
  sample_groups = groups,
  width_ratio = c(5, 1),  # Wider heatmap
  title = "Comprehensive Expression Analysis"
)
```

## Understanding the Output

The function returns an invisible list containing:

```r
result <- plot_heatmap_with_fc(...)

# Access components:
result$heatmap          # pheatmap object with clustering info
result$fc_plot          # ggplot2 object of fold change plot
result$fc_data          # Data frame used for FC plot
result$row_order        # Order of genes after clustering
result$expression_matrix # Original expression matrix
result$log2fc           # Original log2FC values
```

## Tips and Best Practices

1. **Scaling**: Use `scale = "row"` to normalize each gene across samples (recommended for most cases)

2. **Gene names**: For >50 genes, row names are automatically hidden. Override with `show_rownames = TRUE`

3. **Clustering**: The log2FC plot automatically matches the clustered order from the heatmap

4. **Color schemes**:
   - Blue-white-red is standard for expression (down-neutral-up)
   - Customize with `colorRampPalette()` for specific needs

5. **Threshold lines**: Dashed lines in the FC plot indicate the fold change threshold

6. **Sample groups**: Provide `sample_groups` to add color-coded column annotations

## Integration with Differential Expression Tools

### With DESeq2

```r
library(DESeq2)

# After running DESeq2 analysis
res <- results(dds)
res_sig <- subset(res, padj < 0.05)

# Get normalized counts for significant genes
counts_norm <- counts(dds, normalized = TRUE)
counts_sig <- counts_norm[rownames(res_sig), ]

# Create the plot
plot_heatmap_with_fc(
  expression_matrix = counts_sig,
  log2fc = res_sig$log2FoldChange,
  sample_groups = dds$condition,
  title = "DESeq2 Results: Significant Genes"
)
```

### With edgeR

```r
library(edgeR)

# After running edgeR analysis
top_genes <- topTags(et, n = 50)$table

# Get normalized expression
cpm_values <- cpm(dge, log = TRUE)
cpm_sig <- cpm_values[rownames(top_genes), ]

# Create the plot
plot_heatmap_with_fc(
  expression_matrix = cpm_sig,
  log2fc = top_genes$logFC,
  sample_groups = groups,
  title = "edgeR Results: Top 50 Genes"
)
```

## Saving Plots

```r
# Save to PDF
pdf("combined_heatmap_fc.pdf", width = 10, height = 8)
plot_heatmap_with_fc(
  expression_matrix = expression,
  log2fc = log2fc,
  sample_groups = groups
)
dev.off()

# Save to PNG
png("combined_heatmap_fc.png", width = 1200, height = 800, res = 120)
plot_heatmap_with_fc(
  expression_matrix = expression,
  log2fc = log2fc,
  sample_groups = groups
)
dev.off()
```

## Troubleshooting

**Issue**: Gene names overlap
- **Solution**: Set `show_rownames = FALSE` or reduce font size with pheatmap parameters

**Issue**: Colors don't show pattern
- **Solution**: Try different scaling: `scale = "row"` usually works best

**Issue**: FC plot doesn't align with heatmap
- **Solution**: This shouldn't happen as clustering is automatic, but check that log2fc length matches matrix rows

**Issue**: Plot is too wide/narrow
- **Solution**: Adjust `width_ratio` parameter, e.g., `c(3, 1)` or `c(5, 1)`
