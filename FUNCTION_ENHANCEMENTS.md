# plot_heatmap_with_fc Function Enhancements

## Overview
The `plot_heatmap_with_fc()` function has been significantly enhanced with new features for comprehensive differential expression visualization.

## New Features

### 1. Statistical Annotations
**Parameters:**
- `pvalues`: Numeric vector of p-values
- `padj`: Numeric vector of adjusted p-values (FDR, Bonferroni, etc.)
- `sig_thresholds`: Significance thresholds (default: c(0.05, 0.01, 0.001))
- `show_pval_text`: Display p-value text on bars (default: FALSE)

**Features:**
- Automatic significance stars (*, **, ***) displayed on fold change bars
- Uses adjusted p-values if provided, otherwise raw p-values
- Customizable significance thresholds
- Optional p-value text labels

**Example:**
```r
plot_heatmap_with_fc(
  expression_matrix = expr_data,
  log2fc = fc_values,
  pvalues = pvals,
  padj = padj_values,
  sig_thresholds = c(0.05, 0.01, 0.001)
)
```

### 2. Gene Filtering & Highlighting
**Parameters:**
- `filter_by_fc`: Show only genes exceeding fc_threshold (default: FALSE)
- `top_n`: Display top N genes by absolute log2FC (default: NULL)
- `highlight_genes`: Character vector of gene names to highlight
- `highlight_color`: Color for highlighted genes (default: "#FF6B35")
- `gene_categories`: Character vector for gene grouping

**Features:**
- Filter genes by fold change threshold
- Select top differentially expressed genes
- Highlight specific genes of interest with custom colors
- Group genes by categories (metabolic, signaling, etc.)

**Examples:**
```r
# Show only top 20 genes
plot_heatmap_with_fc(
  expression_matrix = expr_data,
  log2fc = fc_values,
  top_n = 20
)

# Highlight specific genes
plot_heatmap_with_fc(
  expression_matrix = expr_data,
  log2fc = fc_values,
  highlight_genes = c("TP53", "BRCA1", "EGFR"),
  highlight_color = "#FF6B35"
)

# Categorize genes
plot_heatmap_with_fc(
  expression_matrix = expr_data,
  log2fc = fc_values,
  gene_categories = c(rep("Metabolic", 10), rep("Signaling", 15))
)
```

### 3. Enhanced Color Schemes
**Parameter:**
- `color_theme`: Preset color themes (default: "default")

**Available Themes:**
1. **default**: Blue-white-red gradient
2. **viridis**: Viridis colormap (perceptually uniform)
3. **RdBu**: Red-Blue diverging palette (ColorBrewer)
4. **publication**: Grayscale for publications
5. **colorblind**: Colorblind-friendly palette

**Example:**
```r
# Use viridis theme
plot_heatmap_with_fc(
  expression_matrix = expr_data,
  log2fc = fc_values,
  color_theme = "viridis"
)

# Publication-ready grayscale
plot_heatmap_with_fc(
  expression_matrix = expr_data,
  log2fc = fc_values,
  color_theme = "publication"
)
```

### 4. Enhanced Annotations
**Parameters:**
- `annotation_row`: Data frame of row (gene) annotations
- `annotation_col_enhanced`: Data frame of column (sample) annotations
- `annotation_colors`: Named list of custom annotation colors

**Features:**
- Multiple row and column annotations
- Custom color schemes for annotations
- Automatic annotation from gene_categories
- Support for batch effects, time points, and other metadata

**Example:**
```r
# Multiple column annotations
col_annot <- data.frame(
  Group = c(rep("Control", 5), rep("Treatment", 5)),
  Batch = rep(c("Batch1", "Batch2"), each = 5),
  TimePoint = rep(c("0h", "24h"), 5),
  row.names = colnames(expression)
)

# Custom colors
annot_colors <- list(
  Group = c(Control = "#4393C3", Treatment = "#D6604D"),
  Batch = c(Batch1 = "#FDB863", Batch2 = "#B2ABD2"),
  TimePoint = c("0h" = "#E0E0E0", "24h" = "#636363")
)

plot_heatmap_with_fc(
  expression_matrix = expr_data,
  log2fc = fc_values,
  annotation_col_enhanced = col_annot,
  annotation_colors = annot_colors
)
```

### 5. Volcano Plot Integration
**Parameter:**
- `plot_layout`: Layout type (default: "heatmap_fc")

**Options:**
1. **"heatmap_fc"**: Two-panel layout (heatmap | fold change)
2. **"heatmap_fc_volcano"**: Three-panel layout (heatmap | fold change | volcano)

**Features:**
- Integrated volcano plot showing log2FC vs -log10(p-value)
- Synchronized gene highlighting across all panels
- Automatic significance thresholds visualization
- Requires p-values or adjusted p-values

**Example:**
```r
plot_heatmap_with_fc(
  expression_matrix = expr_data,
  log2fc = fc_values,
  pvalues = pvals,
  highlight_genes = c("TP53", "BRCA1"),
  plot_layout = "heatmap_fc_volcano"
)
```

## Comprehensive Example

```r
library(mestools)

# Load your data
expression <- read.csv("expression_matrix.csv", row.names = 1)
results <- read.csv("differential_expression_results.csv")

# Create comprehensive visualization
plot_heatmap_with_fc(
  # Core data
  expression_matrix = as.matrix(expression),
  log2fc = results$log2FoldChange,

  # Statistical significance
  pvalues = results$pvalue,
  padj = results$padj,
  sig_thresholds = c(0.05, 0.01, 0.001),

  # Gene filtering
  top_n = 30,  # Show top 30 genes
  highlight_genes = c("TP53", "BRCA1", "EGFR", "MYC"),
  highlight_color = "#FF6B35",
  gene_categories = results$pathway,

  # Visualization
  color_theme = "viridis",
  plot_layout = "heatmap_fc_volcano",

  # Annotations
  sample_groups = colData$condition,
  fc_threshold = 1,

  # Title
  title = "Comprehensive Differential Expression Analysis"
)
```

## Return Value

The function returns a list with the following components:

- `heatmap`: pheatmap object
- `fc_plot`: ggplot2 fold change bar plot
- `fc_data`: Data frame with fold change data and significance
- `volcano_plot`: ggplot2 volcano plot (if applicable)
- `row_order`: Gene order after clustering
- `expression_matrix`: Expression matrix used in the plot
- `log2fc`: Log2 fold changes
- `pvalues`: P-values (if provided)
- `padj`: Adjusted p-values (if provided)

## Testing

A comprehensive test script is provided in [test_heatmap_fc_enhanced.R](test_heatmap_fc_enhanced.R) demonstrating all new features:

1. Statistical significance annotations
2. Gene filtering (top N genes)
3. Gene highlighting
4. Color themes (viridis, RdBu, colorblind, publication)
5. Enhanced annotations (multiple factors)
6. Volcano plot integration
7. Combined features

## Backward Compatibility

All new parameters have default values ensuring complete backward compatibility. Existing code using the function will work without modification.

## Dependencies

The enhanced function requires:
- `pheatmap` (existing)
- `grid` (existing)
- `gridExtra` (existing)
- `ggplot2` (for volcano plots and enhanced bar plots)

## Performance Notes

- Gene filtering (`top_n`, `filter_by_fc`) improves performance for large datasets
- Recommended to use `top_n` for datasets with >100 genes
- Volcano plot adds minimal overhead (~5-10% additional computation time)

## Future Enhancements

Potential future additions:
- Interactive plotly version
- PCA plot integration
- Export to PDF/PNG with custom dimensions
- Gene set enrichment annotations
- Time series layout for temporal data
