# plot_heatmap_with_fc() Quick Reference

## Basic Usage
```r
plot_heatmap_with_fc(
  expression_matrix = my_expr_matrix,
  log2fc = my_log2fc_values
)
```

## Feature Comparison

| Feature | Old Version | New Version |
|---------|-------------|-------------|
| Basic heatmap + FC | ✓ | ✓ |
| Statistical significance | ✗ | ✓ (*, **, ***) |
| Gene filtering | ✗ | ✓ (top_n, filter_by_fc) |
| Gene highlighting | ✗ | ✓ |
| Color themes | 1 (default) | 5 (default, viridis, RdBu, publication, colorblind) |
| Volcano plot | ✗ | ✓ |
| Multi-panel layouts | 1 (heatmap+FC) | 2 (2-panel, 3-panel) |
| Row annotations | Basic | Enhanced (multiple factors) |
| Column annotations | Basic (sample groups) | Enhanced (multiple factors) |

## Quick Examples

### Add Significance Stars
```r
plot_heatmap_with_fc(expr, fc, pvalues = pvals, padj = padj)
```

### Show Top 20 Genes
```r
plot_heatmap_with_fc(expr, fc, top_n = 20)
```

### Highlight Specific Genes
```r
plot_heatmap_with_fc(expr, fc, highlight_genes = c("TP53", "BRCA1"))
```

### Use Different Color Theme
```r
plot_heatmap_with_fc(expr, fc, color_theme = "viridis")
```

### Add Volcano Plot
```r
plot_heatmap_with_fc(expr, fc, pvalues = pvals, plot_layout = "heatmap_fc_volcano")
```

### Everything Combined
```r
plot_heatmap_with_fc(
  expression_matrix = expr,
  log2fc = fc,
  pvalues = pvals,
  padj = padj,
  top_n = 30,
  highlight_genes = c("TP53", "BRCA1", "EGFR"),
  color_theme = "viridis",
  plot_layout = "heatmap_fc_volcano",
  title = "My Analysis"
)
```

## Parameter Quick Reference

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `expression_matrix` | matrix | - | Gene expression data (required) |
| `log2fc` | numeric | - | Log2 fold changes (required) |
| `pvalues` | numeric | NULL | P-values for significance |
| `padj` | numeric | NULL | Adjusted p-values |
| `top_n` | integer | NULL | Show top N genes |
| `filter_by_fc` | logical | FALSE | Filter by FC threshold |
| `highlight_genes` | character | NULL | Genes to highlight |
| `color_theme` | character | "default" | Color theme |
| `plot_layout` | character | "heatmap_fc" | Layout type |
| `sample_groups` | character | NULL | Sample grouping |
| `fc_threshold` | numeric | 1 | FC significance threshold |
| `title` | character | NULL | Plot title |

## Color Themes

```r
# Default: Blue-white-red
color_theme = "default"

# Viridis: Perceptually uniform
color_theme = "viridis"

# RdBu: Red-Blue diverging
color_theme = "RdBu"

# Publication: Grayscale
color_theme = "publication"

# Colorblind-friendly
color_theme = "colorblind"
```

## Common Use Cases

### Publication Figure
```r
plot_heatmap_with_fc(
  expr, fc,
  pvalues = pvals,
  top_n = 25,
  color_theme = "publication",
  title = "Differentially Expressed Genes"
)
```

### Exploratory Analysis
```r
plot_heatmap_with_fc(
  expr, fc,
  pvalues = pvals,
  plot_layout = "heatmap_fc_volcano",
  color_theme = "viridis"
)
```

### Focus on Key Genes
```r
my_genes <- c("TP53", "BRCA1", "EGFR", "MYC", "KRAS")
plot_heatmap_with_fc(
  expr, fc,
  highlight_genes = my_genes,
  highlight_color = "#FF6B35"
)
```

### Large Dataset
```r
plot_heatmap_with_fc(
  expr, fc,
  top_n = 50,  # Reduce to top 50 for clarity
  show_rownames = TRUE
)
```

## Tips

1. **For large datasets**: Use `top_n` to focus on most significant genes
2. **For presentations**: Use `color_theme = "colorblind"` for accessibility
3. **For publications**: Use `color_theme = "publication"` for grayscale
4. **For exploratory analysis**: Use `plot_layout = "heatmap_fc_volcano"`
5. **Always provide p-values**: Get significance stars automatically

## Troubleshooting

**Issue**: Too many genes, plot is crowded
```r
# Solution: Filter to top genes
plot_heatmap_with_fc(expr, fc, top_n = 30)
```

**Issue**: Can't see gene names
```r
# Solution: Force display or reduce gene count
plot_heatmap_with_fc(expr, fc, show_rownames = TRUE, top_n = 40)
```

**Issue**: Volcano plot not showing
```r
# Solution: Make sure you provide p-values
plot_heatmap_with_fc(expr, fc, pvalues = pvals, plot_layout = "heatmap_fc_volcano")
```

**Issue**: Want custom colors
```r
# Solution: Use heatmap_colors parameter
my_colors <- colorRampPalette(c("blue", "yellow", "red"))(100)
plot_heatmap_with_fc(expr, fc, heatmap_colors = my_colors)
```
