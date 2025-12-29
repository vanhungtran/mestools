# plot_heatmap_with_fc() - Quick Reference

## Minimal Example
```r
plot_heatmap_with_fc(
  expression_matrix = my_expression,  # Matrix: genes × samples
  log2fc = my_log2fc                  # Vector: same length as rows
)
```

## Common Patterns

### Pattern 1: Basic DE Visualization
```r
plot_heatmap_with_fc(
  expression_matrix = expr_data,
  log2fc = fc_values,
  sample_groups = c(rep("Control", 5), rep("Treatment", 5)),
  title = "DEG Analysis"
)
```

### Pattern 2: High Stringency
```r
plot_heatmap_with_fc(
  expression_matrix = expr_data,
  log2fc = fc_values,
  fc_threshold = 2,           # |log2FC| > 2
  cluster_rows = TRUE,
  scale = "row"
)
```

### Pattern 3: No Clustering
```r
plot_heatmap_with_fc(
  expression_matrix = expr_data,
  log2fc = fc_values,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  scale = "none"
)
```

### Pattern 4: Custom Colors
```r
plot_heatmap_with_fc(
  expression_matrix = expr_data,
  log2fc = fc_values,
  heatmap_colors = colorRampPalette(c("green", "black", "red"))(100),
  fc_colors = c("darkblue", "darkred")
)
```

### Pattern 5: Wide FC Panel
```r
plot_heatmap_with_fc(
  expression_matrix = expr_data,
  log2fc = fc_values,
  width_ratio = c(3, 1.5)    # Wider FC panel
)
```

## Parameter Quick Lookup

| Parameter | Type | Default | Purpose |
|-----------|------|---------|---------|
| `expression_matrix` | matrix | **required** | Gene expression data |
| `log2fc` | numeric | **required** | Fold change values |
| `sample_groups` | character | NULL | Sample annotations |
| `cluster_rows` | logical | TRUE | Cluster genes |
| `cluster_cols` | logical | TRUE | Cluster samples |
| `scale` | character | "row" | Normalization |
| `show_rownames` | logical | TRUE | Display gene names |
| `show_colnames` | logical | TRUE | Display sample names |
| `fc_threshold` | numeric | 1 | Significance cutoff |
| `width_ratio` | numeric | c(4,1) | Panel proportions |
| `title` | character | NULL | Plot title |

## Scale Options

- `"row"` - Normalize each gene (Z-score per gene) → **Recommended**
- `"column"` - Normalize each sample (Z-score per sample)
- `"none"` - No normalization (raw values)

## Color Schemes

### Diverging (for expression)
```r
# Blue-White-Red (default)
colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)

# Green-Black-Red
colorRampPalette(c("green", "black", "red"))(100)

# Purple-White-Orange
colorRampPalette(c("purple", "white", "orange"))(100)
```

### FC Colors
```r
# Default
fc_colors = c("#2166AC", "#B2182B")  # Blue, Red

# Colorblind-friendly
fc_colors = c("#0072B2", "#D55E00")  # Blue, Orange
```

## Integration Examples

### DESeq2
```r
library(DESeq2)
res <- results(dds, contrast = c("condition", "treated", "control"))
sig <- subset(res, padj < 0.05)

plot_heatmap_with_fc(
  expression_matrix = counts(dds, normalized = TRUE)[rownames(sig), ],
  log2fc = sig$log2FoldChange,
  sample_groups = dds$condition
)
```

### edgeR
```r
library(edgeR)
top <- topTags(et, n = 100)$table

plot_heatmap_with_fc(
  expression_matrix = cpm(dge, log = TRUE)[rownames(top), ],
  log2fc = top$logFC,
  sample_groups = groups
)
```

### limma
```r
library(limma)
top <- topTable(fit, number = 50)

plot_heatmap_with_fc(
  expression_matrix = exprs(eset)[rownames(top), ],
  log2fc = top$logFC,
  sample_groups = pData(eset)$condition
)
```

## Saving Figures

### For Publication (PDF)
```r
pdf("figure1.pdf", width = 10, height = 8)
plot_heatmap_with_fc(expr_data, log2fc, sample_groups = groups)
dev.off()
```

### For Presentation (PNG, high-res)
```r
png("figure1.png", width = 1200, height = 900, res = 150)
plot_heatmap_with_fc(expr_data, log2fc, sample_groups = groups)
dev.off()
```

### For Web (PNG, lower-res)
```r
png("figure1.png", width = 800, height = 600, res = 96)
plot_heatmap_with_fc(expr_data, log2fc, sample_groups = groups)
dev.off()
```

## Troubleshooting

| Problem | Solution |
|---------|----------|
| Gene names overlap | `show_rownames = FALSE` |
| Plot too compressed | Increase figure dimensions |
| Colors too light | Adjust color palette intensity |
| FC bars too small | `width_ratio = c(3, 1.5)` |
| No pattern visible | Try `scale = "row"` |

## Tips

💡 **Best practice**: Always use `scale = "row"` for gene expression

💡 **Gene labels**: Auto-hidden when >50 genes to prevent overlap

💡 **Threshold**: Set based on biology (commonly 1, 1.5, or 2)

💡 **Width ratio**: Adjust based on importance (4:1 default is balanced)

💡 **Clustering**: Reveals patterns; disable for ordered data (e.g., time-series)

## Common Mistakes

❌ **Wrong dimensions**: log2fc length ≠ number of genes
```r
# ERROR: 20 genes but 30 FC values
plot_heatmap_with_fc(matrix[1:20, ], log2fc[1:30])
```

❌ **Mismatched groups**: sample_groups length ≠ number of samples
```r
# ERROR: 10 samples but 8 group labels
plot_heatmap_with_fc(matrix[, 1:10], log2fc, sample_groups = groups[1:8])
```

✅ **Correct usage**:
```r
# Ensure dimensions match
nrow(expr_data) == length(log2fc)  # TRUE
ncol(expr_data) == length(sample_groups)  # TRUE
```

## See Also

- Full documentation: `?plot_heatmap_with_fc`
- Usage guide: `examples/heatmap_fc_usage.md`
- Test examples: `test_heatmap_fc.R`
