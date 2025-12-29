# New Function: plot_heatmap_with_fc()

## Summary

Added a new visualization function to the mestools package that creates a combined plot showing:
- **Heatmap** of gene expression (left panel)
- **Log2 fold change** bar plot (right panel)

## Location

- **Function**: [R/mestools-core.R](R/mestools-core.R) (lines 432-628)
- **Test script**: [test_heatmap_fc.R](test_heatmap_fc.R)
- **Documentation**: [examples/heatmap_fc_usage.md](examples/heatmap_fc_usage.md)

## Features

### Main Capabilities
✅ Hierarchical clustering of genes and samples
✅ Row and column scaling options
✅ Automatic gene name handling (hides if >50 genes)
✅ Sample group annotation (colored columns)
✅ Customizable color schemes
✅ Fold change threshold visualization
✅ Adjustable panel width ratios
✅ Synchronized ordering between heatmap and FC plot

### Technical Details

**Required inputs:**
- `expression_matrix`: Numeric matrix (genes × samples)
- `log2fc`: Numeric vector of log2 fold change values

**Key parameters:**
- `sample_groups`: Sample annotations
- `cluster_rows/cluster_cols`: Enable/disable clustering
- `scale`: Row, column, or no scaling
- `fc_threshold`: Significance threshold for FC
- `heatmap_colors`: Color palette for heatmap
- `fc_colors`: Colors for up/down regulation
- `width_ratio`: Panel width proportions

## Dependencies Added

Updated [DESCRIPTION](DESCRIPTION) file with:
```r
Suggests:
    pheatmap,
    ggplot2,
    grid,
    gridExtra
```

## Usage Example

```r
library(mestools)

# Create expression data
expression <- matrix(rnorm(200), nrow = 20, ncol = 10)
rownames(expression) <- paste0("Gene", 1:20)
colnames(expression) <- paste0("Sample", 1:10)

# Create log2 fold change values
log2fc <- rnorm(20, mean = 0, sd = 2)

# Define groups
groups <- rep(c("Control", "Treatment"), each = 5)

# Generate plot
plot_heatmap_with_fc(
  expression_matrix = expression,
  log2fc = log2fc,
  sample_groups = groups,
  fc_threshold = 1.5,
  title = "Differential Gene Expression"
)
```

## Integration with DE Tools

Works seamlessly with:
- **DESeq2**: Use normalized counts + log2FoldChange
- **edgeR**: Use CPM values + logFC
- **limma**: Use normalized expression + logFC

See [examples/heatmap_fc_usage.md](examples/heatmap_fc_usage.md) for detailed examples.

## Return Value

Returns invisible list containing:
```r
$heatmap            # pheatmap object
$fc_plot            # ggplot2 object
$fc_data            # fold change data frame
$row_order          # gene ordering after clustering
$expression_matrix  # input matrix
$log2fc            # input fold changes
```

## Files Modified

1. ✅ `R/mestools-core.R` - Added function
2. ✅ `DESCRIPTION` - Added dependencies
3. ✅ `test_heatmap_fc.R` - Created test script
4. ✅ `examples/heatmap_fc_usage.md` - Created documentation

## Next Steps

To complete integration:

1. **Generate documentation**:
   ```r
   devtools::document()
   ```

2. **Test the function**:
   ```r
   source("test_heatmap_fc.R")
   ```

3. **Build package**:
   ```r
   devtools::build()
   devtools::check()
   ```

4. **Update NEWS.md** (if exists):
   ```markdown
   # mestools 0.1.0

   ## New Features

   * Added `plot_heatmap_with_fc()` function for combined heatmap and
     log2 fold change visualization
   ```

## Design Decisions

1. **Clustering synchronization**: Log2FC plot automatically matches heatmap gene order
2. **Auto-hide labels**: Prevents overlap when >50 genes
3. **Default colors**: Blue-white-red for expression, blue/red for FC (standard convention)
4. **Width ratio**: Default 4:1 provides good balance
5. **Silent heatmap**: Uses `silent=TRUE` to prevent duplicate display
6. **Grid system**: Uses grid/gridExtra for precise layout control

## Future Enhancements

Potential additions:
- [ ] Add p-value annotations
- [ ] Support for multiple fold change comparisons
- [ ] Interactive plotly version
- [ ] Export to specific figure formats
- [ ] Gene set highlighting
- [ ] Pathway enrichment overlay
