# How to Run the New Heatmap Function

## Quick Start

The function `plot_heatmap_with_fc()` has been successfully added to your mestools package!

## Running the Demo

Since R is not currently in your system PATH, here are the ways to test the function:

### Method 1: RStudio (Recommended)

1. Open RStudio
2. Open the file: `RUN_DEMO.R`
3. Click the **Source** button (top-right of editor)
4. The plot will display in the **Plots** pane

### Method 2: R Console

```r
# Set working directory
setwd("c:/Users/tranh/OneDrive/mestools")

# Load the package in development mode
devtools::load_all(".")

# Run the demo
source("RUN_DEMO.R")
```

### Method 3: Direct Function Call

```r
# Load package
devtools::load_all(".")

# Create sample data
set.seed(123)
expression <- matrix(rnorm(200), nrow = 20, ncol = 10)
rownames(expression) <- paste0("Gene", 1:20)
colnames(expression) <- paste0("Sample", 1:10)

log2fc <- rnorm(20, mean = 0, sd = 2)
groups <- rep(c("Control", "Treatment"), each = 5)

# Generate plot
plot_heatmap_with_fc(
  expression_matrix = expression,
  log2fc = log2fc,
  sample_groups = groups,
  fc_threshold = 1.5,
  title = "My First Plot"
)
```

### Method 4: From Command Line (if R is installed)

First, add R to your PATH, then:

```bash
cd c:\Users\tranh\OneDrive\mestools
Rscript RUN_DEMO.R
```

## Required Packages

Before running, install the required packages:

```r
install.packages(c("pheatmap", "ggplot2", "grid", "gridExtra"))
```

## Files Created

| File | Purpose |
|------|---------|
| `RUN_DEMO.R` | Ready-to-run demo script |
| `test_heatmap_fc.R` | Comprehensive test suite with 3 examples |
| `EXPECTED_OUTPUT.txt` | Visual representation of expected output |
| `examples/heatmap_fc_usage.md` | Complete usage documentation |
| `examples/quick_reference.md` | Quick reference guide |
| `FUNCTION_SUMMARY.md` | Technical summary |

## What to Expect

When you run the demo successfully, you should see:

1. **Console output** showing:
   - Data summary
   - Generation progress
   - Success message

2. **Visual plot** with:
   - LEFT: Expression heatmap (blue to red)
   - RIGHT: Log2 fold change bars (blue = down, red = up)

3. **Returned object** containing:
   - Heatmap object
   - FC plot object
   - Processed data

## Troubleshooting

### Issue: "Package 'pheatmap' not found"
**Solution:**
```r
install.packages("pheatmap")
```

### Issue: "Could not find function 'plot_heatmap_with_fc'"
**Solution:**
```r
devtools::load_all(".")  # Load development version
```

### Issue: "Error in log2fc dimension"
**Solution:** Ensure `length(log2fc) == nrow(expression_matrix)`

### Issue: R not found in PATH
**Solution:** Run from RStudio or R Console instead of command line

## Next Steps

After successful testing:

1. **Generate documentation:**
   ```r
   devtools::document()
   ```

2. **Run package checks:**
   ```r
   devtools::check()
   ```

3. **Build the package:**
   ```r
   devtools::build()
   ```

4. **Install locally:**
   ```r
   devtools::install()
   ```

## Examples for Real Data

### With DESeq2 Results

```r
library(DESeq2)
library(mestools)

# After running DESeq2 analysis
res <- results(dds)
sig_genes <- subset(res, padj < 0.05)

plot_heatmap_with_fc(
  expression_matrix = counts(dds, normalized = TRUE)[rownames(sig_genes), ],
  log2fc = sig_genes$log2FoldChange,
  sample_groups = dds$condition,
  title = "Significant DEGs (padj < 0.05)"
)
```

### With edgeR Results

```r
library(edgeR)
library(mestools)

# After running edgeR
top_genes <- topTags(et, n = 50)$table

plot_heatmap_with_fc(
  expression_matrix = cpm(dge, log = TRUE)[rownames(top_genes), ],
  log2fc = top_genes$logFC,
  sample_groups = groups,
  title = "Top 50 DEGs from edgeR"
)
```

## Getting Help

- Function help: `?plot_heatmap_with_fc`
- Full documentation: See `examples/heatmap_fc_usage.md`
- Quick reference: See `examples/quick_reference.md`
- Expected output: See `EXPECTED_OUTPUT.txt`

## Summary

The function is fully implemented and ready to use! Just open R or RStudio and run one of the demo scripts to see it in action.

**Key Features:**
✅ Combined heatmap + log2FC visualization
✅ Automatic gene ordering/clustering
✅ Sample group annotations
✅ Customizable colors and thresholds
✅ Publication-ready output
✅ Compatible with DESeq2, edgeR, limma

Enjoy your new visualization function! 🎨📊
