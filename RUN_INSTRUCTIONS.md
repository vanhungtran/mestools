# How to Run the plot_heatmap_with_fc() Function

## Quick Start (Choose ONE method)

### Method 1: RStudio (Easiest) ⭐

1. **Open RStudio**
2. **Open the file**: `test_plot_now.R`
3. **Click "Source"** button (top-right of editor)
4. **View the plot** in the Plots pane (bottom-right)

### Method 2: R Console

```r
# Copy and paste this into R console:
setwd("c:/Users/tranh/OneDrive/mestools")
source("test_plot_now.R")
```

### Method 3: Direct Commands

```r
# Load package
devtools::load_all("c:/Users/tranh/OneDrive/mestools")

# Create sample data
set.seed(123)
expression <- matrix(rnorm(160, 10, 3), 20, 8)
expression[1:10, 5:8] <- expression[1:10, 5:8] + 4
expression[15:20, 5:8] <- expression[15:20, 5:8] - 3

rownames(expression) <- c("BRCA1", "TP53", "EGFR", "MYC", "KRAS",
                          "PTEN", "PIK3CA", "AKT1", "BRAF", "ERBB2",
                          "Gene11", "Gene12", "Gene13", "Gene14",
                          "VEGFA", "CDK4", "CDKN2A", "RB1", "NRAS", "FGFR1")
colnames(expression) <- c(paste0("Ctrl_", 1:4), paste0("Treat_", 1:4))

log2fc <- c(rnorm(10, 2.5, 0.4), rnorm(4, 0.1, 0.2), rnorm(6, -2, 0.3))
groups <- c(rep("Control", 4), rep("Treatment", 4))

# Generate plot
plot_heatmap_with_fc(
  expression_matrix = expression,
  log2fc = log2fc,
  sample_groups = groups,
  fc_threshold = 1.5,
  title = "Differential Gene Expression",
  show_rownames = TRUE
)
```

## What You'll See

```
┌────────────────────────────────────────────┬──────────────────┐
│     HEATMAP (Gene Expression)              │   LOG2FC BARS    │
│                                            │                  │
│         Control    Treatment               │                  │
│        ┌────┬────┬────┬────┐               │   -3  -1  1  3   │
│        │Ctrl│Ctrl│Trt │Trt │               │    │   │  │  │   │
│ BRCA1  │ ▒▒ │ ▒▒ │ ██ │ ██ │               │    BRCA1 ████▸  │
│ TP53   │ ▒▒ │ ▒▒ │ ██ │ ██ │               │    TP53  ████▸  │
│ EGFR   │ ▒▒ │ ▒▒ │ ▓▓ │ ▓▓ │               │    EGFR  ███▸   │
│ MYC    │ ▒▒ │ ▒▒ │ ▓▓ │ ▓▓ │               │    MYC   ███▸   │
│ Gene11 │ ▒▒ │ ▒▒ │ ▒▒ │ ▒▒ │               │    Gene11 ▌     │
│ Gene12 │ ▒▒ │ ▒▒ │ ▒▒ │ ▒▒ │               │    Gene12 ▌     │
│ VEGFA  │ ██ │ ██ │ ░░ │ ░░ │               │    VEGFA ◂████  │
│ CDK4   │ ▓▓ │ ▓▓ │ ░░ │ ░░ │               │    CDK4  ◂███   │
│ RB1    │ ▓▓ │ ▓▓ │ ░░ │ ░░ │               │    RB1   ◂███   │
└────────────────────────────────────────────┴──────────────────┘

LEFT: Heatmap colors (blue→white→red) show expression levels
RIGHT: Bars show fold change (blue=down, red=up)
BOTH: Gene names are synchronized and in the same order!
```

## Expected Output

When you run the script, you'll see:

```
===========================================
Generating combined heatmap + log2FC plot
===========================================

Plot generated successfully!
===========================================

What you should see:
  LEFT PANEL: Heatmap with gene names
    - Blue = low expression
    - Red = high expression
    - Genes clustered by similarity

  RIGHT PANEL: Log2FC bars with gene names
    - Red bars = upregulated genes
    - Blue bars = downregulated genes
    - Gene names match left panel

Data Summary:
  Genes: 20
  Samples: 8
  Groups: Control Treatment
  log2FC range: -2.48 2.99

First 5 genes (in clustered order):
       gene   log2fc significant
1    Gene_1  2.87947        TRUE
2    Gene_2  2.78610        TRUE
3    Gene_3  2.64434        TRUE
...
```

## Saving the Plot

### Save as PDF (for publication)
```r
pdf("my_plot.pdf", width = 10, height = 8)
plot_heatmap_with_fc(expression, log2fc, sample_groups = groups)
dev.off()
```

### Save as PNG (for presentation)
```r
png("my_plot.png", width = 1200, height = 900, res = 150)
plot_heatmap_with_fc(expression, log2fc, sample_groups = groups)
dev.off()
```

### Save as TIFF (for journal)
```r
tiff("my_plot.tiff", width = 10, height = 8, units = "in", res = 300)
plot_heatmap_with_fc(expression, log2fc, sample_groups = groups)
dev.off()
```

## Troubleshooting

### Error: "Package 'pheatmap' not found"
```r
install.packages("pheatmap")
```

### Error: "Could not find function plot_heatmap_with_fc"
```r
devtools::load_all(".")
```

### Error: "Dimension mismatch"
Make sure:
```r
length(log2fc) == nrow(expression_matrix)  # Must be TRUE
```

## Files Available

- `test_plot_now.R` - Ready to run script (this file)
- `RUN_DEMO.R` - More detailed demo with 3 examples
- `test_heatmap_fc.R` - Comprehensive test suite
- `examples/heatmap_fc_usage.md` - Full documentation
- `examples/quick_reference.md` - Quick reference guide

## Next Steps

After viewing the plot:
1. Try adjusting parameters (colors, threshold, etc.)
2. Use with your own data
3. Save the plot in your preferred format

Happy plotting! 📊✨
