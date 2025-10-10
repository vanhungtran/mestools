# GEO Dataset Functions in mestools

This document describes the GEO (Gene Expression Omnibus) dataset reading functions added to the mestools package.

## Overview

The mestools package now includes comprehensive functions to download and process GEO datasets from NCBI's Gene Expression Omnibus database. These functions are designed to handle both individual datasets and batch processing of multiple datasets.

## Prerequisites

Before using the GEO functions, you need to install the GEOquery package from Bioconductor:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("GEOquery")
```

## Available Functions

### 1. `read_geo_dataset(gse_id, destdir, getGPL, AnnotGPL)`

Downloads and processes a single GEO dataset.

**Parameters:**
- `gse_id`: Character string specifying the GSE ID (e.g., "GSE102628")
- `destdir`: Directory to save downloaded files (default: tempdir())
- `getGPL`: Whether to download platform annotation data (default: TRUE)
- `AnnotGPL`: Whether to annotate expression data with gene symbols (default: TRUE)

**Returns:** A list containing:
- `gse_object`: The original GEOquery object
- `expression_matrix`: Processed expression matrix
- `phenotype_data`: Sample metadata
- `feature_data`: Probe/gene annotations
- `gse_id`: The GSE identifier
- `n_samples`: Number of samples
- `n_features`: Number of features
- `platform`: Platform information

**Example:**
```r
# Read a single dataset
result <- read_geo_dataset("GSE102628")
expression_data <- result$expression_matrix
sample_info <- result$phenotype_data
```

### 2. `read_multiple_geo_datasets(gse_ids, destdir, getGPL, AnnotGPL, max_parallel, sleep_between)`

Downloads and processes multiple GEO datasets in batch.

**Parameters:**
- `gse_ids`: Character vector of GSE IDs
- `destdir`: Directory to save downloaded files (default: tempdir())
- `getGPL`: Whether to download platform annotation data (default: TRUE)
- `AnnotGPL`: Whether to annotate expression data with gene symbols (default: TRUE)
- `max_parallel`: Maximum parallel downloads (default: 3)
- `sleep_between`: Sleep time between downloads in seconds (default: 2)

**Returns:** Named list where each element contains a processed dataset

**Example:**
```r
# Read multiple datasets
gse_list <- c("GSE102628", "GSE102641", "GSE102725")
results <- read_multiple_geo_datasets(gse_list)

# Access individual datasets
dataset1 <- results$GSE102628$expression_matrix
```

### 3. `get_geo_summary(gse_ids)`

Gets summary information for GEO datasets without downloading full expression data.

**Parameters:**
- `gse_ids`: Character vector of GSE IDs

**Returns:** Data.frame with summary information including title, platform, sample count, etc.

**Example:**
```r
# Get summary info
summary_info <- get_geo_summary(c("GSE102628", "GSE102641"))
print(summary_info)
```

### 4. `get_default_gse_list()`

Returns the predefined list of 127 GSE IDs for batch processing.

**Returns:** Character vector of GSE IDs

**Example:**
```r
# Get the complete list
all_gse <- get_default_gse_list()
length(all_gse)  # 127 datasets
```

### 5. `process_all_geo_datasets(destdir, getGPL, AnnotGPL, subset_size)`

Convenience function to download and process all datasets in the default GSE list.

**Parameters:**
- `destdir`: Directory to save files (default: "geo_data")
- `getGPL`: Whether to download platform annotation data (default: TRUE)
- `AnnotGPL`: Whether to annotate expression data (default: TRUE)
- `subset_size`: If specified, only process first N datasets (for testing)

**Returns:** Named list containing all processed datasets

**Example:**
```r
# Process first 5 datasets for testing
test_results <- process_all_geo_datasets(subset_size = 5)

# Process all datasets (this will take several hours!)
all_results <- process_all_geo_datasets()
```

## Default GSE Dataset List

The package includes a curated list of 127 GEO datasets:

```
GSE102628, GSE102641, GSE102725, GSE103489, GSE104509, GSE106087,
GSE106992, GSE107361, GSE107871, GSE109182, GSE109248, GSE111053,
GSE111054, GSE111055, GSE11307, GSE114286, GSE114729, GSE116486,
... (and 109 more)
```

## Usage Examples

### Basic Usage
```r
# Load the library (after installation)
library(mestools)

# Read a single dataset
result <- read_geo_dataset("GSE121212")
print(dim(result$expression_matrix))

# Get dataset summaries
summaries <- get_geo_summary(c("GSE121212", "GSE102628"))
print(summaries)
```

### Batch Processing
```r
# Process a small subset for testing
test_datasets <- c("GSE121212", "GSE102628", "GSE102641")
results <- read_multiple_geo_datasets(test_datasets, destdir = "my_geo_data")

# Check results
summary(attr(results, "summary"))

# Access individual results
if ("GSE121212" %in% names(results)) {
  expr_matrix <- results$GSE121212$expression_matrix
  sample_data <- results$GSE121212$phenotype_data
}
```

### Large-Scale Processing
```r
# Process all default datasets (WARNING: This will take several hours!)
all_results <- process_all_geo_datasets(
  destdir = "complete_geo_data",
  getGPL = TRUE,
  AnnotGPL = TRUE
)

# Save results for later use
saveRDS(all_results, "all_geo_datasets.rds")
```

## Error Handling

The functions include comprehensive error handling:
- Failed downloads are logged but don't stop batch processing
- Network timeouts are handled gracefully
- Missing or invalid GSE IDs are reported
- Summary information tracks successful vs. failed downloads

## Performance Considerations

- **Respectful downloading**: Functions include delays between downloads to avoid overwhelming NCBI servers
- **Large datasets**: Some datasets can be several GB; ensure sufficient disk space
- **Memory usage**: Large expression matrices require substantial RAM
- **Time requirements**: Processing all 127 datasets can take 2-4 hours depending on network speed

## Troubleshooting

### Common Issues:
1. **GEOquery not installed**: Install with `BiocManager::install("GEOquery")`
2. **Network timeouts**: Increase delay between downloads
3. **Disk space**: Check available space before large batch downloads
4. **Memory errors**: Process datasets in smaller batches

### Support:
If you encounter issues, check:
1. Your internet connection
2. NCBI GEO database status
3. Available disk space and memory
4. R session logs for specific error messages

## Testing

Run the provided test script to verify functionality:
```r
source("test_geo_functions.R")
```

This will test all major functions with small datasets and provide performance feedback.
