# Test script for GEO dataset functions
# Run this script to test the new GEO functionality

# Set CRAN mirror first
options(repos = c(CRAN = "https://cloud.r-project.org"))

# First install required packages if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("GEOquery", quietly = TRUE)) {
  BiocManager::install("GEOquery")
}

# Load the mestools package functions
source("R/mestools-core.R")

# Test 1: Get the default GSE list
message("=== Test 1: Getting default GSE list ===")
gse_list <- get_default_gse_list()
message("Total GSE datasets available: ", length(gse_list))
message("First 10 datasets: ", paste(head(gse_list, 10), collapse = ", "))

# Test 2: Read a single GEO dataset (small one for testing)
message("\n=== Test 2: Reading a single GEO dataset ===")
result <- read_geo_dataset("GSE121212")  # This should be a small dataset
if (!is.null(result)) {
  message("Successfully loaded GSE121212")
  message("Expression matrix dimensions: ", paste(dim(result$expression_matrix), collapse = " x "))
  message("Number of samples: ", result$n_samples)
  message("Number of features: ", result$n_features)
  message("Platform: ", result$platform)
}

# Test 3: Read multiple datasets (just 2 for testing)
message("\n=== Test 3: Reading multiple GEO datasets ===")
test_gse_list <- c("GSE121212", "GSE102628")  # Small test set
multi_results <- read_multiple_geo_datasets(test_gse_list, sleep_between = 1)

if (length(multi_results) > 0) {
  message("Successfully processed ", length(multi_results), " datasets")
  for (gse_id in names(multi_results)) {
    result <- multi_results[[gse_id]]
    message("  ", gse_id, ": ", result$n_features, " features x ", result$n_samples, " samples")
  }
}

# Test 4: Get summary information
message("\n=== Test 4: Getting dataset summaries ===")
summary_info <- get_geo_summary(c("GSE121212", "GSE102628"))
print(summary_info)

message("\n=== All tests completed! ===")
