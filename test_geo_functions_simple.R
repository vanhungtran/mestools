# Simple test for GEO functions (without requiring GEOquery installation)
# Tests basic functionality and structure

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Load the mestools package functions
source("R/mestools-core.R")

message("=== Simple GEO Functions Test ===")

# Test 1: Get the default GSE list (this doesn't require GEOquery)
message("\n=== Test 1: Getting default GSE list ===")
tryCatch({
  gse_list <- get_default_gse_list()
  message("✅ SUCCESS: Total GSE datasets available: ", length(gse_list))
  message("   First 10 datasets: ", paste(head(gse_list, 10), collapse = ", "))
  message("   Last 10 datasets: ", paste(tail(gse_list, 10), collapse = ", "))
}, error = function(e) {
  message("❌ FAILED: ", e$message)
})

# Test 2: Check if functions are properly defined
message("\n=== Test 2: Checking function definitions ===")

functions_to_check <- c("read_geo_dataset", "read_multiple_geo_datasets", 
                       "get_geo_summary", "get_default_gse_list", 
                       "process_all_geo_datasets")

for (func_name in functions_to_check) {
  if (exists(func_name) && is.function(get(func_name))) {
    message("✅ ", func_name, " - Function defined correctly")
  } else {
    message("❌ ", func_name, " - Function not found or not a function")
  }
}

# Test 3: Test function parameters (without actually calling GEOquery)
message("\n=== Test 3: Testing function parameters ===")

# Check read_geo_dataset parameters
tryCatch({
  func_args <- formals(read_geo_dataset)
  expected_args <- c("gse_id", "destdir", "getGPL", "AnnotGPL")
  has_all_args <- all(expected_args %in% names(func_args))
  
  if (has_all_args) {
    message("✅ read_geo_dataset has all expected parameters")
  } else {
    message("❌ read_geo_dataset missing parameters")
  }
}, error = function(e) {
  message("❌ Error checking read_geo_dataset parameters: ", e$message)
})

# Test 4: Test that functions give proper error messages when GEOquery is not available
message("\n=== Test 4: Testing error handling without GEOquery ===")

# This should give a proper error message about GEOquery not being installed
tryCatch({
  result <- read_geo_dataset("GSE121212")
}, error = function(e) {
  if (grepl("GEOquery", e$message)) {
    message("✅ Proper error message for missing GEOquery: ", e$message)
  } else {
    message("❌ Unexpected error: ", e$message)
  }
})

message("\n=== Test Summary ===")
message("✅ Basic function structure tests completed")
message("📝 To run full tests with data download, install GEOquery first:")
message("   if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')")
message("   BiocManager::install('GEOquery')")
message("   Then run the full test_geo_functions.R script")

message("\n=== Test completed! ===")
