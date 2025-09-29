# Setup code that runs before tests
library(testthat)
library(mestools)

# Set test-specific options
options(mestools.test_mode = TRUE)

# Create a temporary directory for file operations
test_temp_dir <- tempfile()
dir.create(test_temp_dir)

# Cleanup function that runs after tests
withr::defer({
  unlink(test_temp_dir, recursive = TRUE)
}, teardown_env())
