# Performance and integration tests for mestools

test_that("batch_apply performance with large datasets", {
  skip_on_cran()
  
  # Test performance with larger datasets
  large_vector <- 1:10000
  
  start_time <- Sys.time()
  result <- batch_apply(large_vector, function(x) x^2, .progress = FALSE)
  end_time <- Sys.time()
  
  # Should complete in reasonable time (less than 5 seconds)
  expect_true(as.numeric(end_time - start_time, units = "secs") < 5)
  expect_equal(length(result), 10000)
  expect_equal(result[[1]], 1)
  expect_equal(result[[100]], 10000)
})

test_that("generate_random_df performance and memory", {
  skip_on_cran()
  
  # Test with moderately large data frame
  df <- generate_random_df(10000, 20)
  
  # Check structure
  expect_equal(nrow(df), 10000)
  expect_equal(ncol(df), 20)
  
  # Check memory usage is reasonable
  mem_size <- as.numeric(object.size(df))
  expect_true(mem_size > 0)
  
  # Check data quality
  summary_result <- quick_summary(df)
  expect_equal(summary_result$dimensions, c(10000, 20))
  expect_true(all(summary_result$missing_values == 0))
})

test_that("create_project_structure integration", {
  # Test complete workflow
  temp_base <- tempfile()
  project_name <- "integration_test_project"
  project_path <- file.path(temp_base, project_name)
  
  # Create project structure
  result <- create_project_structure(project_path)
  expect_true(result)
  
  # Verify all expected directories exist
  expected_dirs <- c("inputs", "scripts", "output", "docs", "ref", "reports")
  for (dir in expected_dirs) {
    full_path <- file.path(project_path, dir)
    expect_true(dir.exists(full_path))
  }
  
  # Test that we can write files to these directories
  test_file <- file.path(project_path, "inputs", "test.txt")
  writeLines("test content", test_file)
  expect_true(file.exists(test_file))
  
  # Test reading the file back
  content <- read_file_safe(test_file)
  expect_equal(content, "test content")
  
  # Clean up
  unlink(temp_base, recursive = TRUE)
})

test_that("file operations integration workflow", {
  temp_dir <- tempfile()
  dir.create(temp_dir)
  
  # Generate test data
  test_df <- generate_random_df(100, 5)
  
  # Save as CSV
  csv_path <- file.path(temp_dir, "test_data.csv")
  write.csv(test_df, csv_path, row.names = FALSE)
  
  # Read back and verify
  read_df <- read_file_safe(csv_path)
  expect_s3_class(read_df, "data.frame")
  expect_equal(nrow(read_df), 100)
  
  # Generate summary
  summary_result <- quick_summary(read_df)
  expect_equal(summary_result$dimensions[1], 100)
  
  # Save as RDS
  rds_path <- file.path(temp_dir, "test_data.rds")
  saveRDS(test_df, rds_path)
  
  # Read RDS back
  read_rds <- read_file_safe(rds_path)
  expect_identical(test_df, read_rds)
  
  # Clean up
  unlink(temp_dir, recursive = TRUE)
})

test_that("dependency checking integration", {
  # Test with known base packages
  base_packages <- c("base", "utils", "stats", "graphics", "grDevices")
  result <- check_dependencies(base_packages, auto_install = FALSE)
  expect_true(result)
  
  # Test workflow with batch processing
  package_lists <- list(
    essential = c("utils", "stats"),
    graphics = c("graphics", "grDevices"),
    fake = c("nonexistent123", "alsofake456")
  )
  
  results <- batch_apply(package_lists, function(pkgs) {
    check_dependencies(pkgs, auto_install = FALSE)
  }, .progress = FALSE)
  
  expect_true(results[[1]])  # essential packages should exist
  expect_true(results[[2]])  # graphics packages should exist  
  expect_false(results[[3]]) # fake packages should not exist
})

test_that("error handling across functions", {
  # Test cascading error handling
  temp_dir <- tempfile()
  # Deliberately not creating the directory
  
  # This should fail gracefully
  expect_error(create_project_structure(file.path(temp_dir, "nonexistent", "path")))
  
  # Test file operations with invalid paths
  expect_error(read_file_safe("/path/that/does/not/exist.csv"))
  
  # Test data operations with invalid input
  expect_error(quick_summary("not a data frame"))
  expect_error(quick_summary(NULL))
  expect_error(quick_summary(list(a = 1, b = 2)))
})

test_that("function parameter consistency", {
  # Test that functions handle similar parameters consistently
  
  # Test path parameters
  temp_dir <- tempdir()
  
  # Both functions should handle absolute paths
  expect_no_error(create_project_structure(temp_dir))
  
  # Test file operations
  test_csv <- file.path(temp_dir, "consistency_test.csv")
  write.csv(mtcars[1:5, ], test_csv, row.names = FALSE)
  
  expect_no_error(read_file_safe(test_csv))
  expect_no_error(quick_summary(read_file_safe(test_csv)))
  
  # Clean up
  unlink(test_csv)
})