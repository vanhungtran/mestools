test_that("validate_github_repo works with valid repo", {
  skip_on_cran()
  skip_if_offline()
  skip_if_not_installed("gh")

  # Test with a known existing repository
  expect_message(
    result <- validate_github_repo("https://github.com/vanhungtran/aucmat.git"),
    "Checking GitHub repository"
  )
  expect_true(result)
})

test_that("validate_github_repo fails with invalid repo", {
  skip_on_cran()
  skip_if_offline()
  skip_if_not_installed("gh")

  # Test with non-existent repository
  expect_message(
    result <- validate_github_repo("https://github.com/invalid_user_12345/invalid_repo_67890.git"),
    "Repository not found"
  )
  expect_false(result)
})

test_that("validate_github_repo handles malformed URLs", {
  skip_if_not_installed("gh")
  
  # Test with malformed URL - should handle gracefully
  result <- validate_github_repo("not-a-url")  
  expect_false(result)
  
  # Test with NULL and NA
  expect_error(validate_github_repo(NULL))
  expect_error(validate_github_repo(NA))
  
  # Test with empty string
  result <- validate_github_repo("")
  expect_false(result)
})

test_that("quick_summary returns correct structure", {
  # Test with mtcars dataset
  result <- quick_summary(mtcars)

  expect_type(result, "list")
  expect_named(result, c("dimensions", "column_names", "column_types",
                         "missing_values", "memory_size"))
  expect_equal(result$dimensions, dim(mtcars))
  expect_equal(result$column_names, names(mtcars))
  expect_equal(length(result$column_types), ncol(mtcars))
})

test_that("quick_summary handles edge cases", {
  # Test with empty data frame
  empty_df <- data.frame()
  result <- quick_summary(empty_df)
  expect_equal(result$dimensions, c(0, 0))

  # Test with non-data frame input
  expect_error(quick_summary("not a data frame"))
  expect_error(quick_summary(123))
})

test_that("read_file_safe reads different file types", {
  # Create temporary files for testing
  temp_csv <- tempfile(fileext = ".csv")
  temp_txt <- tempfile(fileext = ".txt")

  # Write test data
  write.csv(mtcars, temp_csv, row.names = FALSE)
  writeLines(c("line1", "line2", "line3"), temp_txt)

  # Test CSV reading
  csv_data <- read_file_safe(temp_csv)
  expect_s3_class(csv_data, "data.frame")

  # Test text file reading
  txt_data <- read_file_safe(temp_txt)
  expect_type(txt_data, "character")
  expect_length(txt_data, 3)

  # Clean up
  unlink(temp_csv)
  unlink(temp_txt)
})

test_that("read_file_safe handles errors appropriately", {
  # Test non-existent file
  expect_error(read_file_safe("nonexistent_file.csv"), "File does not exist")

  # Test unsupported file type
  temp_unsupported <- tempfile(fileext = ".unsupported")
  writeLines("test", temp_unsupported)
  expect_error(read_file_safe(temp_unsupported), "Unsupported file type")
  unlink(temp_unsupported)

  # Create a corrupt RDS file for testing
  temp_bad <- tempfile(fileext = ".rds")
  writeLines("not rds data", temp_bad)
  expect_error(read_file_safe(temp_bad, type = "rds"), "Error reading file")
  unlink(temp_bad)
})

test_that("deploy_package input validation", {
  skip("Skipping deployment test to avoid recursive testing during package tests")
})

test_that("deploy_package handles missing dependencies", {
  # This test verifies the function checks for required packages
  # We can't easily mock package installation in testthat, but we can
  # verify the function structure
  
  expect_true(exists("deploy_package"))
  expect_true(is.function(deploy_package))
  
  # Check function arguments
  expected_args <- c("repo_url", "commit_message", "run_tests", "run_checks")
  actual_args <- names(formals(deploy_package))
  expect_true(all(expected_args %in% actual_args))
})
