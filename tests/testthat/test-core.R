test_that("validate_github_repo works with valid repo", {
  skip_on_cran()
  skip_if_offline()

  # Test with a known existing repository
  result <- validate_github_repo("https://github.com/vanhungtran/aucmat.git")
  expect_true(result)
})

test_that("validate_github_repo fails with invalid repo", {
  skip_on_cran()
  skip_if_offline()

  # Test with non-existent repository
  result <- validate_github_repo("https://github.com/invalid_user/invalid_repo.git")
  expect_false(result)
})

test_that("validate_github_repo handles malformed URLs", {
  # Test with malformed URL
  expect_error(validate_github_repo("not-a-url"))
  expect_error(validate_github_repo(NULL))
  expect_error(validate_github_repo(NA))
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
  expect_error(read_file_safe("nonexistent_file.csv"))

  # Test unsupported file type
  expect_error(read_file_safe("test.unsupported"))

  # Create a corrupt RDS file for testing
  temp_bad <- tempfile(fileext = ".rds")
  writeLines("not rds data", temp_bad)
  expect_error(read_file_safe(temp_bad, type = "rds"))
  unlink(temp_bad)
})
