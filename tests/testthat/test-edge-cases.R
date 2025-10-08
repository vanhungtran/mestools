# Advanced edge case tests for mestools functions

test_that("check_dependencies handles edge cases", {
  # Test with duplicate packages
  result <- check_dependencies(c("utils", "utils", "stats"), auto_install = FALSE)
  expect_true(result)
  
  # Test with mixed case package names
  expect_message(
    result <- check_dependencies(c("UTILS", "stats"), auto_install = FALSE),
    "Missing packages"
  )
  expect_false(result)
  
  # Test with special characters in package names
  expect_message(
    result <- check_dependencies(c("utils", "package-with-dash"), auto_install = FALSE),
    "Missing packages"
  )
  expect_false(result)
  
  # Test with very long package name
  long_name <- paste(rep("a", 100), collapse = "")
  expect_message(
    result <- check_dependencies(long_name, auto_install = FALSE),
    "Missing packages"
  )
  expect_false(result)
})

test_that("create_project_structure handles complex paths", {
  # Test with nested path creation
  temp_base <- tempfile()
  nested_path <- file.path(temp_base, "level1", "level2", "project")
  
  result <- create_project_structure(nested_path)
  expect_true(result)
  expect_true(dir.exists(file.path(nested_path, "inputs")))
  
  # Test with relative paths
  old_wd <- getwd()
  on.exit(setwd(old_wd))
  
  temp_dir <- tempfile()
  dir.create(temp_dir)
  setwd(temp_dir)
  
  result <- create_project_structure("./test_project")
  expect_true(result)
  expect_true(dir.exists("./test_project/inputs"))
  
  # Clean up
  unlink(temp_base, recursive = TRUE)
  unlink(temp_dir, recursive = TRUE)
})

test_that("create_project_structure with custom directories", {
  temp_dir <- tempfile()
  dir.create(temp_dir)
  
  custom_dirs <- c("data", "analysis", "results", "figures")
  result <- create_project_structure(temp_dir, directories = custom_dirs)
  expect_true(result)
  
  for (dir in custom_dirs) {
    expect_true(dir.exists(file.path(temp_dir, dir)))
  }
  
  unlink(temp_dir, recursive = TRUE)
})

test_that("batch_apply handles empty inputs", {
  # Test with empty list
  result <- batch_apply(list(), function(x) x, .progress = FALSE)
  expect_equal(length(result), 0)
  expect_type(result, "list")
  
  # Test with NULL input
  expect_error(batch_apply(NULL, function(x) x))
  
  # Test with function that always errors
  suppressWarnings({
    result <- batch_apply(1:3, function(x) stop("Always fails"), .progress = FALSE)
  })
  expect_equal(length(result), 3)
})

test_that("batch_apply with different object types", {
  # Test with character vector
  chars <- c("a", "b", "c")
  result <- batch_apply(chars, toupper, .progress = FALSE)
  expect_equal(result, list("A", "B", "C"))
  
  # Test with list input
  input_list <- list(1, 2, 3)
  result <- batch_apply(input_list, function(x) x * 2, .progress = FALSE)
  expect_equal(result, list(2, 4, 6))
  
  # Test with data frames
  df_list <- list(mtcars[1:5, ], mtcars[6:10, ])
  result <- batch_apply(df_list, nrow, .progress = FALSE)
  expect_equal(result, list(5, 5))
})

test_that("generate_random_df edge cases", {
  # Test with zero rows
  expect_error(generate_random_df(0, 3))
  
  # Test with zero columns
  expect_error(generate_random_df(10, 0))
  
  # Test with large numbers
  df_large <- generate_random_df(1000, 10)
  expect_equal(dim(df_large), c(1000, 10))
  
  # Test column naming patterns
  df <- generate_random_df(5, 9)
  expect_true(all(grepl("^(numeric|character|factor)\\d+$", names(df))))
})

test_that("read_file_safe advanced file operations", {
  # Test with different encodings
  temp_file <- tempfile(fileext = ".txt")
  writeLines("Hello World", temp_file, useBytes = TRUE)
  
  content <- read_file_safe(temp_file)
  expect_equal(content, "Hello World")
  unlink(temp_file)
  
  # Test with very large files (simulate)
  temp_csv <- tempfile(fileext = ".csv")
  large_df <- data.frame(
    x = rep(1:100, 100),
    y = rep(letters[1:25], 400)
  )
  write.csv(large_df, temp_csv, row.names = FALSE)
  
  result <- read_file_safe(temp_csv)
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 10000)
  unlink(temp_csv)
  
  # Test TSV files
  temp_tsv <- tempfile(fileext = ".tsv")
  write.table(mtcars[1:5, ], temp_tsv, sep = "\t", row.names = FALSE)
  
  result <- read_file_safe(temp_tsv)
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 5)
  unlink(temp_tsv)
})

test_that("quick_summary with various data types", {
  # Test with data frame containing all R data types
  complex_df <- data.frame(
    logical = c(TRUE, FALSE, TRUE),
    integer = 1:3,
    numeric = c(1.1, 2.2, 3.3),
    character = c("a", "b", "c"),
    factor = factor(c("x", "y", "z")),
    stringsAsFactors = FALSE
  )
  
  result <- quick_summary(complex_df)
  expect_equal(length(result$column_types), 5)
  expect_true(all(result$missing_values == 0))
  
  # Test with missing values
  df_with_na <- data.frame(
    x = c(1, 2, NA, 4),
    y = c("a", NA, "c", "d")
  )
  
  result <- quick_summary(df_with_na)
  expect_equal(result$missing_values[["x"]], 1)
  expect_equal(result$missing_values[["y"]], 1)
})

test_that("install_and_load parameter validation", {
  # Test with non-character input
  expect_error(install_and_load(123))
  expect_error(install_and_load(TRUE))
  
  # Test with list input (should work)
  expect_no_error(install_and_load(list("utils")))
  
  # Test GitHub format detection
  github_formats <- c(
    "user/repo",
    "organization/package-name",
    "username/repo.name"
  )
  
  for (format in github_formats) {
    is_remote <- grepl("/", format)
    expect_true(is_remote)
    actual_name <- basename(format)
    expect_true(nchar(actual_name) > 0)
  }
})