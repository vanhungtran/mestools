test_that("check_dependencies works correctly", {
  # Test with existing packages
  result <- check_dependencies(c("utils", "stats"), auto_install = FALSE)
  expect_true(result)

  # Test with non-existent package (should return FALSE without auto_install)
  result <- check_dependencies("nonexistentpackage123", auto_install = FALSE)
  expect_false(result)
})

test_that("create_project_structure creates directories", {
  # Create temporary directory for testing
  temp_dir <- tempfile()
  dir.create(temp_dir)

  # Test directory creation
  result <- create_project_structure(temp_dir)
  expect_true(result)

  # Check if directories were created
  expected_dirs <- c("inputs", "scripts", "output", "docs", "ref", "reports")
  for (dir in expected_dirs) {
    expect_true(dir.exists(file.path(temp_dir, dir)))
  }

  # Clean up
  unlink(temp_dir, recursive = TRUE)
})

test_that("create_project_structure handles existing directories", {
  temp_dir <- tempfile()
  dir.create(temp_dir)

  # Create one directory in advance
  dir.create(file.path(temp_dir, "inputs"))

  # Should still work without errors
  expect_silent(create_project_structure(temp_dir))

  unlink(temp_dir, recursive = TRUE)
})

test_that("batch_apply processes objects correctly", {
  # Test with simple function
  numbers <- 1:5
  results <- batch_apply(numbers, function(x) x^2, .progress = FALSE)

  expect_type(results, "list")
  expect_length(results, 5)
  expect_equal(results[[1]], 1)
  expect_equal(results[[3]], 9)

  # Test with function that throws errors
  results_with_errors <- batch_apply(
    c(1, "not_number", 3),
    function(x) as.numeric(x)^2,
    .progress = FALSE
  )
  expect_true(is.null(results_with_errors[[2]]))
})

test_that("generate_random_df creates correct structure", {
  # Test default parameters
  df <- generate_random_df()
  expect_equal(dim(df), c(100, 5))
  expect_s3_class(df, "data.frame")

  # Test custom parameters
  df_custom <- generate_random_df(n_rows = 10, n_cols = 3)
  expect_equal(dim(df_custom), c(10, 3))

  # Test reproducibility with seed
  df1 <- generate_random_df(5, 2, seed = 123)
  df2 <- generate_random_df(5, 2, seed = 123)
  expect_identical(df1, df2)

  # Test different seed produces different results
  df3 <- generate_random_df(5, 2, seed = 456)
  expect_false(identical(df1, df3))
})
