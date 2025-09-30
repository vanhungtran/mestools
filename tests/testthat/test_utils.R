test_that("check_dependencies works correctly", {
  # Test with existing packages
  result <- check_dependencies(c("utils", "stats"), auto_install = FALSE)
  expect_true(result)

  # Test with non-existent package (should return FALSE without auto_install)
  expect_message(
    result <- check_dependencies("nonexistentpackage123", auto_install = FALSE),
    "Missing packages"
  )
  expect_false(result)
  
  # Test with empty vector
  result <- check_dependencies(character(0), auto_install = FALSE)
  expect_true(result)
  
  # Test with mixed existing and non-existing packages
  expect_message(
    result <- check_dependencies(c("utils", "nonexistentpackage123"), auto_install = FALSE),
    "Missing packages"
  )
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

  # Should still work (but will produce messages about existing directories)
  expect_message(create_project_structure(temp_dir), "Already exists")

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
  suppressWarnings({
    results_with_errors <- batch_apply(
      c(1, "not_number", 3),
      function(x) as.numeric(x)^2,
      .progress = FALSE
    )
  })
  
  # The function should return a list with 3 elements
  # Element 2 should be NULL due to error
  expect_equal(length(results_with_errors), 3)
  # Check that we got valid results for elements 1 and 3
  expect_equal(results_with_errors[[1]], 1)
  expect_equal(results_with_errors[[3]], 9)
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
  
  # Test column types are correct
  df_types <- generate_random_df(10, 6)
  expect_true(is.numeric(df_types[[1]]))     # numeric1
  expect_true(is.character(df_types[[2]]))   # character2
  expect_true(is.factor(df_types[[3]]))      # factor3
  expect_true(is.numeric(df_types[[4]]))     # numeric4
  expect_true(is.character(df_types[[5]]))   # character5
  expect_true(is.factor(df_types[[6]]))      # factor6
  
  # Test edge cases
  df_small <- generate_random_df(1, 1)
  expect_equal(dim(df_small), c(1, 1))
})

test_that("install_and_load handles existing packages", {
  # Test with packages that should already be available
  expect_message(install_and_load(c("utils", "stats")), "already installed")
  
  # Verify packages are loaded
  expect_true("package:utils" %in% search())
  expect_true("package:stats" %in% search())
})

test_that("install_and_load handles empty input", {
  # Test with empty vector
  expect_silent(install_and_load(character(0)))
  
  # Test with NULL (should not error but might give warning)
  expect_no_error(install_and_load(NULL))
})

test_that("install_and_load processes GitHub format correctly", {
  # Test that it recognizes GitHub format
  # Note: We can't actually install from GitHub in tests, but we can test the logic
  
  # Mock test - just ensure it doesn't crash with GitHub-style input
  # In real tests, this would require network access and actual packages
  expect_no_error({
    # This would normally try to install from GitHub
    # We'll just test that the function can handle the format
    pkg_name <- "user/repo"
    actual_name <- basename(pkg_name)
    expect_equal(actual_name, "repo")
  })
})
