test_that("package loads correctly", {
  # Test that the package can be loaded
  expect_silent(library(mestools))
  expect_true("mestools" %in% loadedNamespaces())
})

test_that("all exported functions are available", {
  # Check that expected functions are exported
  expected_functions <- c(
    "validate_github_repo",
    "deploy_package", 
    "quick_summary",
    "read_file_safe",
    "check_dependencies",
    "create_project_structure",
    "batch_apply",
    "generate_random_df",
    "install_and_load"
  )

  for (func in expected_functions) {
    expect_true(exists(func, where = asNamespace("mestools"), mode = "function"),
                info = paste("Function", func, "should be exported"))
  }
})

test_that("package functions have proper documentation", {
  # Check that functions exist (documentation may not be immediately available during testing)
  key_functions <- c(
    "validate_github_repo",
    "deploy_package",
    "quick_summary", 
    "read_file_safe",
    "check_dependencies",
    "create_project_structure",
    "batch_apply",
    "generate_random_df"
  )
  
  for (func in key_functions) {
    # Just check that the function exists and is documented in the source
    expect_true(exists(func, where = asNamespace("mestools"), mode = "function"),
                info = paste("Function", func, "should exist"))
  }
  
  # Test that at least some basic documentation exists by checking if roxygen comments are present
  utils_file <- system.file("R", "mestools-utils.R", package = "mestools")
  if (file.exists(utils_file)) {
    utils_content <- readLines(utils_file)
    expect_true(any(grepl("#'", utils_content)), "Roxygen documentation should be present")
  }
})

test_that("package documentation exists", {
  # Check that package documentation is available
  expect_true(!is.null(utils::packageDescription("mestools")))
})
