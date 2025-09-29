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
    "generate_random_df"
  )

  for (func in expected_functions) {
    expect_true(exists(func, where = asNamespace("mestools"), mode = "function"),
                info = paste("Function", func, "should be exported"))
  }
})

test_that("package documentation exists", {
  # Check that package documentation is available
  expect_true(!is.null(utils::packageDescription("mestools")))
})
