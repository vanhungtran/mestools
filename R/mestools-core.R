# Package: mestools
# Collection of utility functions for data analysis and package development

NULL

#' Validate GitHub Repository
#'
#' Checks if a GitHub repository exists and is accessible using the GitHub API.
#'
#' @param repo_url GitHub repository URL (e.g., "https://github.com/user/repo.git")
#' @return Logical indicating if repository exists
#' @export
#' @examples
#' \dontrun{
#' validate_github_repo("https://github.com/r-lib/usethis.git")
#' }
validate_github_repo <- function(repo_url) {
  if (!requireNamespace("gh", quietly = TRUE)) {
    stop("Please install the gh package: install.packages(\"gh\")")
  }

  # Extract owner and repo name from URL
  repo_parts <- strsplit(repo_url, "/")[[1]]
  owner <- repo_parts[length(repo_parts)-1]
  repo_name <- gsub("\\.git$", "", repo_parts[length(repo_parts)])

  message("Checking GitHub repository: ", repo_url)

  tryCatch({
    repo_info <- gh::gh(paste0("/repos/", owner, "/", repo_name))
    message("✅ Repository exists: ", repo_info$html_url)
    message("   Description: ", repo_info$description %||% "No description")
    message("   Visibility: ", ifelse(repo_info$private, "Private", "Public"))
    message("   Created: ", repo_info$created_at)
    return(TRUE)
  }, error = function(e) {
    if (grepl("404", e$message)) {
      message("❌ Repository not found: ", repo_url)
    } else {
      message("❌ Error checking repository: ", e$message)
    }
    return(FALSE)
  })
}

#' Deploy R Package to GitHub
#'
#' Complete workflow for rebuilding, testing, and deploying an R package to GitHub.
#' This function automates the entire deployment process including documentation
#' rebuilding, testing, and git operations.
#'
#' @param repo_url GitHub repository URL (optional, defaults to package name-based URL)
#' @param commit_message Commit message (optional, auto-generated if not provided)
#' @param run_tests Logical indicating whether to run tests (default: TRUE)
#' @param run_checks Logical indicating whether to run R CMD check (default: TRUE)
#' @return Logical indicating deployment success
#' @export
#' @examples
#' \dontrun{
#' # Deploy with default settings
#' deploy_package()
#'
#' # Deploy with custom settings
#' deploy_package(
#'   repo_url = "https://github.com/username/mypackage.git",
#'   commit_message = "Custom update message",
#'   run_tests = TRUE,
#'   run_checks = FALSE
#' )
#' }
deploy_package <- function(repo_url = NULL, commit_message = NULL,
                           run_tests = TRUE, run_checks = TRUE) {
  # Load required packages
  required_packages <- c("devtools", "usethis", "gert", "here")
  invisible(lapply(required_packages, function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg)
    }
  }))

  package_name <- basename(here::here())

  if (is.null(repo_url)) {
    repo_url <- paste0("https://github.com/vanhungtran/", package_name, ".git")
  }

  if (is.null(commit_message)) {
    commit_message <- paste("Daily update for", package_name, "-", Sys.Date())
  }

  tryCatch({
    message("🚀 Starting deployment for: ", package_name)
    message("📦 Target repository: ", repo_url)

    # Validate GitHub repo first
    if (!validate_github_repo(repo_url)) {
      stop("GitHub repository validation failed. Deployment aborted.")
    }

    # Package rebuilding steps
    message("1. 📝 Cleaning and rebuilding documentation...")
    unlink("man", recursive = TRUE)
    dir.create("man")
    devtools::document()

    message("2. 📚 Rebuilding vignettes...")
    devtools::build_vignettes()

    if (run_tests) {
      message("3. 🧪 Running tests...")
      test_results <- devtools::test()
      if (any(as.data.frame(test_results)$failed > 0)) {
        stop("Tests failed. Deployment aborted.")
      }
    } else {
      message("3. ⚠️  Skipping tests...")
    }

    if (run_checks) {
      message("4. ✅ Running package check...")
      check_results <- devtools::check(error_on = "never")
      if (length(check_results$errors) > 0) {
        stop("Package check failed with errors.")
      }
      if (length(check_results$warnings) > 0) {
        warning("Package check completed with warnings - review manually")
      }
    } else {
      message("4. ⚠️  Skipping package check...")
    }

    # Git operations
    message("5. 🔗 Setting up git remote...")
    usethis::use_git_remote("origin", repo_url, overwrite = TRUE)

    message("6. 💾 Committing changes...")

    # Stage all changes
    gert::git_add(".")

    # Check if there are changes to commit
    status <- gert::git_status()
    if (nrow(status) == 0) {
      message("No changes to commit.")
      return(invisible(FALSE))
    }

    # Commit changes
    gert::git_commit(commit_message)

    message("7. 🚀 Pushing to GitHub...")
    gert::git_push(remote = "origin", set_upstream = TRUE)

    message("🎉 Deployment completed successfully!")
    message("📊 View your package at: ", sub(".git$", "", repo_url))
    return(invisible(TRUE))

  }, error = function(e) {
    message("❌ Deployment failed: ", e$message)
    return(invisible(FALSE))
  })
}

#' Quick Data Summary
#'
#' Provides a comprehensive summary of a data frame including dimensions,
#' column types, missing values, and memory usage.
#'
#' @param data Data frame to summarize
#' @return List containing summary information
#' @export
#' @examples
#' data_summary <- quick_summary(mtcars)
#' print(data_summary$dimensions)
quick_summary <- function(data) {
  if (!is.data.frame(data)) {
    stop("Input must be a data frame")
  }

  list(
    dimensions = dim(data),
    column_names = names(data),
    column_types = sapply(data, class),
    missing_values = sapply(data, function(x) sum(is.na(x))),
    memory_size = format(utils::object.size(data), units = "auto")
  )
}

#' Safe File Reader
#'
#' Safely reads various file formats with comprehensive error handling
#' and automatic format detection.
#'
#' @param file_path Path to the file
#' @param type File type (auto-detected from extension if NULL)
#' @param ... Additional arguments passed to read functions
#' @return File contents
#' @export
#' @examples
#' \dontrun{
#' data <- read_file_safe("data.csv")
#' text <- read_file_safe("notes.txt")
#' }
read_file_safe <- function(file_path, type = NULL, ...) {
  if (!file.exists(file_path)) {
    stop("File does not exist: ", file_path)
  }

  if (is.null(type)) {
    type <- tools::file_ext(file_path)
  }

  tryCatch({
    switch(tolower(type),
           csv = read.csv(file_path, ...),
           tsv = read.delim(file_path, ...),
           rds = readRDS(file_path),
           txt = readLines(file_path, ...),
           stop("Unsupported file type: ", type)
    )
  }, error = function(e) {
    stop("Error reading file ", file_path, ": ", e$message)
  })
}