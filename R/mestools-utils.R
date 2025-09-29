# Utility functions for mestools package

#' Check Package Dependencies
#'
#' Verifies if all required packages are installed and available.
#'
#' @param packages Character vector of package names
#' @param auto_install Logical, if TRUE missing packages will be installed automatically
#' @return Logical indicating if all packages are available
#' @export
#' @examples
#' check_dependencies(c("ggplot2", "dplyr"))
check_dependencies <- function(packages, auto_install = FALSE) {
  missing_pkgs <- character(0)

  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      missing_pkgs <- c(missing_pkgs, pkg)
    }
  }

  if (length(missing_pkgs) > 0) {
    message("Missing packages: ", paste(missing_pkgs, collapse = ", "))

    if (auto_install) {
      message("Installing missing packages...")
      utils::install.packages(missing_pkgs)
      return(TRUE)
    } else {
      message("Use auto_install = TRUE to install missing packages automatically")
      return(FALSE)
    }
  }

  message("All packages are available!")
  return(TRUE)
}

#' Create Project Directory Structure
#'
#' Creates a standardized directory structure for data analysis projects
#' with common folders for organization.
#'
#' @param path Base path for project (default: current directory)
#' @param directories Character vector of directory names to create
#' @return Invisible TRUE, creates directories as side effect
#' @export
#' @examples
#' \dontrun{
#' create_project_structure("my_analysis_project")
#' }
create_project_structure <- function(path = ".",
                                     directories = c("inputs",
                                                     "scripts",
                                                     "output",
                                                     "docs",
                                                     "ref",
                                                     "reports")) {

  created_dirs <- character(0)

  for (dir in directories) {
    dir_path <- file.path(path, dir)
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
      created_dirs <- c(created_dirs, dir_path)
      message("✅ Created: ", dir_path)
    } else {
      message("⚠️  Already exists: ", dir_path)
    }
  }

  if (length(created_dirs) > 0) {
    message("🎉 Project structure created successfully!")
    message("📁 Created directories: ", length(created_dirs))
  } else {
    message("📁 All directories already exist.")
  }

  invisible(TRUE)
}

#' Batch Apply Function with Progress Reporting
#'
#' Applies a function to multiple objects with progress reporting and
#' comprehensive error handling.
#'
#' @param objects List or vector of objects to process
#' @param fun Function to apply to each object
#' @param ... Additional arguments passed to fun
#' @param .progress Logical, whether to show progress (default: TRUE)
#' @return List of results, with NULL for failed operations
#' @export
#' @examples
#' results <- batch_apply(1:5, function(x) x^2)
batch_apply <- function(objects, fun, ..., .progress = TRUE) {
  results <- list()
  errors <- list()

  if (.progress) {
    message("Processing ", length(objects), " items...")
    pb <- utils::txtProgressBar(min = 0, max = length(objects), style = 3)
  }

  for (i in seq_along(objects)) {
    tryCatch({
      results[[i]] <- fun(objects[[i]], ...)
    }, error = function(e) {
      results[[i]] <- NULL
      errors[[i]] <- e$message
      warning("Error processing item ", i, ": ", e$message)
    })

    if (.progress) {
      utils::setTxtProgressBar(pb, i)
    }
  }

  if (.progress) {
    close(pb)
  }

  if (length(errors) > 0) {
    message("Completed with ", length(errors), " errors")
  } else {
    message("All items processed successfully!")
  }

  return(results)
}

#' Generate Random Data Frame
#'
#' Creates a random data frame for testing and demonstration purposes.
#'
#' @param n_rows Number of rows (default: 100)
#' @param n_cols Number of columns (default: 5)
#' @param seed Random seed for reproducibility (default: 123)
#' @return Data frame with random data
#' @export
#' @examples
#' df <- generate_random_df(50, 3)
generate_random_df <- function(n_rows = 100, n_cols = 5, seed = 123) {
  set.seed(seed)

  data_list <- list()
  for (i in 1:n_cols) {
    if (i %% 3 == 0) {
      data_list[[paste0("factor", i)]] <- factor(sample(letters[1:5], n_rows, replace = TRUE))
    } else if (i %% 3 == 1) {
      data_list[[paste0("numeric", i)]] <- rnorm(n_rows)
    } else {
      data_list[[paste0("character", i)]] <- replicate(n_rows, paste(sample(letters, 5), collapse = ""))
    }
  }

  as.data.frame(data_list)
}
