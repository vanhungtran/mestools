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






# ==============================================================================
# R Package Installation and Loading Function (with GitHub Search)
# ==============================================================================

#' Install and Load R Packages from Multiple Sources
#'
#' This function iterates through a vector of package names. If a package is
#' not installed, it attempts to install it from CRAN, then Bioconductor. If
#' still not found, it will search GitHub for a repository with that name and
#' prompt the user to interactively select from the top results for installation.
#'
#' @param packages A character vector or list of package names.
#'   To install directly from GitHub without searching, use the 'username/reponame' format.
#'
#' @return Loads the requested packages into the current session. Prints status
#'   messages and continues processing even if one package fails.
#'
#' @export
#' @examples
#' # packages <- c("ggplot2", "limma", "tidyverse/dplyr", "hrbrthemes")
#' # install_and_load(packages)
install_and_load <- function(packages) {

  # Loop through each package in the provided vector/list
  for (pkg in packages) {
    message(paste("\n--- Processing package:", pkg, "---"))

    actual_name <- basename(pkg)

    # 1. Check if the package is already installed and load it
    if (requireNamespace(actual_name, quietly = TRUE)) {
      message(paste("✅ Package '", actual_name, "' is already installed. Loading now.", sep = ""))
      library(actual_name, character.only = TRUE)
      next # Skip to the next package in the loop
    }

    message(paste("-> Package '", actual_name, "' not found. Attempting installation...", sep = ""))
    is_remote <- grepl("/", pkg)

    if (is_remote) {
      # Direct installation from GitHub/GitLab
      tryCatch({
        if (!requireNamespace("remotes", quietly = TRUE)) {
          message("--> Installing 'remotes' to handle remote installation...")
          install.packages("remotes", quiet = TRUE, repos = "https://cran.rstudio.com/")
        }
        message(paste("--> Installing '", pkg, "' from remote source...", sep = ""))
        remotes::install_github(pkg, quiet = TRUE)
      }, error = function(e) {
        warning(paste("--> Remote installation failed for '", pkg, "'.", sep = ""))
      })

    } else {
      # Installation from repositories (CRAN, Bioconductor, then search GitHub)
      # 2a. Try CRAN
      tryCatch({
        message(paste("--> Attempt 1/3: Installing '", pkg, "' from CRAN...", sep = ""))
        install.packages(pkg, quiet = TRUE, repos = "https://cran.rstudio.com/")
      }, error = function(e) {})

      # 2b. Try Bioconductor if not found on CRAN
      if (!requireNamespace(pkg, quietly = TRUE)) {
        tryCatch({
          message(paste("--> Attempt 2/3: Installing '", pkg, "' from Bioconductor...", sep = ""))
          if (!requireNamespace("BiocManager", quietly = TRUE)) {
            install.packages("BiocManager", quiet = TRUE, repos = "https://cran.rstudio.com/")
          }
          # Suppress updates to avoid lengthy builds during script run
          BiocManager::install(pkg, ask = FALSE, update = FALSE)
        }, error = function(e) {})
      }

      # 2c. Search GitHub if not found on CRAN or Bioconductor
      if (!requireNamespace(pkg, quietly = TRUE)) {
        message(paste("--> Attempt 3/3: Searching for '", pkg, "' on GitHub...", sep = ""))
        tryCatch({
          # Ensure 'gh' package is installed for API access
          if (!requireNamespace("gh", quietly = TRUE)) {
            message("--> Installing 'gh' package to search GitHub API...")
            install.packages("gh", quiet = TRUE, repos = "https://cran.rstudio.com/")
          }

          # Search GitHub repositories, sorted by stars
          search_query <- paste0(pkg, " in:name language:R")
          results <- gh::gh("/search/repositories", q = search_query, per_page = 5, sort = "stars", order = "desc")

          if (results$total_count > 0) {
            repos <- sapply(results$items, function(item) item$full_name)
            choice_msg <- paste("Found repositories for '", pkg, "'. Please choose one to install:", sep="")

            # Use menu() for interactive selection in the console
            choice <- utils::menu(c(repos, "None of the above"), title = choice_msg)

            if (choice > 0 && choice <= length(repos)) {
              chosen_repo <- repos[choice]
              message(paste("--> User selected:", chosen_repo))

              if (!requireNamespace("remotes", quietly = TRUE)) {
                install.packages("remotes", quiet = TRUE, repos = "https://cran.rstudio.com/")
              }
              remotes::install_github(chosen_repo, quiet = TRUE)
            } else {
              message("--> No package selected. Skipping GitHub installation.")
            }
          } else {
            message("--> No R packages found on GitHub matching that name.")
          }
        }, error = function(e) {
          warning(paste("--> GitHub search failed. Error:", e$message))
        })
      }
    }

    # --- Final Loading Attempt ---
    if (requireNamespace(actual_name, quietly = TRUE)) {
      message(paste("✅ Successfully installed and loaded '", actual_name, "'.", sep = ""))
      library(actual_name, character.only = TRUE)
    } else {
      warning(paste("❌ Failed to install package '", actual_name, "' from any source.", sep = ""))
    }
  } # End of for loop
}

# ==============================================================================
# Example Usage
# ==============================================================================
# To run these examples, source this file and then call the function in your console.
#
# source("install_and_load.R")
#
# # Example: Install a vector of packages, including one to be found on GitHub
# message("\n--- Testing a vector of packages ---")
#
 # packages_to_install <- c(
 #   "ggplot2",                                 # From CRAN
 #   "limma",                                   # From Bioconductor
 #   "tidyverse/dplyr",                         # Direct from GitHub
 #   "hrbrthemes",                              # Will trigger GitHub search
 #   "exampleRPackage"        # Non-existent, will fail gracefully
 # )

 # install_and_load(packages_to_install)

