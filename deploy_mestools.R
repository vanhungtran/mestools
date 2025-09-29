#' Deploy mestools Package to GitHub
#'
#' Complete workflow for deploying mestools package to GitHub, handling both
#' initial deployment and subsequent updates with modifications.
#'
#' @param commit_message Custom commit message (optional)
#' @param force_initial Force reinitialization of Git repository
#' @param run_checks Whether to run devtools::check() (default: TRUE)
#' @param skip_document Whether to skip documentation rebuild (default: FALSE)
#' @return Invisible TRUE if successful
#' @export
deploy_mestools <- function(commit_message = NULL,
                            force_initial = FALSE,
                            run_checks = TRUE,
                            skip_document = FALSE) {

  # Load required packages
  required_packages <- c("usethis", "gert", "devtools", "here")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg)
      library(pkg, character.only = TRUE)
    } else {
      library(pkg, character.only = TRUE)
    }
  }

  # Set working directory to mestools package
  pkg_dir <- file.path(dirname(here::here()), "mestools")
  if (!dir.exists(pkg_dir)) {
    stop("Package directory not found: ", pkg_dir)
  }
  setwd(pkg_dir)

  message("🚀 Starting mestools package deployment...")
  message("📁 Package directory: ", pkg_dir)

  # Step 1: Update package documentation
  if (!skip_document) {
    message("1. 📝 Updating package documentation...")
    devtools::document()
  } else {
    message("1. ⚠️  Skipping documentation update...")
  }

  # Step 2: Run package checks
  if (run_checks) {
    message("2. ✅ Running package checks...")
    check_results <- devtools::check(error_on = "never")
    if (length(check_results$errors) > 0) {
      warning("Package check completed with errors. Review before deployment.")
    } else if (length(check_results$warnings) > 0) {
      message("Package check completed with warnings - review manually")
    } else {
      message("✅ Package check passed!")
    }
  } else {
    message("2. ⚠️  Skipping package checks...")
  }

  # Step 3: Git repository setup
  git_initialized <- dir.exists(".git")

  if (force_initial || !git_initialized) {
    message("3. 🔧 Setting up Git repository (initial deployment)...")

    # Remove existing .git if force_initial is TRUE
    if (force_initial && git_initialized) {
      message("   🧹 Removing existing Git repository...")
      unlink(".git", recursive = TRUE, force = TRUE)
    }

    # Initialize new Git repository
    usethis::use_git()

    # Configure Git user
    gert::git_config_set("user.name", "Van Hung Tran")
    gert::git_config_set("user.email", "tranhungydhcm@gmail.com")

    # Set up .gitignore
    usethis::use_git_ignore(c(
      ".Rhistory", ".RData", ".Rproj.user",
      "creat_pkg_name.R", "deploy_package.R", "Rfunctions_and_test.R",
      ".DS_Store", "*.Rproj", "*.tar.gz", "*.pdf"
    ))

  } else {
    message("3. 🔧 Git repository already initialized")
  }

  # Step 4: GitHub remote setup
  remote_list <- gert::git_remote_list()
  has_remote <- nrow(remote_list) > 0

  if (!has_remote) {
    message("4. 🐙 Setting up GitHub remote...")
    tryCatch({
      usethis::use_github()
    }, error = function(e) {
      message("   ⚠️  use_github() failed, setting remote manually...")
      gert::git_remote_add(
        name = "origin",
        url = "https://github.com/vanhungtran/mestools.git"
      )
    })
  } else {
    message("4. 🐙 GitHub remote already configured: ", remote_list$url[1])
  }

  # Step 5: Commit changes
  message("5. 💾 Committing changes...")

  # Check git status
  status <- gert::git_status()
  if (nrow(status) == 0) {
    message("   📭 No changes to commit")
    message("🎉 Deployment completed (no changes to push)")
    return(invisible(TRUE))
  }

  # Add all files
  gert::git_add(".")

  # Generate commit message if not provided
  if (is.null(commit_message)) {
    if (force_initial || !git_initialized) {
      commit_message <- "Initial commit: mestools package with utility functions for data analysis and package development"
    } else {
      commit_message <- paste(
        "Update:", format(Sys.time(), "%Y-%m-%d %H:%M"),
        "- Package improvements and documentation updates"
      )
    }
  }

  # Commit changes
  gert::git_commit(commit_message)
  message("   ✅ Committed: ", commit_message)

  # Step 6: Push to GitHub
  message("6. 🚀 Pushing to GitHub...")
  gert::git_push(remote = "origin", set_upstream = TRUE)

  message("🎉 mestools package successfully deployed to GitHub!")
  message("📊 Repository: https://github.com/vanhungtran/mestools")

  # Final package build and install
  message("7. 📦 Final package build and install...")
  devtools::build()
  devtools::install()

  message("✅ Package installed locally!")

  return(invisible(TRUE))
}

#' Quick Update Function
#'
#' Fast deployment for quick updates without running checks
#'
#' @param commit_message Custom commit message
#' @return Invisible TRUE if successful
#' @export
deploy_mestools_quick <- function(commit_message = NULL) {
  deploy_mestools(
    commit_message = commit_message,
    run_checks = FALSE,
    skip_document = FALSE
  )
}

#' Force Reinitialization
#'
#' Completely reinitialize the Git repository and redeploy
#'
#' @param commit_message Custom commit message
#' @return Invisible TRUE if successful
#' @export
deploy_mestools_clean <- function(commit_message = NULL) {
  deploy_mestools(
    commit_message = commit_message,
    force_initial = TRUE,
    run_checks = TRUE,
    skip_document = FALSE
  )
}

# Additional utility functions for package management

#' Check Package Status
#'
#' Check the current status of the mestools package and Git repository
#'
#' @return List with status information
#' @export
check_mestools_status <- function() {
  pkg_dir <- file.path(dirname(here::here()), "mestools")
  setwd(pkg_dir)

  status <- list(
    package_dir = pkg_dir,
    dir_exists = dir.exists(pkg_dir),
    git_initialized = dir.exists(".git"),
    has_remote = nrow(gert::git_remote_list()) > 0,
    uncommitted_changes = nrow(gert::git_status()) > 0,
    last_commit = if (dir.exists(".git")) {
      log <- gert::git_log(max = 1)
      if (nrow(log) > 0) log else NULL
    } else NULL
  )

  # Print status summary
  message("📊 mestools Package Status:")
  message("  📁 Package directory: ", status$package_dir)
  message("  📁 Directory exists: ", status$dir_exists)
  message("  🔧 Git initialized: ", status$git_initialized)
  message("  🐙 Remote configured: ", status$has_remote)
  message("  📝 Uncommitted changes: ", status$uncommitted_changes)

  if (!is.null(status$last_commit)) {
    message("  ⏰ Last commit: ", status$last_commit$time)
    message("  💬 Last message: ", status$last_commit$message)
  }

  invisible(status)
}

#' Build and Install Only
#'
#' Just build and install the package without Git operations
#'
#' @return Invisible TRUE if successful
#' @export
build_install_mestools <- function() {
  pkg_dir <- file.path(dirname(here::here()), "mestools")
  setwd(pkg_dir)

  message("📦 Building mestools package...")
  devtools::document()
  devtools::build()
  devtools::install()
  message("✅ Package built and installed successfully!")

  invisible(TRUE)
}

# Main execution block (run when script is sourced)
if (sys.nframe() == 0) {
  # This block runs only when the script is executed directly,
  # not when it's sourced as part of a package

  message("🤖 mestools Deployment Script")
  message("==============================")
  message("Available functions:")
  message("• deploy_mestools() - Full deployment")
  message("• deploy_mestools_quick() - Quick update")
  message("• deploy_mestools_clean() - Clean reinitialization")
  message("• check_mestools_status() - Check current status")
  message("• build_install_mestools() - Build and install only")
  message("")
  message("Usage examples:")
  message('deploy_mestools("Custom commit message")')
  message('deploy_mestools_quick("Quick fix")')
  message("check_mestools_status()")
  message("")

  # Check status by default
  check_mestools_status()
}



