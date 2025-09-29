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

  # Step 5: Check and setup upstream branch
  message("5. 🔄 Checking upstream branch configuration...")
  current_branch <- gert::git_branch()

  # Check if upstream is configured
  upstream_configured <- tryCatch({
    gert::git_branch_remote()
    TRUE
  }, error = function(e) FALSE)

  if (!upstream_configured && has_remote) {
    message("   ⚙️  Setting upstream branch for '", current_branch, "'...")
    tryCatch({
      # Try to set upstream
      gert::git_push(remote = "origin", set_upstream = TRUE)
      message("   ✅ Upstream branch configured")
    }, error = function(e) {
      message("   ⚠️  Could not set upstream automatically: ", e$message)
    })
  } else if (upstream_configured) {
    message("   ✅ Upstream branch already configured")
  }

  # Step 6: Sync with remote (only if upstream is configured)
  if (has_remote && git_initialized && upstream_configured) {
    message("6. 🔄 Syncing with remote repository...")
    tryCatch({
      # Fetch latest changes from remote
      gert::git_fetch()

      # Check if we're behind the remote
      ahead_behind <- gert::git_ahead_behind()
      if (ahead_behind$behind > 0) {
        message("   📥 Pulling ", ahead_behind$behind, " commit(s) from remote...")
        gert::git_pull()
        message("   ✅ Successfully pulled remote changes")
      } else {
        message("   ✅ Local repository is up to date with remote")
      }
    }, error = function(e) {
      message("   ⚠️  Sync failed: ", e$message)
    })
  } else {
    message("6. ⚠️  Skipping sync (no upstream configured)")
  }

  # Step 7: Commit changes
  message("7. 💾 Committing changes...")

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

  # Step 8: Push to GitHub
  message("8. 🚀 Pushing to GitHub...")
  tryCatch({
    # Check if we have commits to push
    ahead_behind <- tryCatch({
      gert::git_ahead_behind()
    }, error = function(e) {
      list(ahead = 1, behind = 0) # Assume we're ahead if we can't check
    })

    if (ahead_behind$ahead > 0 || !upstream_configured) {
      message("   Pushing to remote...")
      gert::git_push(remote = "origin", set_upstream = !upstream_configured)
      message("   ✅ Successfully pushed to GitHub!")
    } else {
      message("   📭 No new commits to push")
    }
  }, error = function(e) {
    if (grepl("contains commits that are not present locally", e$message)) {
      message("   ❌ Push failed: Remote has commits not present locally")
      message("   💡 Run: gert::git_pull() to merge remote changes first")
    } else if (grepl("No upstream", e$message)) {
      message("   ❌ Push failed: No upstream branch configured")
      message("   💡 Run: gert::git_push(remote = 'origin', set_upstream = TRUE)")
    } else {
      message("   ❌ Push failed: ", e$message)
    }
    return(invisible(FALSE))
  })

  message("🎉 mestools package deployment completed!")
  message("📊 Repository: https://github.com/vanhungtran/mestools")

  # Final package build and install
  message("9. 📦 Final package build and install...")
  devtools::build()
  devtools::install()

  message("✅ Package installed locally!")

  return(invisible(TRUE))
}

#' Safe Deployment with Upstream Handling
#'
#' Pull remote changes and set upstream before pushing
#'
#' @param commit_message Custom commit message
#' @return Invisible TRUE if successful
#' @export
deploy_mestools_safe <- function(commit_message = NULL) {
  pkg_dir <- file.path(dirname(here::here()), "mestools")
  setwd(pkg_dir)

  message("🔄 Safe deployment - ensuring upstream configuration...")

  # Step 1: Ensure upstream is configured
  current_branch <- gert::git_branch()

  # Check if upstream exists
  upstream_exists <- tryCatch({
    gert::git_branch_remote()
    TRUE
  }, error = function(e) FALSE)

  if (!upstream_exists) {
    message("1. ⚙️  Setting upstream branch...")
    tryCatch({
      gert::git_push(remote = "origin", set_upstream = TRUE)
      message("   ✅ Upstream configured")
    }, error = function(e) {
      message("   ❌ Failed to set upstream: ", e$message)
      message("   💡 Trying manual upstream setup...")

      # Manual upstream setup
      system(paste("git branch --set-upstream-to=origin/", current_branch, " ", current_branch, sep = ""))
    })
  }

  # Step 2: Pull changes
  message("2. 📥 Pulling remote changes...")
  tryCatch({
    gert::git_pull()
    message("   ✅ Successfully pulled remote changes")
  }, error = function(e) {
    message("   ⚠️  Pull failed: ", e$message)
  })

  # Step 3: Continue with normal deployment
  deploy_mestools(
    commit_message = commit_message,
    run_checks = FALSE,
    skip_document = FALSE
  )
}

#' Manual Upstream Setup
#'
#' Manually set upstream branch for the current branch
#'
#' @return Invisible TRUE if successful
#' @export
setup_upstream <- function() {
  pkg_dir <- file.path(dirname(here::here()), "mestools")
  setwd(pkg_dir)

  current_branch <- gert::git_branch()
  message("Setting upstream for branch: ", current_branch)

  tryCatch({
    # Method 1: Using gert
    gert::git_push(remote = "origin", set_upstream = TRUE)
    message("✅ Upstream configured successfully")
  }, error = function(e) {
    message("❌ Method 1 failed: ", e$message)

    # Method 2: Using system command
    message("Trying system command...")
    cmd <- paste("git push --set-upstream origin", current_branch)
    system(cmd)
    message("✅ Upstream configuration attempted with system command")
  })

  return(invisible(TRUE))
}

# Other functions remain similar but with improved upstream handling...

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

#' Check Package Status with Upstream Info
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
    has_remote = nrow(gert::git_remote_list()) > 0
  )

  if (status$git_initialized) {
    status$branch <- gert::git_branch()

    # Check upstream
    status$upstream_configured <- tryCatch({
      gert::git_branch_remote()
      TRUE
    }, error = function(e) FALSE)

    status$uncommitted_changes = nrow(gert::git_status()) > 0

    if (status$upstream_configured) {
      status$ahead_behind <- tryCatch({
        gert::git_ahead_behind()
      }, error = function(e) {
        list(ahead = 0, behind = 0)
      })
    }

    status$last_commit <- tryCatch({
      log <- gert::git_log(max = 1)
      if (nrow(log) > 0) log else NULL
    }, error = function(e) NULL)
  }

  # Print status summary
  message("📊 mestools Package Status:")
  message("  📁 Package directory: ", status$package_dir)
  message("  📁 Directory exists: ", status$dir_exists)
  message("  🔧 Git initialized: ", status$git_initialized)
  message("  🐙 Remote configured: ", status$has_remote)

  if (status$git_initialized) {
    message("  🌿 Current branch: ", status$branch)
    message("  🔗 Upstream configured: ", status$upstream_configured)
    message("  📝 Uncommitted changes: ", status$uncommitted_changes)

    if (status$upstream_configured && !is.null(status$ahead_behind)) {
      message("  📊 Ahead/Behind: ", status$ahead_behind$ahead, " ahead, ",
              status$ahead_behind$behind, " behind")
    }

    if (!is.null(status$last_commit)) {
      message("  ⏰ Last commit: ", status$last_commit$time)
      message("  💬 Last message: ", status$last_commit$message)
    }
  }

  invisible(status)
}

# Main execution block
if (sys.nframe() == 0) {
  message("🤖 mestools Deployment Script")
  message("==============================")
  message("Available functions:")
  message("• deploy_mestools() - Full deployment")
  message("• deploy_mestools_safe() - Safe deployment with upstream setup")
  message("• setup_upstream() - Manual upstream configuration")
  message("• deploy_mestools_quick() - Quick update")
  message("• check_mestools_status() - Check current status")
  message("")
  message("To fix your current error, run:")
  message("setup_upstream()")
  message("Then:")
  message('deploy_mestools_safe("Your commit message")')
  message("")

  # Check status by default
  check_mestools_status()
}
