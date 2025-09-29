
# Load required packages
library(usethis)
library(gert)
library(devtools)
library(here)

# Set working directory to your mestools package
pkg_dir <- file.path(dirname(here::here()), "mestools")
setwd(pkg_dir)

# Remove existing .git folder if it exists (to start fresh)
if (dir.exists(".git")) {
  unlink(".git", recursive = TRUE, force = TRUE)
  message("🧹 Removed existing Git repository")
}

# Step 1: Initialize new Git repository
message("1. 🔧 Initializing new Git repository...")
usethis::use_git()

# Step 2: Configure Git user (replace with your info)
message("2. 👤 Configuring Git user...")
gert::git_config_set("user.name", "Van Hung Tran")
gert::git_config_set("user.email", "tranhungydhcm@gmail.com")

# Step 3: Add all files to Git
message("3. 📁 Adding files to Git...")
gert::git_add(".")

# Step 4: Make initial commit
message("4. 💾 Making initial commit...")
gert::git_commit("Initial commit: mestools package with utility functions for data analysis and package development")

# Step 5: Create GitHub repository remotely
message("5. 🐙 Creating GitHub repository...")
# Note: You'll need a GitHub Personal Access Token (PAT)
# If you haven't set it up, run: usethis::create_github_token()
usethis::use_github()

# Alternative if use_github() doesn't work - manual remote setup:
if (!length(gert::git_remote_list()$url) > 0) {
  message("Setting up GitHub remote manually...")
  gert::git_remote_add(
    name = "origin",
    url = "https://github.com/vanhungtran/mestools.git"
  )
}

# Step 6: Push to GitHub
message("6. 🚀 Pushing to GitHub...")
gert::git_push(remote = "origin", set_upstream = TRUE)

message("🎉 Successfully initialized new Git repository and pushed to GitHub!")




here::here("mestools")

devtools::document()

devtools::check()



# Initialize Git repository
usethis::use_git()

# Create GitHub repository and set remote
usethis::use_github()

# Alternatively, if we already have a remote, we can set it manually
# usethis::use_git_remote("origin", "https://github.com/vanhungtran/mestools.git")

# Commit all changes
gert::git_add(".")
gert::git_commit("Initial commit: mestools package with utility functions")

# Push to GitHub
gert::git_push()
