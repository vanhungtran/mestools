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




library(httr)

#' Check if a specific NCBI-processed RNA-seq file exists for a GEO Series
#'
#' @param gse_id A character string with the GEO Series ID (e.g., "GSE121212").
#' @return A logical value: `TRUE` if the file exists and is accessible, `FALSE` otherwise.

check_geo_file_exists <- function(gse_id) {

  # Construct the exact URL based on the pattern in your script
  file_name <- paste0(gse_id, "_raw_counts_GRCh38.p13_NCBI.tsv.gz")
  url <- paste0("https://www.ncbi.nlm.nih.gov/geo/download/",
                "?format=file&type=rnaseq_counts",
                "&acc=", gse_id,
                "&file=", file_name)

  message(paste("Checking:", gse_id))

  # Use tryCatch to handle potential network errors gracefully
  result <- tryCatch({

    # Send a HEAD request to get headers without downloading the body
    response <- httr::HEAD(url)

    # A status code of 200 means the request was successful and the file exists
    return(httr::status_code(response) == 200)

  }, error = function(e) {

    # If any error occurs (e.g., no internet connection), return FALSE
    warning(paste("An error occurred for", gse_id, ":", e$message))
    return(FALSE)
  })

  return(result)
}




#' Read GEO Dataset
#'
#' Downloads and processes a single GEO dataset using GEOquery.
#'
#' @param gse_id A character string specifying the GSE ID (e.g., "GSE102628")
#' @param destdir Directory to save downloaded files. Default is tempdir().
#' @param getGPL Logical. Whether to download platform annotation data.
#' @param AnnotGPL Logical. Whether to annotate the expression data with gene symbols.
#' @return A list containing the GEOquery object and processed expression matrix
#' @export
#' @examples
#' \dontrun{
#' # Read a single GEO dataset
#' result <- read_geo_dataset("GSE102628")
#' expression_matrix <- result$expression_matrix
#' metadata <- result$phenotype_data
#' }
read_geo_dataset <- function(gse_id, destdir = tempdir(), getGPL = TRUE, AnnotGPL = TRUE) {
  if (!requireNamespace("GEOquery", quietly = TRUE)) {
    stop("Please install GEOquery: BiocManager::install('GEOquery')")
  }

  message("📥 Downloading GEO dataset: ", gse_id)

  tryCatch({
    # Download the GEO dataset
    gse <- GEOquery::getGEO(gse_id, destdir = destdir, getGPL = getGPL, AnnotGPL = AnnotGPL)

    # Extract the first (and usually only) dataset
    if (length(gse) > 1) {
      message("⚠️  Multiple platforms found. Using the first one.")
    }
    eset <- gse[[1]]

    # Extract expression matrix
    expression_matrix <- GEOquery::exprs(eset)

    # Extract phenotype data (sample metadata)
    phenotype_data <- GEOquery::pData(eset)

    # Extract feature data (probe/gene annotations)
    feature_data <- GEOquery::fData(eset)

    # Basic information
    n_samples <- ncol(expression_matrix)
    n_features <- nrow(expression_matrix)

    message("✅ Successfully processed ", gse_id)
    message("   📊 Features: ", n_features)
    message("   🧪 Samples: ", n_samples)
    message("   📋 Platform: ", GEOquery::annotation(eset))

    return(list(
      gse_object = eset,
      expression_matrix = expression_matrix,
      phenotype_data = phenotype_data,
      feature_data = feature_data,
      gse_id = gse_id,
      n_samples = n_samples,
      n_features = n_features,
      platform = GEOquery::annotation(eset)
    ))

  }, error = function(e) {
    message("❌ Error processing ", gse_id, ": ", e$message)
    return(NULL)
  })
}

#' Read Multiple GEO Datasets
#'
#' Downloads and processes multiple GEO datasets in batch.
#'
#' @param gse_ids A character vector of GSE IDs
#' @param destdir Directory to save downloaded files. Default is tempdir().
#' @param getGPL Logical. Whether to download platform annotation data.
#' @param AnnotGPL Logical. Whether to annotate the expression data with gene symbols.
#' @param sleep_between Sleep time in seconds between downloads to be respectful to NCBI servers.
#' @return A named list where each element contains the processed dataset
#' @export
#' @examples
#' \dontrun{
#' # Read multiple GEO datasets
#' gse_list <- c("GSE102628", "GSE102641", "GSE102725")
#' results <- read_multiple_geo_datasets(gse_list)
#'
#' # Access individual datasets
#' dataset1 <- results$GSE102628$expression_matrix
#' }
read_multiple_geo_datasets <- function(gse_list) {
  results <- list()

  for (gse_id in gse_list) {
    print(paste("📦 Processing dataset:", gse_id))

    # Use a tryCatch block to handle errors gracefully
    gse_data <- tryCatch({

      # Download the GEO dataset
      gse <- getGEO(gse_id, GSEMatrix = TRUE, AnnotGPL = TRUE)

      # IMPORTANT: Use assay() instead of exprs()
      # If multiple platforms exist, it returns a list, so we take the first element gse[[1]]
      if (length(gse) > 1) {
        print("⚠️ Multiple platforms found. Using the first one.")
      }
      expression_matrix <- assay(gse[[1]])

      # Return the expression matrix (or the full object if you prefer)
      expression_matrix

    }, error = function(e) {
      print(paste("❌ Error processing", gse_id, ":", e$message))
      return(NULL) # Return NULL if an error occurs
    })

    results[[gse_id]] <- gse_data
  }

  return(results)
}

# Now you can run your code again
gse_list <- c("GSE102628", "GSE102641", "GSE102725")
results <- read_multiple_geo_datasets(gse_list)

#' Get Default GSE Dataset List
#'
#' Returns the predefined list of GSE IDs for batch processing.
#'
#' @return A character vector of GSE IDs
#' @export
#' @examples
#' # Get the default GSE list
#' gse_list <- get_default_gse_list()
#' length(gse_list)
get_default_gse_list <- function() {
  c("GSE102628", "GSE102641", "GSE102725", "GSE103489", "GSE104509", "GSE106087",
    "GSE106992", "GSE107361", "GSE107871", "GSE109182", "GSE109248", "GSE111053",
    "GSE111054", "GSE111055", "GSE11307", "GSE114286", "GSE114729", "GSE116486",
    "GSE117239", "GSE117405", "GSE11903", "GSE120721", "GSE120899", "GSE121212",
    "GSE123785", "GSE123786", "GSE123787", "GSE124700", "GSE124701", "GSE130588",
    "GSE133385", "GSE133477", "GSE13355", "GSE14905", "GSE16161", "GSE18686",
    "GSE18948", "GSE20264", "GSE24767", "GSE26866", "GSE26952", "GSE2737",
    "GSE27887", "GSE30355", "GSE30768", "GSE30999", "GSE31652", "GSE32407",
    "GSE32473", "GSE32620", "GSE32924", "GSE34248", "GSE36381", "GSE36387",
    "GSE36842", "GSE38039", "GSE40033", "GSE40263", "GSE41662", "GSE41663",
    "GSE41664", "GSE41745", "GSE41905", "GSE42305", "GSE42632", "GSE47598",
    "GSE47751", "GSE47944", "GSE47965", "GSE48586", "GSE50598", "GSE50614",
    "GSE50790", "GSE51440", "GSE52361", "GSE52471", "GSE53431", "GSE53552",
    "GSE54456", "GSE55201", "GSE5667", "GSE57225", "GSE57376", "GSE57383",
    "GSE57386", "GSE57405", "GSE58121", "GSE58558", "GSE58749", "GSE59294",
    "GSE60481", "GSE60709", "GSE60971", "GSE61281", "GSE62408", "GSE63079",
    "GSE63741", "GSE63979", "GSE63980", "GSE121212", "GSE141570", "GSE65832",
    "GSE6601", "GSE66511", "GSE6710", "GSE67785", "GSE67853", "GSE68923",
    "GSE68924", "GSE68939", "GSE69967", "GSE72246", "GSE74697", "GSE75343",
    "GSE75890", "GSE77719", "GSE78023", "GSE78097", "GSE79704", "GSE80047",
    "GSE80429", "GSE82140", "GSE83582", "GSE83645", "GSE85034", "GSE86451",
    "GSE89725", "GSE92472", "GSE93423", "GSE99802", "GSE224783", "GSE140684",
    "GSE15719", "GSE95759", "GSE67785", "GSE157194", "GSE182740", "GSE261704",
    "GSE283265")
}
