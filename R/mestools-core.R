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




# Improved GEO Data Processing Code

# Set up environment
options('download.file.method.GEOquery' = 'libcurl')
library(GEOquery)
library(tidyverse)
library(janitor)
library(here)

# Define GEO datasets
gses <- unique(c(
  "GSE32924", "GSE34248", "GSE36842", "GSE120721", "GSE16161", 
  "GSE107361", "GSE5667", "GSE26952", "GSE75890", "GSE60709", "GSE133477", 
  "GSE58558", "GSE130588", "GSE99802", "GSE133385", "GSE140684", "GSE121212", 
  "GSE137430", "GSE65832", "GSE206391", "GSE147424", "GSE277961", "GSE280220", 
  "GSE157194", "GSE289784", "GSE232127", "GSE186063", "GSE141570", "GSE176279", 
  "GSE223799"
))

# Configuration
data_dir <- paste0(here::here(), "/data")
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

# Helper function to extract metadata
extract_metadata <- function(meta, sep = ": ") {
  purrr::map(meta@gsms, ~.x@header$characteristics_ch1) %>%
    stack() %>%
    tidyr::separate(values, into = c("feature", "value"), sep = sep) %>%
    pivot_wider(names_from = feature, values_from = value) %>%
    janitor::clean_names()
}

# Download and process GEO data with error handling
download_geo_data <- function(gse_list) {
  geo_data <- list()
  exps <- list()
  
  for (gse_id in gse_list) {
    tryCatch({
      message("Processing ", gse_id)
      
      # Download GEO data
      gse_data <- getGEO(gse_id, GSEMatrix = TRUE, getGPL = FALSE, destdir = data_dir)
      geo_data[[gse_id]] <- gse_data
      
      # Extract expression data
      if (length(gse_data) > 0) {
        exps[[gse_id]] <- exprs(gse_data[[1]])
      }
      
    }, error = function(e) {
      warning("Failed to process ", gse_id, ": ", e$message)
    })
  }
  
  return(list(geo_data = geo_data, expressions = exps))
}

# Process metadata with standardized transformations
process_metadata <- function(gse_list) {
  meta_list <- list()
  
  for (gse_id in gse_list) {
    tryCatch({
      message("Processing metadata for ", gse_id)
      
      # Get metadata
      meta <- getGEO(GEO = gse_id, GSEMatrix = FALSE, getGPL = TRUE, destdir = data_dir)
      
      # Extract base metadata
      if (gse_id == "GSE5667") {
        current_meta <- extract_metadata(meta, sep = ":")
      } else {
        current_meta <- extract_metadata(meta)
      }
      
      # Add project identifier
      current_meta$project <- gse_id
      
      # Apply dataset-specific transformations
      current_meta <- apply_dataset_transformations(gse_id, current_meta)
      
      meta_list[[gse_id]] <- current_meta
      
    }, error = function(e) {
      warning("Failed to process metadata for ", gse_id, ": ", e$message)
    })
  }
  
  return(meta_list)
}

# Dataset-specific transformations
apply_dataset_transformations <- function(gse_id, meta_df) {
  transformations <- list(
    GSE34248 = function(df) {
      df %>% mutate(tissue_type = ifelse(tissue_type == "Lesional skin", "LS", "NL"))
    },
    
    GSE32924 = function(df) {
      df %>% mutate(
        status = ifelse(condition == "Normal", "HC", "AD"),
        tissue_type = paste(condition, tissue)
      )
    },
    
    GSE36842 = function(df) {
      df %>% mutate(
        status = ifelse(group == "Normal", "HC", "AD"),
        tissue_type = group
      )
    },
    
    GSE120721 = function(df) {
      df %>% mutate(
        status = ifelse(diagnosis == "Normal", "HC", "AD"),
        tissue_type = paste0(tissue, region),
        id = patient_id
      )
    },
    
    GSE16161 = function(df) {
      df %>% mutate(
        status = case_when(
          condition == "Normal" ~ "HC",
          condition == "Psoriasis" ~ "Psoriasis",
          TRUE ~ "AD"
        ),
        tissue_type = tissue,
        id = indivisual
      )
    },
    
    GSE107361 = function(df) {
      df %>% mutate(
        lesion = ifelse(lesion == "N", "Normal", lesion),
        status = ifelse(lesion == "Normal", "HC", "AD"),
        tissue_type = paste(lesion, tissue),
        id = indivisual
      )
    },
    
    GSE140684 = function(df) {
      df %>% mutate(
        status = ifelse(diagnosis == "atopic dermatitis", "AD", ""),
        tissue_type = paste(tissue_type, tissue),
        id = subject_id
      )
    },
    
    GSE75890 = function(df) {
      df %>% mutate(
        id = ind,
        status = ifelse(diagnosis == "Atopic Dermatitis", "AD", "HC"),
        tissue_type = tissue
      )
    },
    
    GSE26952 = function(df) {
      df %>% mutate(
        id = ind,
        tissue_type = "NL-epidermal",
        status = recode(disease_group, 
                       'Nonatopic Control' = 'HC', 
                       "Psoriasis" = "Psoriasis", 
                       "Atopic Dermatitis" = "AD")
      )
    },
    
    GSE60709 = function(df) {
      df %>% mutate(
        id = ind,
        tissue_type = ifelse(tissue == "lesional skin", "LS", "NL"),
        sex = gender,
        status = recode(phenotype, 'Control' = 'HC', "Atopic Dermatitis" = "AD")
      )
    },
    
    GSE133477 = function(df) {
      df %>% mutate(
        id = patient_id,
        timepoint = week,
        batch = batch_date,
        status = "AD",
        tissue_type = ifelse(tissue == "LESIONAL", "LS", "NL")
      )
    },
    
    GSE58558 = function(df) {
      df %>% mutate(
        status = "AD",
        tissue_type = ifelse(lesions == "lesional", "LS", "NL"),
        id = ind,
        timepoint = time,
        sex = gender,
        treatment = "cyclosporine A"
      )
    },
    
    GSE130588 = function(df) {
      df %>% mutate(
        status = "AD",
        tissue_type = tissue,
        id = subject_id,
        timepoint = time
      )
    },
    
    GSE99802 = function(df) {
      df %>% mutate(
        status = "AD",
        tissue_type = tissue,
        id = patient_id,
        timepoint = week,
        batch = batch_date,
        treatment = ifelse(treatment == "Drug", "Fezakinumab", "Placebo")
      )
    },
    
    GSE133385 = function(df) {
      df %>% mutate(
        status = "AD",
        tissue_type = tissue,
        id = patient_id,
        timepoint = week,
        batch = batch_date
      )
    },
    
    GSE137430 = function(df) {
      df %>% mutate(
        status = "AD",
        tissue_type = skin_type,
        id = ptid
      )
    },
    
    GSE121212 = function(df) {
      df %>% mutate(
        status = patients_condition,
        tissue_type = skin_type,
        id = ind,
        tissue_type = recode(tissue_type, 
                           'lesional' = 'LS',
                           "non-lesional" = "NLS")
      )
    },
    
    GSE65832 = function(df) {
      df %>% mutate(
        status = "AD",
        sex = gender,
        id = ind,
        tissue_type = recode(tissue, 
                           'lesional skin' = 'LS',
                           "non-lesional skin" = "NLS")
      )
    },
    
    GSE206391 = function(df) {
      df %>% mutate(
        id = patient,
        tissue_type = recode(tissue_type, 
                           'lesional' = 'LS',
                           "non lesional" = "NLS"),
        status = recode(disease, 'atopic dermatitis' = 'AD')
      )
    },
    
    GSE147424 = function(df) {
      df %>% mutate(
        id = ind,
        tissue_type = recode(tissue_type, 
                           'Lesional' = 'LS',
                           "Non-lesional" = "NLS"),
        status = recode(disease, 
                       'Atopic dermatitis' = 'AD',
                       "Healthy" = "HC")
      )
    }
  )
  
  # Apply transformation if exists, otherwise return original
  if (gse_id %in% names(transformations)) {
    return(transformations[[gse_id]](meta_df))
  } else {
    return(meta_df)
  }
}

# Main execution
if(TRUE) {
  # Download and process data
  geo_results <- download_geo_data(gses)
  geo_data <- geo_results$geo_data
  exps <- geo_results$expressions
  
  # Process metadata
  META <- process_metadata(gses)
  
  # Summary
  message("Processing complete:")
  message("  - GEO datasets processed: ", length(geo_data))
  message("  - Expression matrices: ", length(exps))
  message("  - Metadata tables: ", length(META))
  
  # Remove failed downloads
  geo_data <- geo_data[!sapply(geo_data, is.null)]
  exps <- exps[!sapply(exps, is.null)]
  META <- META[!sapply(META, is.null)]
}

# Optional: Save results
save_results <- function() {
  saveRDS(geo_data, file.path(data_dir, "geo_data.rds"))
  saveRDS(exps, file.path(data_dir, "expression_data.rds"))
  saveRDS(META, file.path(data_dir, "metadata.rds"))
  message("Results saved to ", data_dir)
}

# Uncomment to save results
# save_results()




# IMPROVED GEO RNA-SEQ DOWNLOAD FUNCTION AND EXAMPLE

# ==============================================================================
# 1. SETUP AND LIBRARIES
# ==============================================================================
# Set options for GEOquery download method
options('download.file.method.GEOquery' = 'libcurl')

# Load necessary libraries
# Note: Ensure these packages are installed: install.packages(c("tidyverse", "janitor", "data.table"))
# Bioconductor packages: BiocManager::install(c("GEOquery", "here", "Biobase"))
library(GEOquery)
library(Biobase)      # For exprs()
library(tidyverse)    # For general data manipulation
library(janitor)      # For data cleaning (though not explicitly used in this snippet)
library(here)         # For creating platform-independent paths
library(data.table)   # For fast reading of large files (fread)

# ==============================================================================
# 2. IMPROVED FUNCTION DEFINITION
# ==============================================================================
download_rnaseq_data_improved <- function(gse_list, data_dir) {
  # Base URL for GEO file downloads (common structure for NCBI-provided counts)
  urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
  exps_list <- list()

  # Common RNA-seq counts file name patterns to try for robustness
  rnaseq_patterns <- c(
    "_raw_counts_GRCh38.p13_NCBI.tsv.gz", # Original, most specific pattern
    "_raw_counts_GRCh38.tsv.gz",
    "_raw_counts_GRCh37.tsv.gz",
    "_counts.tsv.gz",
    "_counts.txt.gz"
  )

  for (gse_id in gse_list) {
    message("Attempting to download RNA-seq data for ", gse_id)
    download_successful <- FALSE

    # 1. Try the direct GEO RNA-seq counts link with multiple patterns
    for (pattern in rnaseq_patterns) {
      if (download_successful) break

      tryCatch({
        filename <- paste0(gse_id, pattern)
        # Construct the download URL using the consistent parts and the variable filename
        download_url <- paste(urld, paste0("acc=", gse_id),
                              paste0("file=", filename),
                              sep = "&")

        # Download and read counts data using fread (handles gzipped URLs)
        counts_data <- as.matrix(
          data.table::fread(download_url, header = TRUE, colClasses = "integer"),
          rownames = 1
        )

        exps_list[[gse_id]] <- counts_data
        message("Successfully downloaded RNA-seq data for ", gse_id, " using pattern: ", pattern)
        download_successful <- TRUE

      }, error = function(e) {
        # Fail silently for pattern attempts to allow trying the next pattern
      })
    }

    # 2. Fallback: If direct download fails, try the standard GEOquery approach
    if (!download_successful) {
      message("Direct RNA-seq download failed for ", gse_id, ". Falling back to getGEO.")
      tryCatch({
        # Download the standard series matrix file
        gse_data <- GEOquery::getGEO(gse_id, GSEMatrix = TRUE, getGPL = FALSE, destdir = data_dir)
        if (length(gse_data) > 0) {
          # Extract the expression data (which might contain counts)
          exps_list[[gse_id]] <- Biobase::exprs(gse_data[[1]])
          message("Successfully downloaded and extracted data for ", gse_id, " using GEOquery fallback.")
          download_successful <- TRUE
        }
      }, error = function(e) {
        warning("Failed to process ", gse_id, " using GEOquery fallback: ", e$message)
      })
    }

    if (!download_successful) {
      warning("Failed to download RNA-seq data for ", gse_id, " using all methods.")
    }
  }

  return(exps_list)
}

# ==============================================================================
# 3. EXAMPLE USAGE
# ==============================================================================

# Define a list of RNA-seq datasets
rnaseq_datasets <- c("GSE121212", "GSE137430", "GSE65832", "GSE157194",
                     "GSE186063", "GSE141570", "GSE176279", "GSE223799")

# Define data directory and ensure it exists
data_dir <- paste0(here::here(), "/data_geo")
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
message("Data will be downloaded to: ", data_dir)


# Subset for a quick test (GSE121212 is known to work with the direct NCBI link)
# NOTE: This execution might take a few minutes as it downloads large files.
test_datasets <- c("GSE121212", "GSE157194")
message("\n--- Starting download for test datasets: ", paste(test_datasets, collapse = ", "), " ---")

# Run the improved function
rnaseq_exps_test <- download_rnaseq_data_improved(test_datasets, data_dir)

# Check the results
message("\n--- Results Summary ---")
message("Number of datasets successfully processed: ", length(rnaseq_exps_test))
if ("GSE121212" %in% names(rnaseq_exps_test)) {
  message("GSE121212 expression matrix dimensions (Genes x Samples): ", paste(dim(rnaseq_exps_test$GSE121212), collapse = " x "))
}
if ("GSE157194" %in% names(rnaseq_exps_test)) {
  message("GSE157194 expression matrix dimensions (Genes x Samples): ", paste(dim(rnaseq_exps_test$GSE157194), collapse = " x "))
}
# Clean up the test directory (optional)
# unlink(data_dir, recursive = TRUE)