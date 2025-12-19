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

# Set up environment for GEO operations (when needed)
# options('download.file.method.GEOquery' = 'libcurl')

# Define GEO datasets (commented out to prevent execution during package load)
# gses <- unique(c(
#   "GSE32924", "GSE34248", "GSE36842", "GSE120721", "GSE16161", 
#   "GSE107361", "GSE5667", "GSE26952", "GSE75890", "GSE60709", "GSE133477", 
#   "GSE58558", "GSE130588", "GSE99802", "GSE133385", "GSE140684", "GSE121212", 
#   "GSE137430", "GSE65832", "GSE206391", "GSE147424", "GSE277961", "GSE280220", 
#   "GSE157194", "GSE289784", "GSE232127", "GSE186063", "GSE141570", "GSE176279", 
#   "GSE223799"
# ))

# Configuration (commented out to prevent execution during package load)
# data_dir <- paste0(here::here(), "/data")
# dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

#' Extract Metadata from GEO Object
#'
#' Helper function to extract metadata from GEO series matrix objects.
#'
#' @param meta GEO metadata object
#' @param sep Separator for characteristic parsing
#' @return Data frame with cleaned metadata
#' @export
extract_metadata <- function(meta, sep = ": ") {
  if (!requireNamespace("purrr", quietly = TRUE)) {
    stop("Please install purrr: install.packages('purrr')")
  }
  if (!requireNamespace("tidyr", quietly = TRUE)) {
    stop("Please install tidyr: install.packages('tidyr')")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Please install dplyr: install.packages('dplyr')")
  }
  if (!requireNamespace("janitor", quietly = TRUE)) {
    stop("Please install janitor: install.packages('janitor')")
  }
  
  # Extract characteristics without using pipes to avoid namespace issues
  characteristics_list <- purrr::map(meta@gsms, ~.x@header$characteristics_ch1)
  stacked_data <- utils::stack(characteristics_list)
  separated_data <- tidyr::separate(stacked_data, "values", into = c("feature", "value"), sep = sep)
  wide_data <- tidyr::pivot_wider(separated_data, names_from = "feature", values_from = "value")
  janitor::clean_names(wide_data)
}

#' Download and Process GEO Data
#'
#' Downloads and processes multiple GEO datasets with error handling.
#'
#' @param gse_list Character vector of GSE IDs to download
#' @param data_dir Directory to save downloaded files
#' @return List containing geo_data and expressions lists
#' @export
download_geo_data <- function(gse_list, data_dir = tempdir()) {
  if (!requireNamespace("GEOquery", quietly = TRUE)) {
    stop("Please install GEOquery: BiocManager::install('GEOquery')")
  }
  if (!requireNamespace("Biobase", quietly = TRUE)) {
    stop("Please install Biobase: BiocManager::install('Biobase')")
  }
  geo_data <- list()
  exps <- list()
  
  for (gse_id in gse_list) {
    tryCatch({
      message("Processing ", gse_id)
      
      # Download GEO data
      gse_data <- GEOquery::getGEO(gse_id, GSEMatrix = TRUE, getGPL = FALSE, destdir = data_dir)
      geo_data[[gse_id]] <- gse_data
      
      # Extract expression data
      if (length(gse_data) > 0) {
        exps[[gse_id]] <- Biobase::exprs(gse_data[[1]])
      }
      
    }, error = function(e) {
      warning("Failed to process ", gse_id, ": ", e$message)
    })
  }
  
  return(list(geo_data = geo_data, expressions = exps))
}

#' Process Metadata with Standardized Transformations
#'
#' Processes metadata for multiple GEO datasets with dataset-specific transformations.
#'
#' @param gse_list Character vector of GSE IDs to process
#' @param data_dir Directory containing downloaded files
#' @return List of processed metadata data frames
#' @export
process_metadata <- function(gse_list, data_dir = tempdir()) {
  if (!requireNamespace("GEOquery", quietly = TRUE)) {
    stop("Please install GEOquery: BiocManager::install('GEOquery')")
  }
  meta_list <- list()
  
  for (gse_id in gse_list) {
    tryCatch({
      message("Processing metadata for ", gse_id)
      
      # Get metadata
      meta <- GEOquery::getGEO(GEO = gse_id, GSEMatrix = FALSE, getGPL = TRUE, destdir = data_dir)
      
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

# Main execution (commented out to prevent execution during package load)
# if(TRUE) {
#   # Download and process data
#   geo_results <- download_geo_data(gses)
#   geo_data <- geo_results$geo_data
#   exps <- geo_results$expressions
#   
#   # Process metadata
#   META <- process_metadata(gses)
#   
#   # Summary
#   message("Processing complete:")
#   message("  - GEO datasets processed: ", length(geo_data))
#   message("  - Expression matrices: ", length(exps))
#   message("  - Metadata tables: ", length(META))
#   
#   # Remove failed downloads
#   geo_data <- geo_data[!sapply(geo_data, is.null)]
#   exps <- exps[!sapply(exps, is.null)]
#   META <- META[!sapply(META, is.null)]
# }

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
# Set options for GEOquery download method (when needed)
# options('download.file.method.GEOquery' = 'libcurl')

# Note: Required packages are loaded conditionally in functions:
# install.packages(c("tidyverse", "janitor", "data.table"))
# BiocManager::install(c("GEOquery", "here", "Biobase"))
# library(data.table)   # For fast reading of large files (fread) - load conditionally

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

# Define a list of RNA-seq datasets (commented out to prevent execution during package load)
# rnaseq_datasets <- c("GSE121212", "GSE137430", "GSE65832", "GSE157194",
#                      "GSE186063", "GSE141570", "GSE176279", "GSE223799")

# Define data directory and ensure it exists (commented out to prevent execution during package load)
# data_dir <- paste0(here::here(), "/data_geo")
# dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
# message("Data will be downloaded to: ", data_dir)


# Subset for a quick test (commented out to prevent execution during package load)
# # NOTE: This execution might take a few minutes as it downloads large files.
# test_datasets <- c("GSE121212", "GSE157194")
# message("\n--- Starting download for test datasets: ", paste(test_datasets, collapse = ", "), " ---")
# 
# # Run the improved function
# rnaseq_exps_test <- download_rnaseq_data_improved(test_datasets, data_dir)
# 
# # Check the results
# message("\n--- Results Summary ---")
# message("Number of datasets successfully processed: ", length(rnaseq_exps_test))
# if ("GSE121212" %in% names(rnaseq_exps_test)) {
#   message("GSE121212 expression matrix dimensions (Genes x Samples): ", paste(dim(rnaseq_exps_test$GSE121212), collapse = " x "))
# }
# if ("GSE157194" %in% names(rnaseq_exps_test)) {
#   message("GSE157194 expression matrix dimensions (Genes x Samples): ", paste(dim(rnaseq_exps_test$GSE157194), collapse = " x "))
# }
# # Clean up the test directory (optional)
# # unlink(data_dir, recursive = TRUE)


# --- GEO Dataset Functions ---

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
#' @param max_parallel Maximum number of datasets to process in parallel. Default is 3.
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
read_multiple_geo_datasets <- function(gse_ids, destdir = tempdir(), getGPL = TRUE, 
                                     AnnotGPL = TRUE, max_parallel = 3, sleep_between = 2) {
  
  message("🚀 Starting batch download of ", length(gse_ids), " GEO datasets...")
  
  # Initialize results list
  results <- list()
  failed_downloads <- character()
  
  # Process datasets
  for (i in seq_along(gse_ids)) {
    gse_id <- gse_ids[i]
    
    message("\n📦 Processing dataset ", i, "/", length(gse_ids), ": ", gse_id)
    
    # Download and process the dataset
    result <- read_geo_dataset(gse_id, destdir = destdir, getGPL = getGPL, AnnotGPL = AnnotGPL)
    
    if (!is.null(result)) {
      results[[gse_id]] <- result
    } else {
      failed_downloads <- c(failed_downloads, gse_id)
    }
    
    # Sleep between downloads to be respectful
    if (i < length(gse_ids)) {
      message("⏳ Waiting ", sleep_between, " seconds before next download...")
      Sys.sleep(sleep_between)
    }
  }
  
  # Summary
  successful <- length(results)
  failed <- length(failed_downloads)
  
  message("\n📊 Batch download completed!")
  message("   ✅ Successful: ", successful, "/", length(gse_ids))
  message("   ❌ Failed: ", failed, "/", length(gse_ids))
  
  if (failed > 0) {
    message("   📝 Failed datasets: ", paste(failed_downloads, collapse = ", "))
  }
  
  # Add summary information
  attr(results, "summary") <- list(
    total_requested = length(gse_ids),
    successful = successful,
    failed = failed,
    failed_ids = failed_downloads,
    download_time = Sys.time()
  )
  
  return(results)
}

#' Get GEO Dataset Summary
#'
#' Get a quick summary of one or more GEO datasets without downloading expression data.
#'
#' @param gse_ids A character vector of GSE IDs
#' @return A data.frame with summary information for each dataset
#' @export
#' @examples
#' \dontrun{
#' # Get summary information
#' summary_info <- get_geo_summary(c("GSE102628", "GSE102641"))
#' print(summary_info)
#' }
get_geo_summary <- function(gse_ids) {
  if (!requireNamespace("GEOquery", quietly = TRUE)) {
    stop("Please install GEOquery: BiocManager::install('GEOquery')")
  }
  
  message("📋 Fetching summary information for ", length(gse_ids), " datasets...")
  
  summary_df <- data.frame(
    GSE_ID = character(),
    Title = character(),
    Platform = character(),
    Samples = integer(),
    Features = integer(),
    Submission_Date = character(),
    Last_Update = character(),
    Status = character(),
    stringsAsFactors = FALSE
  )
  
  for (gse_id in gse_ids) {
    tryCatch({
      message("   Fetching ", gse_id, "...")
      gse <- GEOquery::getGEO(gse_id, getGPL = FALSE)
      eset <- gse[[1]]
      
      summary_df <- rbind(summary_df, data.frame(
        GSE_ID = gse_id,
        Title = GEOquery::experimentData(eset)@title,
        Platform = GEOquery::annotation(eset),
        Samples = ncol(GEOquery::exprs(eset)),
        Features = nrow(GEOquery::exprs(eset)),
        Submission_Date = GEOquery::experimentData(eset)@pubMedIds,
        Last_Update = as.character(GEOquery::experimentData(eset)@other$last_update_date),
        Status = "Success",
        stringsAsFactors = FALSE
      ))
      
    }, error = function(e) {
      summary_df <<- rbind(summary_df, data.frame(
        GSE_ID = gse_id,
        Title = "Error",
        Platform = "Error",
        Samples = NA,
        Features = NA,
        Submission_Date = "Error",
        Last_Update = "Error",
        Status = paste("Failed:", e$message),
        stringsAsFactors = FALSE
      ))
    })
  }
  
  message("✅ Summary completed!")
  return(summary_df)
}

#' Default GSE Dataset List
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

#' Process All Default GEO Datasets
#'
#' Convenience function to download and process all datasets in the default GSE list.
#'
#' @param destdir Directory to save downloaded files. Default creates a 'geo_data' folder.
#' @param getGPL Logical. Whether to download platform annotation data.
#' @param AnnotGPL Logical. Whether to annotate the expression data with gene symbols.
#' @param subset_size If specified, only process the first N datasets (for testing).
#' @return A named list containing all processed datasets
#' @export
#' @examples
#' \dontrun{
#' # Process first 5 datasets for testing
#' test_results <- process_all_geo_datasets(subset_size = 5)
#' 
#' # Process all datasets (this will take a while!)
#' all_results <- process_all_geo_datasets()
#' }
process_all_geo_datasets <- function(destdir = "geo_data", getGPL = TRUE, 
                                   AnnotGPL = TRUE, subset_size = NULL) {
  
  # Create destination directory if it doesn't exist
  if (!dir.exists(destdir)) {
    dir.create(destdir, recursive = TRUE)
    message("📁 Created directory: ", destdir)
  }
  
  # Get the GSE list
  gse_list <- get_default_gse_list()
  
  # Subset if requested
  if (!is.null(subset_size)) {
    gse_list <- gse_list[1:min(subset_size, length(gse_list))]
    message("🔍 Processing subset of ", length(gse_list), " datasets for testing")
  }
  
  message("🎯 Target directory: ", normalizePath(destdir))
  message("📊 Total datasets to process: ", length(gse_list))
  
  # Process the datasets
  results <- read_multiple_geo_datasets(
    gse_ids = gse_list,
    destdir = destdir,
    getGPL = getGPL,
    AnnotGPL = AnnotGPL
  )
  
  return(results)
}