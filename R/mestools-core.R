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

# Example usage (commented out - run manually if needed):
# gse_list <- c("GSE102628", "GSE102641", "GSE102725")
# results <- read_multiple_geo_datasets(gse_list)

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

#' Plot Combined Heatmap and Log2 Fold Change
#'
#' Creates a combined visualization with a heatmap of gene expression on the left
#' and a bar plot of log2 fold change on the right. Gene names and plot heights are
#' perfectly synchronized between both plots. Useful for visualizing differential
#' expression results.
#'
#' @param expression_matrix Numeric matrix of gene expression values (genes as rows, samples as columns)
#' @param log2fc Numeric vector of log2 fold change values (must match row order of expression_matrix)
#' @param sample_groups Optional character vector indicating sample groups for column annotation
#' @param cluster_rows Logical. Whether to cluster rows (genes). Default: TRUE
#' @param cluster_cols Logical. Whether to cluster columns (samples). Default: TRUE
#' @param scale Character. Scale the data: "row", "column", or "none". Default: "row"
#' @param show_rownames Logical. Show gene names on both heatmap and FC plot. Default: TRUE (auto-hidden if > 50 genes)
#' @param show_colnames Logical. Show sample names. Default: TRUE
#' @param heatmap_colors Vector of colors for heatmap gradient. Default: blue-white-red
#' @param fc_colors Vector of 2 colors for fold change bars (down, up). Default: blue and red
#' @param fc_threshold Numeric. Highlight genes with |log2FC| above this threshold. Default: 1
#' @param title Character. Main plot title. Default: NULL
#' @param gene_labels Optional character vector of gene names (uses rownames if NULL)
#' @param width_ratio Numeric vector of length 2 specifying width ratio of heatmap to FC plot. Default: c(4, 1)
#' @param ... Additional arguments passed to pheatmap
#'
#' @return Invisibly returns a list containing the heatmap object and processed data
#' @export
#' @examples
#' \dontrun{
#' # Create example data
#' set.seed(123)
#' expression <- matrix(rnorm(200), nrow = 20, ncol = 10)
#' rownames(expression) <- paste0("Gene", 1:20)
#' colnames(expression) <- paste0("Sample", 1:10)
#' log2fc <- rnorm(20, mean = 0, sd = 2)
#' groups <- rep(c("Control", "Treatment"), each = 5)
#'
#' # Plot combined heatmap and fold change
#' plot_heatmap_with_fc(
#'   expression_matrix = expression,
#'   log2fc = log2fc,
#'   sample_groups = groups,
#'   fc_threshold = 1.5,
#'   title = "Differential Gene Expression"
#' )
#' }
plot_heatmap_with_fc <- function(expression_matrix,
                                 log2fc,
                                 sample_groups = NULL,
                                 cluster_rows = TRUE,
                                 cluster_cols = TRUE,
                                 scale = c("row", "column", "none"),
                                 show_rownames = TRUE,
                                 show_colnames = TRUE,
                                 heatmap_colors = NULL,
                                 fc_colors = c("#2166AC", "#B2182B"),
                                 fc_threshold = 1,
                                 title = NULL,
                                 gene_labels = NULL,
                                 width_ratio = c(4, 1),
                                 # Statistical annotations
                                 pvalues = NULL,
                                 padj = NULL,
                                 sig_thresholds = c(0.05, 0.01, 0.001),
                                 show_pval_text = FALSE,
                                 # Gene filtering and highlighting
                                 filter_by_fc = FALSE,
                                 top_n = NULL,
                                 highlight_genes = NULL,
                                 highlight_color = "#FF6B35",
                                 gene_categories = NULL,
                                 # Color themes
                                 color_theme = c("default", "viridis", "RdBu", "publication", "colorblind"),
                                 # Enhanced annotations
                                 annotation_row = NULL,
                                 annotation_col_enhanced = NULL,
                                 annotation_colors = NULL,
                                 # Layout options
                                 plot_layout = c("heatmap_fc", "heatmap_fc_volcano"),
                                 ...) {

  # Match arguments
  scale <- match.arg(scale)
  color_theme <- match.arg(color_theme)
  plot_layout <- match.arg(plot_layout)

  # Check required packages
  required_packages <- c("pheatmap", "grid", "gridExtra")
  if (plot_layout == "heatmap_fc_volcano") {
    required_packages <- c(required_packages, "ggplot2")
  }
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' is required. Install with: install.packages('", pkg, "')")
    }
  }

  # Validate inputs
  if (!is.matrix(expression_matrix)) {
    expression_matrix <- as.matrix(expression_matrix)
  }

  if (length(log2fc) != nrow(expression_matrix)) {
    stop("Length of log2fc must match number of rows in expression_matrix")
  }

  # Validate p-values if provided
  if (!is.null(pvalues) && length(pvalues) != nrow(expression_matrix)) {
    stop("Length of pvalues must match number of rows in expression_matrix")
  }
  if (!is.null(padj) && length(padj) != nrow(expression_matrix)) {
    stop("Length of padj must match number of rows in expression_matrix")
  }

  # Set gene labels
  if (is.null(gene_labels)) {
    gene_labels <- rownames(expression_matrix)
    if (is.null(gene_labels)) {
      gene_labels <- paste0("Gene", seq_len(nrow(expression_matrix)))
    }
  }

  # Store original indices for tracking
  original_indices <- seq_len(nrow(expression_matrix))

  # Gene filtering
  if (filter_by_fc || !is.null(top_n)) {
    keep_genes <- rep(TRUE, nrow(expression_matrix))

    # Filter by FC threshold
    if (filter_by_fc) {
      keep_genes <- abs(log2fc) > fc_threshold
    }

    # Select top N up/down regulated genes
    if (!is.null(top_n)) {
      # Sort by absolute log2FC
      fc_order <- order(abs(log2fc), decreasing = TRUE)
      top_indices <- fc_order[1:min(top_n, length(fc_order))]
      keep_genes <- original_indices %in% top_indices
    }

    # Apply filtering
    expression_matrix <- expression_matrix[keep_genes, , drop = FALSE]
    log2fc <- log2fc[keep_genes]
    gene_labels <- gene_labels[keep_genes]
    original_indices <- original_indices[keep_genes]

    if (!is.null(pvalues)) pvalues <- pvalues[keep_genes]
    if (!is.null(padj)) padj <- padj[keep_genes]
    if (!is.null(gene_categories)) gene_categories <- gene_categories[keep_genes]

    message(sprintf("Filtered to %d genes (from %d total)",
                    sum(keep_genes), length(keep_genes)))
  }

  # Auto-hide row names if too many genes
  if (show_rownames && nrow(expression_matrix) > 50) {
    show_rownames <- FALSE
    message("Note: Row names hidden (>50 genes). Set show_rownames=TRUE to override.")
  }

  # Apply color theme
  if (is.null(heatmap_colors)) {
    heatmap_colors <- switch(color_theme,
      "default" = grDevices::colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
      "viridis" = grDevices::colorRampPalette(c("#440154", "#31688E", "#35B779", "#FDE724"))(100),
      "RdBu" = grDevices::colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                                             "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                                             "#4393C3", "#2166AC", "#053061"))(100),
      "publication" = grDevices::colorRampPalette(c("black", "grey90", "white"))(100),
      "colorblind" = grDevices::colorRampPalette(c("#0173B2", "white", "#DE8F05"))(100)
    )
  }

  # Prepare column annotations
  annotation_col <- NULL
  if (!is.null(annotation_col_enhanced)) {
    # Use enhanced annotations if provided
    annotation_col <- annotation_col_enhanced
  } else if (!is.null(sample_groups)) {
    # Fall back to simple sample groups
    if (length(sample_groups) != ncol(expression_matrix)) {
      warning("Length of sample_groups doesn't match number of columns. Ignoring sample_groups.")
    } else {
      annotation_col <- data.frame(Group = sample_groups)
      rownames(annotation_col) <- colnames(expression_matrix)
    }
  }

  # Prepare row annotations
  if (!is.null(annotation_row)) {
    if (nrow(annotation_row) != nrow(expression_matrix)) {
      warning("Row number of annotation_row doesn't match expression_matrix rows. Ignoring annotation_row.")
      annotation_row <- NULL
    }
  } else if (!is.null(gene_categories)) {
    # Create annotation from gene categories
    annotation_row <- data.frame(Category = gene_categories)
    rownames(annotation_row) <- gene_labels
  }

  # Add highlighted genes to row annotation
  if (!is.null(highlight_genes)) {
    highlight_status <- ifelse(gene_labels %in% highlight_genes, "Highlighted", "Other")
    if (!is.null(annotation_row)) {
      annotation_row$Highlight <- highlight_status
    } else {
      annotation_row <- data.frame(Highlight = highlight_status)
      rownames(annotation_row) <- gene_labels
    }

    # Add highlight colors to annotation_colors
    if (is.null(annotation_colors)) {
      annotation_colors <- list()
    }
    if (!"Highlight" %in% names(annotation_colors)) {
      annotation_colors$Highlight <- c("Highlighted" = highlight_color, "Other" = "white")
    }
  }

  # Create the main heatmap
  hm <- pheatmap::pheatmap(
    expression_matrix,
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    scale = scale,
    show_rownames = show_rownames,
    show_colnames = show_colnames,
    color = heatmap_colors,
    annotation_col = annotation_col,
    annotation_row = annotation_row,
    annotation_colors = annotation_colors,
    main = title,
    silent = TRUE,
    ...
  )

  # Get the row order from clustering
  if (cluster_rows) {
    row_order <- hm$tree_row$order
  } else {
    row_order <- seq_len(nrow(expression_matrix))
  }

  # Reorder log2fc to match heatmap
  log2fc_ordered <- log2fc[row_order]
  genes_ordered <- gene_labels[row_order]

  # Reorder p-values if provided
  pval_ordered <- NULL
  padj_ordered <- NULL
  if (!is.null(pvalues)) pval_ordered <- pvalues[row_order]
  if (!is.null(padj)) padj_ordered <- padj[row_order]

  # Determine significance level for each gene
  sig_stars <- rep("", length(log2fc_ordered))
  if (!is.null(padj_ordered) || !is.null(pval_ordered)) {
    pval_to_use <- if (!is.null(padj_ordered)) padj_ordered else pval_ordered
    sig_stars <- ifelse(pval_to_use < sig_thresholds[3], "***",
                 ifelse(pval_to_use < sig_thresholds[2], "**",
                 ifelse(pval_to_use < sig_thresholds[1], "*", "")))
  }

  # Create fold change data frame
  fc_data <- data.frame(
    gene = factor(genes_ordered, levels = genes_ordered),
    log2fc = log2fc_ordered,
    significant = abs(log2fc_ordered) > fc_threshold,
    sig_stars = sig_stars,
    pval = pval_ordered,
    stringsAsFactors = FALSE
  )

  # Create the fold change bar plot
  # Show gene names on FC plot only if they're shown on heatmap
  show_fc_labels <- show_rownames

  # Reverse the factor levels so genes appear in same order as heatmap (top to bottom)
  fc_data$gene <- factor(fc_data$gene, levels = rev(levels(fc_data$gene)))

  # Build the base plot
  fc_plot <- ggplot2::ggplot(fc_data, ggplot2::aes(x = gene, y = log2fc)) +
    ggplot2::geom_hline(yintercept = 0, color = "gray40", linewidth = 0.5) +
    ggplot2::geom_hline(yintercept = c(-fc_threshold, fc_threshold),
                       color = "gray70", linetype = "dashed", linewidth = 0.3) +
    ggplot2::geom_bar(
      stat = "identity",
      ggplot2::aes(fill = log2fc > 0),
      width = 0.8,
      show.legend = FALSE
    ) +
    ggplot2::scale_fill_manual(values = fc_colors) +
    ggplot2::coord_flip() +
    ggplot2::scale_x_discrete(drop = FALSE) +
    ggplot2::labs(
      x = NULL,
      y = "log2FC",
      title = NULL
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      axis.text.y = if (show_fc_labels) ggplot2::element_text(size = 9, hjust = 1) else ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(t = 5, r = 10, b = 5, l = 5),
      axis.title.x = ggplot2::element_text(size = 9)
    )

  # Add numeric FC values
  fc_plot <- fc_plot +
    ggplot2::geom_text(
      ggplot2::aes(
        label = sprintf("%.2f", log2fc),
        hjust = ifelse(log2fc > 0, -0.1, 1.1)
      ),
      size = 3,
      color = "black"
    )

  # Add significance stars if available
  if (any(fc_data$sig_stars != "")) {
    fc_plot <- fc_plot +
      ggplot2::geom_text(
        data = fc_data[fc_data$sig_stars != "", ],
        ggplot2::aes(
          label = sig_stars,
          y = 0,
          hjust = 0.5
        ),
        size = 4,
        color = "black",
        fontface = "bold"
      )
  }

  # Add p-value text if requested
  if (show_pval_text && !is.null(fc_data$pval)) {
    fc_plot <- fc_plot +
      ggplot2::geom_text(
        ggplot2::aes(
          label = ifelse(!is.na(pval), sprintf("p=%.3f", pval), ""),
          hjust = ifelse(log2fc > 0, -0.1, 1.1),
          y = log2fc
        ),
        size = 2.5,
        color = "gray30",
        vjust = -0.5
      )
  }

  # Convert fold change plot to grob
  fc_grob <- ggplot2::ggplotGrob(fc_plot)

  # Create volcano plot if requested
  volcano_plot <- NULL
  volcano_grob <- NULL
  if (plot_layout == "heatmap_fc_volcano") {
    if (is.null(pval_ordered) && is.null(padj_ordered)) {
      warning("Cannot create volcano plot without p-values. Falling back to heatmap_fc layout.")
      plot_layout <- "heatmap_fc"
    } else {
      # Prepare volcano data
      volcano_data <- data.frame(
        gene = genes_ordered,
        log2fc = log2fc_ordered,
        neglog10p = -log10(if (!is.null(padj_ordered)) padj_ordered else pval_ordered),
        significant = (abs(log2fc_ordered) > fc_threshold) &
                      ((if (!is.null(padj_ordered)) padj_ordered else pval_ordered) < sig_thresholds[1]),
        highlighted = genes_ordered %in% highlight_genes,
        stringsAsFactors = FALSE
      )

      # Remove infinite values
      volcano_data <- volcano_data[is.finite(volcano_data$neglog10p), ]

      # Create volcano plot
      volcano_plot <- ggplot2::ggplot(volcano_data,
                                      ggplot2::aes(x = log2fc, y = neglog10p)) +
        ggplot2::geom_point(ggplot2::aes(color = significant), alpha = 0.6, size = 2) +
        ggplot2::geom_vline(xintercept = c(-fc_threshold, fc_threshold),
                           linetype = "dashed", color = "gray50") +
        ggplot2::geom_hline(yintercept = -log10(sig_thresholds[1]),
                           linetype = "dashed", color = "gray50") +
        ggplot2::scale_color_manual(values = c("FALSE" = "gray60", "TRUE" = "#D73027"),
                                   name = "Significant") +
        ggplot2::labs(
          x = "log2 Fold Change",
          y = "-log10(p-value)",
          title = "Volcano Plot"
        ) +
        ggplot2::theme_minimal(base_size = 10) +
        ggplot2::theme(
          legend.position = "bottom",
          plot.margin = ggplot2::margin(t = 5, r = 5, b = 5, l = 5)
        )

      # Highlight specific genes if provided
      if (!is.null(highlight_genes) && any(volcano_data$highlighted)) {
        volcano_plot <- volcano_plot +
          ggplot2::geom_point(data = volcano_data[volcano_data$highlighted, ],
                             color = highlight_color, size = 3, shape = 21,
                             fill = highlight_color, stroke = 1.5)
      }

      volcano_grob <- ggplot2::ggplotGrob(volcano_plot)
    }
  }

  # Extract heatmap grob
  hm_grob <- hm$gtable

  # Combine plots based on layout
  grid::grid.newpage()

  # Match the plot area heights
  # Find the panel in heatmap gtable
  hm_panel_pos <- grep("panel", hm_grob$layout$name)
  if (length(hm_panel_pos) > 0) {
    panel_row <- hm_grob$layout$t[hm_panel_pos[1]]
    panel_height <- hm_grob$heights[panel_row]

    # Set fc_grob panel to match heatmap panel height
    fc_panel_pos <- grep("panel", fc_grob$layout$name)
    if (length(fc_panel_pos) > 0) {
      fc_panel_row <- fc_grob$layout$t[fc_panel_pos[1]]
      fc_grob$heights[fc_panel_row] <- panel_height
    }
  }

  # Create layout based on plot_layout parameter
  if (plot_layout == "heatmap_fc_volcano" && !is.null(volcano_grob)) {
    # Three panel layout: heatmap, fc, volcano
    grid::pushViewport(grid::viewport(
      layout = grid::grid.layout(1, 3,
        widths = grid::unit(c(4, 1, 2), c("null", "null", "null")))
    ))

    # Draw heatmap on the left
    grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
    grid::grid.draw(hm_grob)
    grid::popViewport()

    # Draw fold change plot in the middle
    grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
    grid::grid.draw(fc_grob)
    grid::popViewport()

    # Draw volcano plot on the right
    grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 3))
    grid::grid.draw(volcano_grob)
    grid::popViewport()
  } else {
    # Two panel layout: heatmap and fc (default)
    grid::pushViewport(grid::viewport(
      layout = grid::grid.layout(1, 2,
        widths = grid::unit(width_ratio, c("null", "null")))
    ))

    # Draw heatmap on the left
    grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
    grid::grid.draw(hm_grob)
    grid::popViewport()

    # Draw fold change plot on the right
    grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
    grid::grid.draw(fc_grob)
    grid::popViewport()
  }

  # Return results invisibly
  result <- list(
    heatmap = hm,
    fc_plot = fc_plot,
    fc_data = fc_data,
    row_order = row_order,
    expression_matrix = expression_matrix,
    log2fc = log2fc,
    volcano_plot = volcano_plot,
    pvalues = pval_ordered,
    padj = padj_ordered
  )

  invisible(result)
}
