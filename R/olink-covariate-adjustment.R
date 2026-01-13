#' Adjust Covariates for Olink Count Matrix Data
#'
#' Performs covariate adjustment on Olink protein expression count matrix data
#' using linear regression. This function removes the effects of specified
#' covariates (e.g., age, sex, batch effects) from the protein expression data
#' while preserving the biological signal of interest.
#'
#' @param count_matrix A numeric matrix or data.frame with proteins as rows and
#'   samples as columns. Row names should be protein IDs and column names should
#'   be sample IDs.
#' @param metadata A data.frame with sample metadata. Must contain a column matching
#'   sample IDs that correspond to column names in count_matrix. Row names or a
#'   specified ID column should match count_matrix column names.
#' @param covariates Character vector of covariate names to adjust for. These should
#'   be column names in the metadata data.frame (e.g., c("age", "sex", "batch")).
#' @param sample_id_col Name of the column in metadata containing sample IDs. If NULL,
#'   uses row names of metadata (default: NULL).
#' @param method Adjustment method to use. Options are:
#'   \itemize{
#'     \item "residuals" - Returns residuals from linear model (default)
#'     \item "adjusted" - Returns adjusted values (original - covariate effects + grand mean)
#'   }
#' @param log_transform Logical indicating whether to log-transform the data before
#'   adjustment (default: FALSE). If TRUE, uses log2(x + 1) transformation.
#' @param return_log Logical indicating whether to return log-transformed values when
#'   log_transform = TRUE (default: TRUE). If FALSE, back-transforms to original scale.
#' @param center_covariates Logical indicating whether to center numeric covariates
#'   before adjustment (default: TRUE).
#' @param verbose Logical indicating whether to print progress messages (default: TRUE).
#'
#' @return A list containing:
#'   \itemize{
#'     \item adjusted_matrix - The covariate-adjusted count matrix
#'     \item model_info - List of model information for each protein
#'     \item covariates_used - Character vector of covariates that were adjusted
#'     \item samples_used - Character vector of samples included in adjustment
#'     \item proteins_failed - Character vector of proteins where adjustment failed
#'   }
#'
#' @export
#' @examples
#' \dontrun{
#' # Create example data
#' count_matrix <- matrix(rnorm(1000, mean = 100, sd = 20),
#'                        nrow = 20, ncol = 50)
#' rownames(count_matrix) <- paste0("Protein_", 1:20)
#' colnames(count_matrix) <- paste0("Sample_", 1:50)
#'
#' metadata <- data.frame(
#'   sample_id = paste0("Sample_", 1:50),
#'   age = rnorm(50, mean = 45, sd = 10),
#'   sex = sample(c("M", "F"), 50, replace = TRUE),
#'   batch = sample(c("Batch1", "Batch2", "Batch3"), 50, replace = TRUE)
#' )
#'
#' # Adjust for age, sex, and batch
#' result <- adjust_olink_covariates(
#'   count_matrix = count_matrix,
#'   metadata = metadata,
#'   covariates = c("age", "sex", "batch"),
#'   sample_id_col = "sample_id",
#'   method = "adjusted"
#' )
#'
#' adjusted_data <- result$adjusted_matrix
#' }
adjust_olink_covariates <- function(count_matrix,
                                     metadata,
                                     covariates,
                                     sample_id_col = NULL,
                                     method = c("residuals", "adjusted"),
                                     log_transform = FALSE,
                                     return_log = TRUE,
                                     center_covariates = TRUE,
                                     verbose = TRUE) {

  # Input validation
  method <- match.arg(method)

  if (!is.matrix(count_matrix) && !is.data.frame(count_matrix)) {
    stop("count_matrix must be a matrix or data.frame")
  }

  if (!is.data.frame(metadata)) {
    stop("metadata must be a data.frame")
  }

  if (length(covariates) == 0) {
    stop("At least one covariate must be specified")
  }

  # Convert count_matrix to matrix if data.frame
  if (is.data.frame(count_matrix)) {
    count_matrix <- as.matrix(count_matrix)
  }

  # Get sample IDs from metadata
  if (is.null(sample_id_col)) {
    sample_ids <- rownames(metadata)
    if (is.null(sample_ids)) {
      stop("metadata must have row names or sample_id_col must be specified")
    }
  } else {
    if (!sample_id_col %in% colnames(metadata)) {
      stop(paste("Column", sample_id_col, "not found in metadata"))
    }
    sample_ids <- metadata[[sample_id_col]]
    rownames(metadata) <- sample_ids
  }

  # Check that covariates exist in metadata
  missing_covariates <- setdiff(covariates, colnames(metadata))
  if (length(missing_covariates) > 0) {
    stop(paste("Covariates not found in metadata:",
               paste(missing_covariates, collapse = ", ")))
  }

  # Match samples between count_matrix and metadata
  common_samples <- intersect(colnames(count_matrix), sample_ids)
  if (length(common_samples) == 0) {
    stop("No common samples found between count_matrix and metadata")
  }

  if (verbose) {
    message(paste0("Processing ", nrow(count_matrix), " proteins across ",
                   length(common_samples), " samples"))
    message(paste0("Adjusting for covariates: ", paste(covariates, collapse = ", ")))
  }

  # Subset and align data
  count_matrix <- count_matrix[, common_samples, drop = FALSE]
  metadata <- metadata[common_samples, , drop = FALSE]

  # Log transformation if requested
  if (log_transform) {
    if (verbose) message("Applying log2(x + 1) transformation...")
    count_matrix <- log2(count_matrix + 1)
  }

  # Prepare covariate data
  covariate_data <- metadata[, covariates, drop = FALSE]

  # Center numeric covariates if requested
  if (center_covariates) {
    numeric_cols <- sapply(covariate_data, is.numeric)
    if (any(numeric_cols)) {
      covariate_data[, numeric_cols] <- scale(covariate_data[, numeric_cols],
                                                center = TRUE, scale = FALSE)
      if (verbose) message("Centered numeric covariates")
    }
  }

  # Initialize output
  adjusted_matrix <- matrix(NA, nrow = nrow(count_matrix),
                            ncol = ncol(count_matrix))
  rownames(adjusted_matrix) <- rownames(count_matrix)
  colnames(adjusted_matrix) <- colnames(count_matrix)

  model_info <- list()
  proteins_failed <- character(0)

  # Progress bar
  if (verbose) {
    message("Adjusting protein expression values...")
    pb <- txtProgressBar(min = 0, max = nrow(count_matrix), style = 3)
  }

  # Adjust each protein
  for (i in seq_len(nrow(count_matrix))) {
    protein_name <- rownames(count_matrix)[i]
    protein_values <- count_matrix[i, ]

    # Create data frame for modeling
    model_data <- data.frame(
      expression = protein_values,
      covariate_data,
      stringsAsFactors = FALSE
    )

    # Remove samples with missing data
    complete_cases <- complete.cases(model_data)
    if (sum(complete_cases) < 3) {
      proteins_failed <- c(proteins_failed, protein_name)
      if (verbose) setTxtProgressBar(pb, i)
      next
    }

    model_data <- model_data[complete_cases, , drop = FALSE]

    # Fit linear model
    tryCatch({
      formula_str <- paste("expression ~", paste(covariates, collapse = " + "))
      fit <- lm(as.formula(formula_str), data = model_data)

      # Extract adjusted values based on method
      if (method == "residuals") {
        adjusted_values <- residuals(fit)
      } else if (method == "adjusted") {
        # Adjusted = Original - Covariate Effects + Grand Mean
        covariate_effects <- predict(fit) - mean(model_data$expression)
        adjusted_values <- model_data$expression - covariate_effects
      }

      # Store adjusted values
      adjusted_matrix[i, complete_cases] <- adjusted_values

      # Store model info
      model_info[[protein_name]] <- list(
        r_squared = summary(fit)$r.squared,
        adj_r_squared = summary(fit)$adj.r.squared,
        n_samples = sum(complete_cases),
        coefficients = coef(fit)
      )

    }, error = function(e) {
      proteins_failed <<- c(proteins_failed, protein_name)
    })

    if (verbose) setTxtProgressBar(pb, i)
  }

  if (verbose) {
    close(pb)
    message(paste0("Successfully adjusted ", nrow(count_matrix) - length(proteins_failed),
                   " proteins"))
    if (length(proteins_failed) > 0) {
      message(paste0("Failed to adjust ", length(proteins_failed), " proteins"))
    }
  }

  # Back-transform if needed
  if (log_transform && !return_log) {
    if (verbose) message("Back-transforming to original scale...")
    adjusted_matrix <- 2^adjusted_matrix - 1
    adjusted_matrix[adjusted_matrix < 0] <- 0
  }

  # Return results
  return(list(
    adjusted_matrix = adjusted_matrix,
    model_info = model_info,
    covariates_used = covariates,
    samples_used = common_samples,
    proteins_failed = proteins_failed
  ))
}


#' Plot Covariate Adjustment Results
#'
#' Creates diagnostic plots to visualize the effects of covariate adjustment
#' on Olink protein expression data.
#'
#' @param original_matrix Original count matrix before adjustment
#' @param adjusted_matrix Adjusted count matrix after adjustment
#' @param metadata Sample metadata data.frame
#' @param covariate Name of a specific covariate to visualize
#' @param proteins Character vector of protein names to plot (default: first 6)
#'
#' @return A ggplot object showing before/after comparison
#' @export
#' @examples
#' \dontrun{
#' plot_adjustment_results(
#'   original_matrix = original_data,
#'   adjusted_matrix = result$adjusted_matrix,
#'   metadata = metadata,
#'   covariate = "age"
#' )
#' }
plot_adjustment_results <- function(original_matrix,
                                     adjusted_matrix,
                                     metadata,
                                     covariate,
                                     proteins = NULL) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting. Please install it.")
  }

  if (!requireNamespace("tidyr", quietly = TRUE)) {
    stop("Package 'tidyr' is required for plotting. Please install it.")
  }

  # Select proteins to plot
  if (is.null(proteins)) {
    proteins <- head(rownames(original_matrix), 6)
  }

  # Prepare data for plotting
  plot_data <- data.frame()

  for (protein in proteins) {
    original_vals <- original_matrix[protein, ]
    adjusted_vals <- adjusted_matrix[protein, ]

    common_samples <- intersect(names(original_vals),
                                 names(adjusted_vals))

    temp_df <- data.frame(
      sample = common_samples,
      protein = protein,
      original = original_vals[common_samples],
      adjusted = adjusted_vals[common_samples],
      covariate_value = metadata[common_samples, covariate]
    )

    plot_data <- rbind(plot_data, temp_df)
  }

  # Reshape for plotting
  plot_data_long <- tidyr::pivot_longer(
    plot_data,
    cols = c("original", "adjusted"),
    names_to = "type",
    values_to = "expression"
  )

  # Create plot
  p <- ggplot2::ggplot(plot_data_long,
                       ggplot2::aes(x = covariate_value, y = expression,
                                    color = type)) +
    ggplot2::geom_point(alpha = 0.5) +
    ggplot2::geom_smooth(method = "lm", se = TRUE) +
    ggplot2::facet_wrap(~ protein, scales = "free_y", ncol = 2) +
    ggplot2::labs(
      title = paste("Covariate Adjustment Results"),
      subtitle = paste("Covariate:", covariate),
      x = covariate,
      y = "Expression",
      color = "Data Type"
    ) +
    ggplot2::theme_bw() +
    ggplot2::scale_color_manual(
      values = c("original" = "#E74C3C", "adjusted" = "#3498DB"),
      labels = c("Original", "Adjusted")
    )

  return(p)
}
