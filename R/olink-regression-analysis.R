#' Perform Univariate and Multivariate Regression Analysis on Olink Data
#'
#' Conducts comprehensive regression analysis (linear or logistic) on Olink protein
#' expression data. Performs both univariate (one protein at a time) and multivariate
#' (adjusted for covariates) analyses. Returns beta coefficients for continuous outcomes
#' or odds ratios for binary outcomes.
#'
#' @param count_matrix A numeric matrix or data.frame with proteins as rows and
#'   samples as columns. Row names should be protein IDs.
#' @param metadata A data.frame with sample metadata including outcome and covariates.
#' @param outcome Character string specifying the outcome variable name in metadata.
#' @param outcome_type Type of outcome: "binary" for logistic regression (returns OR),
#'   "continuous" for linear regression (returns beta). Auto-detected if NULL.
#' @param covariates Character vector of covariate names for multivariate adjustment.
#'   If NULL, only univariate analysis is performed.
#' @param sample_id_col Name of the column in metadata containing sample IDs.
#'   If NULL, uses row names of metadata (default: NULL).
#' @param log_transform Logical indicating whether to log-transform protein data
#'   before analysis (default: FALSE). Uses log2(x + 1) transformation.
#' @param conf_level Confidence level for confidence intervals (default: 0.95).
#' @param p_adjust_method Method for multiple testing correction. Options: "none",
#'   "bonferroni", "holm", "hochberg", "hommel", "BH", "BY", "fdr" (default: "BH").
#' @param min_samples Minimum number of complete samples required per protein
#'   (default: 10).
#' @param verbose Logical indicating whether to print progress messages (default: TRUE).
#'
#' @return A list containing:
#'   \itemize{
#'     \item univariate_results - Data frame with univariate analysis results
#'     \item multivariate_results - Data frame with multivariate analysis results (if covariates provided)
#'     \item summary_stats - Summary statistics of the analysis
#'     \item model_type - Type of regression model used ("logistic" or "linear")
#'   }
#'
#' @details
#' For binary outcomes (logistic regression):
#'   - Returns odds ratios (OR) with 95% CI
#'   - Outcome should be coded as 0/1 or factor
#'
#' For continuous outcomes (linear regression):
#'   - Returns beta coefficients with 95% CI
#'   - Outcome should be numeric
#'
#' @export
#' @examples
#' \dontrun{
#' # Binary outcome example (case-control study)
#' count_matrix <- matrix(rnorm(1000), nrow = 50, ncol = 20)
#' rownames(count_matrix) <- paste0("Protein_", 1:50)
#' colnames(count_matrix) <- paste0("Sample_", 1:20)
#'
#' metadata <- data.frame(
#'   sample_id = paste0("Sample_", 1:20),
#'   disease_status = sample(c(0, 1), 20, replace = TRUE),
#'   age = rnorm(20, 50, 10),
#'   sex = sample(c("M", "F"), 20, replace = TRUE),
#'   bmi = rnorm(20, 25, 4)
#' )
#'
#' # Run analysis
#' results <- olink_regression_analysis(
#'   count_matrix = count_matrix,
#'   metadata = metadata,
#'   outcome = "disease_status",
#'   outcome_type = "binary",
#'   covariates = c("age", "sex", "bmi"),
#'   sample_id_col = "sample_id"
#' )
#'
#' # View top results
#' head(results$univariate_results)
#' head(results$multivariate_results)
#'
#' # Continuous outcome example
#' metadata$severity_score <- rnorm(20, 50, 15)
#'
#' results_cont <- olink_regression_analysis(
#'   count_matrix = count_matrix,
#'   metadata = metadata,
#'   outcome = "severity_score",
#'   outcome_type = "continuous",
#'   covariates = c("age", "sex")
#' )
#' }
olink_regression_analysis <- function(count_matrix,
                                       metadata,
                                       outcome,
                                       outcome_type = NULL,
                                       covariates = NULL,
                                       sample_id_col = NULL,
                                       log_transform = FALSE,
                                       conf_level = 0.95,
                                       p_adjust_method = "BH",
                                       min_samples = 10,
                                       verbose = TRUE) {

  # Input validation
  if (!is.matrix(count_matrix) && !is.data.frame(count_matrix)) {
    stop("count_matrix must be a matrix or data.frame")
  }

  if (!is.data.frame(metadata)) {
    stop("metadata must be a data.frame")
  }

  if (!outcome %in% colnames(metadata)) {
    stop(paste("Outcome variable", outcome, "not found in metadata"))
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

  # Match samples between count_matrix and metadata
  common_samples <- intersect(colnames(count_matrix), sample_ids)
  if (length(common_samples) == 0) {
    stop("No common samples found between count_matrix and metadata")
  }

  if (length(common_samples) < min_samples) {
    stop(paste("Insufficient samples:", length(common_samples),
               "< minimum required:", min_samples))
  }

  # Subset and align data
  count_matrix <- count_matrix[, common_samples, drop = FALSE]
  metadata <- metadata[common_samples, , drop = FALSE]

  # Auto-detect outcome type if not specified
  if (is.null(outcome_type)) {
    outcome_vals <- metadata[[outcome]]
    unique_vals <- length(unique(na.omit(outcome_vals)))

    if (is.factor(outcome_vals) || unique_vals == 2) {
      outcome_type <- "binary"
    } else if (is.numeric(outcome_vals) && unique_vals > 2) {
      outcome_type <- "continuous"
    } else {
      stop("Cannot auto-detect outcome type. Please specify outcome_type.")
    }

    if (verbose) {
      message(paste("Auto-detected outcome type:", outcome_type))
    }
  }

  # Validate outcome type
  outcome_type <- match.arg(outcome_type, c("binary", "continuous"))

  # Prepare outcome variable
  if (outcome_type == "binary") {
    metadata[[outcome]] <- as.numeric(as.factor(metadata[[outcome]])) - 1
    model_family <- "binomial"
    model_type <- "logistic"
  } else {
    metadata[[outcome]] <- as.numeric(metadata[[outcome]])
    model_family <- "gaussian"
    model_type <- "linear"
  }

  # Check covariates
  if (!is.null(covariates)) {
    missing_covariates <- setdiff(covariates, colnames(metadata))
    if (length(missing_covariates) > 0) {
      stop(paste("Covariates not found in metadata:",
                 paste(missing_covariates, collapse = ", ")))
    }
  }

  # Log transformation if requested
  if (log_transform) {
    if (verbose) message("Applying log2(x + 1) transformation...")
    count_matrix <- log2(count_matrix + 1)
  }

  if (verbose) {
    message(paste0("Analyzing ", nrow(count_matrix), " proteins across ",
                   length(common_samples), " samples"))
    message(paste0("Model type: ", model_type, " regression"))
    message(paste0("Outcome: ", outcome))
    if (!is.null(covariates)) {
      message(paste0("Covariates: ", paste(covariates, collapse = ", ")))
    }
  }

  # Perform univariate analysis
  if (verbose) message("\nPerforming univariate analysis...")
  univariate_results <- perform_univariate_analysis(
    count_matrix = count_matrix,
    metadata = metadata,
    outcome = outcome,
    model_family = model_family,
    model_type = model_type,
    conf_level = conf_level,
    min_samples = min_samples,
    verbose = verbose
  )

  # Apply multiple testing correction
  univariate_results$p_adjusted <- p.adjust(univariate_results$p_value,
                                              method = p_adjust_method)

  # Perform multivariate analysis if covariates provided
  multivariate_results <- NULL
  if (!is.null(covariates)) {
    if (verbose) message("\nPerforming multivariate analysis...")
    multivariate_results <- perform_multivariate_analysis(
      count_matrix = count_matrix,
      metadata = metadata,
      outcome = outcome,
      covariates = covariates,
      model_family = model_family,
      model_type = model_type,
      conf_level = conf_level,
      min_samples = min_samples,
      verbose = verbose
    )

    # Apply multiple testing correction
    multivariate_results$p_adjusted <- p.adjust(multivariate_results$p_value,
                                                  method = p_adjust_method)
  }

  # Generate summary statistics
  summary_stats <- list(
    n_proteins_total = nrow(count_matrix),
    n_proteins_analyzed_univariate = nrow(univariate_results),
    n_samples = length(common_samples),
    outcome_variable = outcome,
    model_type = model_type,
    p_adjust_method = p_adjust_method,
    n_significant_univariate_raw = sum(univariate_results$p_value < 0.05, na.rm = TRUE),
    n_significant_univariate_adj = sum(univariate_results$p_adjusted < 0.05, na.rm = TRUE)
  )

  if (!is.null(multivariate_results)) {
    summary_stats$n_proteins_analyzed_multivariate <- nrow(multivariate_results)
    summary_stats$n_significant_multivariate_raw <- sum(multivariate_results$p_value < 0.05, na.rm = TRUE)
    summary_stats$n_significant_multivariate_adj <- sum(multivariate_results$p_adjusted < 0.05, na.rm = TRUE)
    summary_stats$covariates_used <- covariates
  }

  if (verbose) {
    message("\n=== Analysis Summary ===")
    message(paste("Total proteins analyzed:", summary_stats$n_proteins_analyzed_univariate))
    message(paste("Significant (p < 0.05, univariate):", summary_stats$n_significant_univariate_raw))
    message(paste("Significant (adjusted p < 0.05, univariate):", summary_stats$n_significant_univariate_adj))
    if (!is.null(multivariate_results)) {
      message(paste("Significant (p < 0.05, multivariate):", summary_stats$n_significant_multivariate_raw))
      message(paste("Significant (adjusted p < 0.05, multivariate):", summary_stats$n_significant_multivariate_adj))
    }
  }

  # Return results
  return(list(
    univariate_results = univariate_results,
    multivariate_results = multivariate_results,
    summary_stats = summary_stats,
    model_type = model_type
  ))
}


#' Perform Univariate Analysis
#'
#' Internal function to perform univariate regression for each protein.
#'
#' @keywords internal
perform_univariate_analysis <- function(count_matrix, metadata, outcome,
                                         model_family, model_type,
                                         conf_level, min_samples, verbose) {

  results_list <- list()
  n_proteins <- nrow(count_matrix)

  if (verbose) {
    pb <- txtProgressBar(min = 0, max = n_proteins, style = 3)
  }

  for (i in seq_len(n_proteins)) {
    protein_name <- rownames(count_matrix)[i]
    protein_values <- count_matrix[i, ]

    # Create analysis data frame
    analysis_data <- data.frame(
      outcome = metadata[[outcome]],
      protein = protein_values,
      stringsAsFactors = FALSE
    )

    # Remove missing values
    analysis_data <- analysis_data[complete.cases(analysis_data), , drop = FALSE]

    # Check minimum samples
    if (nrow(analysis_data) < min_samples) {
      if (verbose) setTxtProgressBar(pb, i)
      next
    }

    # Fit model
    tryCatch({
      if (model_family == "binomial") {
        fit <- glm(outcome ~ protein, data = analysis_data, family = binomial())
      } else {
        fit <- lm(outcome ~ protein, data = analysis_data)
      }

      # Extract results
      coef_summary <- summary(fit)$coefficients
      protein_coef <- coef_summary["protein", "Estimate"]
      protein_se <- coef_summary["protein", "Std. Error"]

      if (model_family == "binomial") {
        protein_pval <- coef_summary["protein", "Pr(>|z|)"]
      } else {
        protein_pval <- coef_summary["protein", "Pr(>|t|)"]
      }

      # Calculate confidence intervals
      ci <- confint.default(fit, parm = "protein", level = conf_level)

      # For logistic regression, convert to OR
      if (model_type == "logistic") {
        or <- exp(protein_coef)
        or_lower <- exp(ci[1])
        or_upper <- exp(ci[2])

        results_list[[protein_name]] <- data.frame(
          protein = protein_name,
          OR = or,
          OR_lower = or_lower,
          OR_upper = or_upper,
          log_OR = protein_coef,
          SE = protein_se,
          p_value = protein_pval,
          n_samples = nrow(analysis_data),
          stringsAsFactors = FALSE
        )
      } else {
        # Linear regression - return beta
        results_list[[protein_name]] <- data.frame(
          protein = protein_name,
          beta = protein_coef,
          beta_lower = ci[1],
          beta_upper = ci[2],
          SE = protein_se,
          p_value = protein_pval,
          n_samples = nrow(analysis_data),
          stringsAsFactors = FALSE
        )
      }

    }, error = function(e) {
      # Skip proteins that fail
    })

    if (verbose) setTxtProgressBar(pb, i)
  }

  if (verbose) close(pb)

  # Combine results
  if (length(results_list) == 0) {
    stop("No proteins were successfully analyzed")
  }

  results_df <- do.call(rbind, results_list)
  rownames(results_df) <- NULL

  # Sort by p-value
  results_df <- results_df[order(results_df$p_value), ]

  return(results_df)
}


#' Perform Multivariate Analysis
#'
#' Internal function to perform multivariate regression for each protein
#' adjusted for covariates.
#'
#' @keywords internal
perform_multivariate_analysis <- function(count_matrix, metadata, outcome,
                                           covariates, model_family, model_type,
                                           conf_level, min_samples, verbose) {

  results_list <- list()
  n_proteins <- nrow(count_matrix)

  if (verbose) {
    pb <- txtProgressBar(min = 0, max = n_proteins, style = 3)
  }

  for (i in seq_len(n_proteins)) {
    protein_name <- rownames(count_matrix)[i]
    protein_values <- count_matrix[i, ]

    # Create analysis data frame
    analysis_data <- data.frame(
      outcome = metadata[[outcome]],
      protein = protein_values,
      metadata[, covariates, drop = FALSE],
      stringsAsFactors = FALSE
    )

    # Remove missing values
    analysis_data <- analysis_data[complete.cases(analysis_data), , drop = FALSE]

    # Check minimum samples
    if (nrow(analysis_data) < min_samples) {
      if (verbose) setTxtProgressBar(pb, i)
      next
    }

    # Build formula
    formula_str <- paste("outcome ~ protein +", paste(covariates, collapse = " + "))

    # Fit model
    tryCatch({
      if (model_family == "binomial") {
        fit <- glm(as.formula(formula_str), data = analysis_data, family = binomial())
      } else {
        fit <- lm(as.formula(formula_str), data = analysis_data)
      }

      # Extract results for protein term
      coef_summary <- summary(fit)$coefficients
      protein_coef <- coef_summary["protein", "Estimate"]
      protein_se <- coef_summary["protein", "Std. Error"]

      if (model_family == "binomial") {
        protein_pval <- coef_summary["protein", "Pr(>|z|)"]
      } else {
        protein_pval <- coef_summary["protein", "Pr(>|t|)"]
      }

      # Calculate confidence intervals
      ci <- confint.default(fit, parm = "protein", level = conf_level)

      # For logistic regression, convert to OR
      if (model_type == "logistic") {
        or <- exp(protein_coef)
        or_lower <- exp(ci[1])
        or_upper <- exp(ci[2])

        results_list[[protein_name]] <- data.frame(
          protein = protein_name,
          OR = or,
          OR_lower = or_lower,
          OR_upper = or_upper,
          log_OR = protein_coef,
          SE = protein_se,
          p_value = protein_pval,
          n_samples = nrow(analysis_data),
          stringsAsFactors = FALSE
        )
      } else {
        # Linear regression - return beta
        results_list[[protein_name]] <- data.frame(
          protein = protein_name,
          beta = protein_coef,
          beta_lower = ci[1],
          beta_upper = ci[2],
          SE = protein_se,
          p_value = protein_pval,
          n_samples = nrow(analysis_data),
          stringsAsFactors = FALSE
        )
      }

    }, error = function(e) {
      # Skip proteins that fail
    })

    if (verbose) setTxtProgressBar(pb, i)
  }

  if (verbose) close(pb)

  # Combine results
  if (length(results_list) == 0) {
    stop("No proteins were successfully analyzed in multivariate models")
  }

  results_df <- do.call(rbind, results_list)
  rownames(results_df) <- NULL

  # Sort by p-value
  results_df <- results_df[order(results_df$p_value), ]

  return(results_df)
}


#' Create Forest Plot for Regression Results
#'
#' Generates a forest plot to visualize odds ratios or beta coefficients
#' from Olink regression analysis.
#'
#' @param results Results object from olink_regression_analysis()
#' @param analysis_type Type of results to plot: "univariate" or "multivariate"
#' @param top_n Number of top proteins to display (default: 20)
#' @param sort_by Variable to sort by: "p_value" or "effect" (default: "p_value")
#' @param show_pval Logical, whether to show p-values on plot (default: TRUE)
#' @param title Plot title (optional)
#'
#' @return A ggplot object
#' @export
#' @examples
#' \dontrun{
#' plot_forest_results(results, analysis_type = "univariate", top_n = 15)
#' }
plot_forest_results <- function(results,
                                 analysis_type = c("univariate", "multivariate"),
                                 top_n = 20,
                                 sort_by = c("p_value", "effect"),
                                 show_pval = TRUE,
                                 title = NULL) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting. Please install it.")
  }

  analysis_type <- match.arg(analysis_type)
  sort_by <- match.arg(sort_by)

  # Get appropriate results
  if (analysis_type == "univariate") {
    plot_data <- results$univariate_results
  } else {
    if (is.null(results$multivariate_results)) {
      stop("No multivariate results available")
    }
    plot_data <- results$multivariate_results
  }

  # Limit to top N
  plot_data <- head(plot_data, top_n)

  # Determine effect column and label
  model_type <- results$model_type
  if (model_type == "logistic") {
    effect_col <- "OR"
    lower_col <- "OR_lower"
    upper_col <- "OR_upper"
    xlab <- "Odds Ratio (95% CI)"
    null_value <- 1
  } else {
    effect_col <- "beta"
    lower_col <- "beta_lower"
    upper_col <- "beta_upper"
    xlab <- "Beta Coefficient (95% CI)"
    null_value <- 0
  }

  # Sort data
  if (sort_by == "p_value") {
    plot_data <- plot_data[order(plot_data$p_value, decreasing = TRUE), ]
  } else {
    plot_data <- plot_data[order(abs(plot_data[[effect_col]]), decreasing = FALSE), ]
  }

  # Create factor for ordering
  plot_data$protein <- factor(plot_data$protein, levels = plot_data$protein)

  # Create p-value labels
  if (show_pval) {
    plot_data$pval_label <- ifelse(plot_data$p_value < 0.001,
                                    "p < 0.001",
                                    sprintf("p = %.3f", plot_data$p_value))
  }

  # Create plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data[[effect_col]], y = protein)) +
    ggplot2::geom_vline(xintercept = null_value, linetype = "dashed", color = "gray50") +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = .data[[lower_col]],
                                          xmax = .data[[upper_col]]),
                            height = 0.2, color = "gray30") +
    ggplot2::geom_point(ggplot2::aes(color = p_value < 0.05), size = 3) +
    ggplot2::scale_color_manual(values = c("TRUE" = "#E74C3C", "FALSE" = "#95A5A6"),
                                labels = c("TRUE" = "p < 0.05", "FALSE" = "p ≥ 0.05"),
                                name = "Significance") +
    ggplot2::labs(
      x = xlab,
      y = "Protein",
      title = title %||% paste(tools::toTitleCase(analysis_type), "Analysis Results")
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )

  # Add p-value labels if requested
  if (show_pval) {
    p <- p + ggplot2::geom_text(ggplot2::aes(label = pval_label),
                                hjust = -0.1, size = 3)
  }

  return(p)
}
