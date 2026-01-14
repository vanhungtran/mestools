# Test file for Olink regression analysis functions

test_that("olink_regression_analysis works with binary outcome (logistic)", {
  set.seed(123)

  # Create test data with binary outcome
  count_matrix <- matrix(rnorm(200, mean = 100, sd = 20),
                         nrow = 10, ncol = 20)
  rownames(count_matrix) <- paste0("Protein_", 1:10)
  colnames(count_matrix) <- paste0("Sample_", 1:20)

  metadata <- data.frame(
    sample_id = paste0("Sample_", 1:20),
    disease_status = sample(c(0, 1), 20, replace = TRUE),
    age = rnorm(20, mean = 45, sd = 10),
    sex = sample(c("M", "F"), 20, replace = TRUE),
    stringsAsFactors = FALSE
  )

  # Run analysis
  result <- olink_regression_analysis(
    count_matrix = count_matrix,
    metadata = metadata,
    outcome = "disease_status",
    outcome_type = "binary",
    covariates = NULL,
    sample_id_col = "sample_id",
    verbose = FALSE
  )

  # Check structure
  expect_type(result, "list")
  expect_named(result, c("univariate_results", "multivariate_results",
                         "summary_stats", "model_type"))

  # Check model type
  expect_equal(result$model_type, "logistic")

  # Check univariate results columns
  expect_true(all(c("protein", "OR", "OR_lower", "OR_upper", "log_OR",
                    "SE", "p_value", "n_samples", "p_adjusted") %in%
                    colnames(result$univariate_results)))

  # Check that OR values are positive
  expect_true(all(result$univariate_results$OR > 0))

  # Check that p-values are between 0 and 1
  expect_true(all(result$univariate_results$p_value >= 0 &
                    result$univariate_results$p_value <= 1))

  # Check number of proteins analyzed
  expect_true(nrow(result$univariate_results) <= nrow(count_matrix))

  # Check that multivariate_results is NULL when no covariates
  expect_null(result$multivariate_results)
})

test_that("olink_regression_analysis works with continuous outcome (linear)", {
  set.seed(456)

  # Create test data with continuous outcome
  count_matrix <- matrix(rnorm(150, mean = 80, sd = 15),
                         nrow = 10, ncol = 15)
  rownames(count_matrix) <- paste0("Protein_", 1:10)
  colnames(count_matrix) <- paste0("Sample_", 1:15)

  metadata <- data.frame(
    sample_id = paste0("Sample_", 1:15),
    severity_score = rnorm(15, mean = 50, sd = 20),
    age = rnorm(15, mean = 50, sd = 12),
    stringsAsFactors = FALSE
  )

  # Run analysis
  result <- olink_regression_analysis(
    count_matrix = count_matrix,
    metadata = metadata,
    outcome = "severity_score",
    outcome_type = "continuous",
    covariates = NULL,
    sample_id_col = "sample_id",
    verbose = FALSE
  )

  # Check structure
  expect_type(result, "list")
  expect_equal(result$model_type, "linear")

  # Check univariate results columns for continuous outcome
  expect_true(all(c("protein", "beta", "beta_lower", "beta_upper",
                    "SE", "p_value", "n_samples", "p_adjusted") %in%
                    colnames(result$univariate_results)))

  # Check that confidence intervals make sense
  expect_true(all(result$univariate_results$beta_lower <=
                    result$univariate_results$beta))
  expect_true(all(result$univariate_results$beta <=
                    result$univariate_results$beta_upper))
})

test_that("olink_regression_analysis performs multivariate analysis", {
  set.seed(789)

  count_matrix <- matrix(rnorm(200, mean = 100, sd = 25),
                         nrow = 10, ncol = 20)
  rownames(count_matrix) <- paste0("Protein_", 1:10)
  colnames(count_matrix) <- paste0("Sample_", 1:20)

  metadata <- data.frame(
    sample_id = paste0("Sample_", 1:20),
    outcome = sample(c(0, 1), 20, replace = TRUE),
    age = rnorm(20, mean = 45, sd = 10),
    sex = sample(c("M", "F"), 20, replace = TRUE),
    bmi = rnorm(20, mean = 25, sd = 4),
    stringsAsFactors = FALSE
  )

  # Run with covariates
  result <- olink_regression_analysis(
    count_matrix = count_matrix,
    metadata = metadata,
    outcome = "outcome",
    outcome_type = "binary",
    covariates = c("age", "sex", "bmi"),
    sample_id_col = "sample_id",
    verbose = FALSE
  )

  # Check that multivariate results exist
  expect_false(is.null(result$multivariate_results))
  expect_s3_class(result$multivariate_results, "data.frame")

  # Check multivariate results have correct columns
  expect_true(all(c("protein", "OR", "OR_lower", "OR_upper", "log_OR",
                    "SE", "p_value", "n_samples", "p_adjusted") %in%
                    colnames(result$multivariate_results)))

  # Check summary stats include multivariate info
  expect_true("n_proteins_analyzed_multivariate" %in% names(result$summary_stats))
  expect_true("n_significant_multivariate_raw" %in% names(result$summary_stats))
  expect_true("covariates_used" %in% names(result$summary_stats))
  expect_equal(result$summary_stats$covariates_used, c("age", "sex", "bmi"))
})

test_that("olink_regression_analysis auto-detects outcome type", {
  set.seed(111)

  count_matrix <- matrix(rnorm(100, mean = 50, sd = 10),
                         nrow = 5, ncol = 20)
  rownames(count_matrix) <- paste0("Protein_", 1:5)
  colnames(count_matrix) <- paste0("Sample_", 1:20)

  # Binary outcome (auto-detect)
  metadata_binary <- data.frame(
    sample_id = paste0("Sample_", 1:20),
    outcome = sample(c(0, 1), 20, replace = TRUE),
    stringsAsFactors = FALSE
  )

  result_binary <- olink_regression_analysis(
    count_matrix = count_matrix,
    metadata = metadata_binary,
    outcome = "outcome",
    outcome_type = NULL,  # Auto-detect
    sample_id_col = "sample_id",
    verbose = FALSE
  )

  expect_equal(result_binary$model_type, "logistic")

  # Continuous outcome (auto-detect)
  metadata_cont <- data.frame(
    sample_id = paste0("Sample_", 1:20),
    outcome = rnorm(20, mean = 100, sd = 30),
    stringsAsFactors = FALSE
  )

  result_cont <- olink_regression_analysis(
    count_matrix = count_matrix,
    metadata = metadata_cont,
    outcome = "outcome",
    outcome_type = NULL,  # Auto-detect
    sample_id_col = "sample_id",
    verbose = FALSE
  )

  expect_equal(result_cont$model_type, "linear")
})

test_that("olink_regression_analysis handles log transformation", {
  set.seed(222)

  count_matrix <- matrix(abs(rnorm(100, mean = 100, sd = 30)),
                         nrow = 5, ncol = 20)
  rownames(count_matrix) <- paste0("Protein_", 1:5)
  colnames(count_matrix) <- paste0("Sample_", 1:20)

  metadata <- data.frame(
    sample_id = paste0("Sample_", 1:20),
    outcome = sample(c(0, 1), 20, replace = TRUE),
    stringsAsFactors = FALSE
  )

  # With log transformation
  result_log <- olink_regression_analysis(
    count_matrix = count_matrix,
    metadata = metadata,
    outcome = "outcome",
    outcome_type = "binary",
    log_transform = TRUE,
    sample_id_col = "sample_id",
    verbose = FALSE
  )

  # Without log transformation
  result_no_log <- olink_regression_analysis(
    count_matrix = count_matrix,
    metadata = metadata,
    outcome = "outcome",
    outcome_type = "binary",
    log_transform = FALSE,
    sample_id_col = "sample_id",
    verbose = FALSE
  )

  expect_type(result_log, "list")
  expect_type(result_no_log, "list")

  # Results should differ
  expect_false(identical(result_log$univariate_results$OR,
                        result_no_log$univariate_results$OR))
})

test_that("olink_regression_analysis applies multiple testing correction", {
  set.seed(333)

  count_matrix <- matrix(rnorm(500, mean = 100, sd = 20),
                         nrow = 25, ncol = 20)
  rownames(count_matrix) <- paste0("Protein_", 1:25)
  colnames(count_matrix) <- paste0("Sample_", 1:20)

  metadata <- data.frame(
    sample_id = paste0("Sample_", 1:20),
    outcome = sample(c(0, 1), 20, replace = TRUE),
    stringsAsFactors = FALSE
  )

  # Test different correction methods
  result_bh <- olink_regression_analysis(
    count_matrix = count_matrix,
    metadata = metadata,
    outcome = "outcome",
    outcome_type = "binary",
    p_adjust_method = "BH",
    sample_id_col = "sample_id",
    verbose = FALSE
  )

  result_bonf <- olink_regression_analysis(
    count_matrix = count_matrix,
    metadata = metadata,
    outcome = "outcome",
    outcome_type = "binary",
    p_adjust_method = "bonferroni",
    sample_id_col = "sample_id",
    verbose = FALSE
  )

  # Check p_adjusted column exists
  expect_true("p_adjusted" %in% colnames(result_bh$univariate_results))
  expect_true("p_adjusted" %in% colnames(result_bonf$univariate_results))

  # Bonferroni should be more conservative
  expect_true(mean(result_bonf$univariate_results$p_adjusted) >=
              mean(result_bh$univariate_results$p_adjusted))

  # Check correction method in summary
  expect_equal(result_bh$summary_stats$p_adjust_method, "BH")
  expect_equal(result_bonf$summary_stats$p_adjust_method, "bonferroni")
})

test_that("olink_regression_analysis handles missing data", {
  set.seed(444)

  count_matrix <- matrix(rnorm(100, mean = 80, sd = 15),
                         nrow = 5, ncol = 20)
  rownames(count_matrix) <- paste0("Protein_", 1:5)
  colnames(count_matrix) <- paste0("Sample_", 1:20)

  # Add some NA values
  count_matrix[1, 1:3] <- NA
  count_matrix[2, 5:7] <- NA

  metadata <- data.frame(
    sample_id = paste0("Sample_", 1:20),
    outcome = sample(c(0, 1), 20, replace = TRUE),
    age = rnorm(20, mean = 50, sd = 12),
    stringsAsFactors = FALSE
  )

  # Add NA in outcome
  metadata$outcome[c(2, 4)] <- NA

  result <- olink_regression_analysis(
    count_matrix = count_matrix,
    metadata = metadata,
    outcome = "outcome",
    outcome_type = "binary",
    covariates = "age",
    sample_id_col = "sample_id",
    min_samples = 5,
    verbose = FALSE
  )

  # Should still produce results
  expect_type(result, "list")
  expect_true(nrow(result$univariate_results) > 0)

  # Sample counts should reflect complete cases
  expect_true(all(result$univariate_results$n_samples <= 20))
})

test_that("olink_regression_analysis handles partial sample overlap", {
  set.seed(555)

  count_matrix <- matrix(rnorm(60, mean = 100, sd = 20),
                         nrow = 6, ncol = 10)
  rownames(count_matrix) <- paste0("Protein_", 1:6)
  colnames(count_matrix) <- paste0("Sample_", 1:10)

  # Metadata with only partial overlap (samples 5-15)
  metadata <- data.frame(
    sample_id = paste0("Sample_", 5:15),
    outcome = rnorm(11, mean = 50, sd = 20),
    stringsAsFactors = FALSE
  )

  result <- olink_regression_analysis(
    count_matrix = count_matrix,
    metadata = metadata,
    outcome = "outcome",
    outcome_type = "continuous",
    sample_id_col = "sample_id",
    min_samples = 3,
    verbose = FALSE
  )

  # Should only use overlapping samples (5-10)
  expect_equal(result$summary_stats$n_samples, 6)
})

test_that("olink_regression_analysis validates inputs correctly", {
  count_matrix <- matrix(1:20, nrow = 4, ncol = 5)
  rownames(count_matrix) <- paste0("Protein_", 1:4)
  colnames(count_matrix) <- paste0("Sample_", 1:5)

  metadata <- data.frame(
    sample_id = paste0("Sample_", 1:5),
    outcome = c(0, 1, 0, 1, 0)
  )

  # Should error on invalid count_matrix type
  expect_error(
    olink_regression_analysis(
      count_matrix = "not a matrix",
      metadata = metadata,
      outcome = "outcome"
    ),
    "count_matrix must be a matrix or data.frame"
  )

  # Should error on invalid metadata type
  expect_error(
    olink_regression_analysis(
      count_matrix = count_matrix,
      metadata = "not a dataframe",
      outcome = "outcome"
    ),
    "metadata must be a data.frame"
  )

  # Should error on non-existent outcome
  expect_error(
    olink_regression_analysis(
      count_matrix = count_matrix,
      metadata = metadata,
      outcome = "nonexistent_variable",
      sample_id_col = "sample_id"
    ),
    "Outcome variable nonexistent_variable not found"
  )

  # Should error on non-existent covariates
  expect_error(
    olink_regression_analysis(
      count_matrix = count_matrix,
      metadata = metadata,
      outcome = "outcome",
      covariates = c("age", "sex"),
      sample_id_col = "sample_id"
    ),
    "Covariates not found in metadata"
  )
})

test_that("olink_regression_analysis handles insufficient samples", {
  count_matrix <- matrix(1:12, nrow = 3, ncol = 4)
  rownames(count_matrix) <- paste0("Protein_", 1:3)
  colnames(count_matrix) <- paste0("Sample_", 1:4)

  metadata <- data.frame(
    sample_id = paste0("Sample_", 1:4),
    outcome = c(0, 1, 0, 1)
  )

  # Should error when samples < min_samples
  expect_error(
    olink_regression_analysis(
      count_matrix = count_matrix,
      metadata = metadata,
      outcome = "outcome",
      sample_id_col = "sample_id",
      min_samples = 10,
      verbose = FALSE
    ),
    "Insufficient samples"
  )
})

test_that("olink_regression_analysis results are sorted by p-value", {
  set.seed(666)

  count_matrix <- matrix(rnorm(200, mean = 100, sd = 20),
                         nrow = 10, ncol = 20)
  rownames(count_matrix) <- paste0("Protein_", 1:10)
  colnames(count_matrix) <- paste0("Sample_", 1:20)

  metadata <- data.frame(
    sample_id = paste0("Sample_", 1:20),
    outcome = sample(c(0, 1), 20, replace = TRUE),
    stringsAsFactors = FALSE
  )

  result <- olink_regression_analysis(
    count_matrix = count_matrix,
    metadata = metadata,
    outcome = "outcome",
    outcome_type = "binary",
    sample_id_col = "sample_id",
    verbose = FALSE
  )

  # Check that results are sorted by p-value (ascending)
  expect_true(all(diff(result$univariate_results$p_value) >= 0))
})

test_that("olink_regression_analysis uses rownames when sample_id_col is NULL", {
  set.seed(777)

  count_matrix <- matrix(rnorm(60, mean = 100, sd = 20),
                         nrow = 6, ncol = 10)
  rownames(count_matrix) <- paste0("Protein_", 1:6)
  colnames(count_matrix) <- paste0("Sample_", 1:10)

  metadata <- data.frame(
    outcome = sample(c(0, 1), 10, replace = TRUE),
    age = rnorm(10, mean = 45, sd = 10)
  )
  rownames(metadata) <- paste0("Sample_", 1:10)

  result <- olink_regression_analysis(
    count_matrix = count_matrix,
    metadata = metadata,
    outcome = "outcome",
    outcome_type = "binary",
    sample_id_col = NULL,
    verbose = FALSE
  )

  expect_type(result, "list")
  expect_equal(result$summary_stats$n_samples, 10)
})

test_that("summary_stats contains expected information", {
  set.seed(888)

  count_matrix <- matrix(rnorm(150, mean = 100, sd = 20),
                         nrow = 10, ncol = 15)
  rownames(count_matrix) <- paste0("Protein_", 1:10)
  colnames(count_matrix) <- paste0("Sample_", 1:15)

  metadata <- data.frame(
    sample_id = paste0("Sample_", 1:15),
    outcome = sample(c(0, 1), 15, replace = TRUE),
    age = rnorm(15, mean = 45, sd = 10),
    sex = sample(c("M", "F"), 15, replace = TRUE)
  )

  result <- olink_regression_analysis(
    count_matrix = count_matrix,
    metadata = metadata,
    outcome = "outcome",
    outcome_type = "binary",
    covariates = c("age", "sex"),
    sample_id_col = "sample_id",
    verbose = FALSE
  )

  # Check summary stats structure
  expect_true("n_proteins_total" %in% names(result$summary_stats))
  expect_true("n_proteins_analyzed_univariate" %in% names(result$summary_stats))
  expect_true("n_samples" %in% names(result$summary_stats))
  expect_true("outcome_variable" %in% names(result$summary_stats))
  expect_true("model_type" %in% names(result$summary_stats))
  expect_true("p_adjust_method" %in% names(result$summary_stats))
  expect_true("n_significant_univariate_raw" %in% names(result$summary_stats))
  expect_true("n_significant_univariate_adj" %in% names(result$summary_stats))

  # Check values
  expect_equal(result$summary_stats$n_proteins_total, 10)
  expect_equal(result$summary_stats$n_samples, 15)
  expect_equal(result$summary_stats$outcome_variable, "outcome")
  expect_equal(result$summary_stats$model_type, "logistic")
})

test_that("confidence intervals respect conf_level parameter", {
  set.seed(999)

  count_matrix <- matrix(rnorm(60, mean = 100, sd = 20),
                         nrow = 6, ncol = 10)
  rownames(count_matrix) <- paste0("Protein_", 1:6)
  colnames(count_matrix) <- paste0("Sample_", 1:10)

  metadata <- data.frame(
    sample_id = paste0("Sample_", 1:10),
    outcome = rnorm(10, mean = 50, sd = 15)
  )

  # 95% CI
  result_95 <- olink_regression_analysis(
    count_matrix = count_matrix,
    metadata = metadata,
    outcome = "outcome",
    outcome_type = "continuous",
    conf_level = 0.95,
    sample_id_col = "sample_id",
    verbose = FALSE
  )

  # 99% CI
  result_99 <- olink_regression_analysis(
    count_matrix = count_matrix,
    metadata = metadata,
    outcome = "outcome",
    outcome_type = "continuous",
    conf_level = 0.99,
    sample_id_col = "sample_id",
    verbose = FALSE
  )

  # 99% CI should be wider than 95% CI
  width_95 <- result_95$univariate_results$beta_upper -
              result_95$univariate_results$beta_lower
  width_99 <- result_99$univariate_results$beta_upper -
              result_99$univariate_results$beta_lower

  expect_true(all(width_99 >= width_95))
})
