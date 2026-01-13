# Test file for Olink covariate adjustment functions

test_that("adjust_olink_covariates works with basic inputs", {
  # Create test data
  set.seed(123)
  count_matrix <- matrix(rnorm(200, mean = 100, sd = 20),
                         nrow = 10, ncol = 20)
  rownames(count_matrix) <- paste0("Protein_", 1:10)
  colnames(count_matrix) <- paste0("Sample_", 1:20)

  metadata <- data.frame(
    sample_id = paste0("Sample_", 1:20),
    age = rnorm(20, mean = 45, sd = 10),
    sex = sample(c("M", "F"), 20, replace = TRUE),
    batch = sample(c("Batch1", "Batch2"), 20, replace = TRUE),
    stringsAsFactors = FALSE
  )

  # Run adjustment
  result <- adjust_olink_covariates(
    count_matrix = count_matrix,
    metadata = metadata,
    covariates = c("age", "sex"),
    sample_id_col = "sample_id",
    method = "residuals",
    verbose = FALSE
  )

  # Check structure
  expect_type(result, "list")
  expect_named(result, c("adjusted_matrix", "model_info", "covariates_used",
                         "samples_used", "proteins_failed"))

  # Check dimensions
  expect_equal(dim(result$adjusted_matrix), dim(count_matrix))
  expect_equal(rownames(result$adjusted_matrix), rownames(count_matrix))
  expect_equal(colnames(result$adjusted_matrix), colnames(count_matrix))

  # Check covariates used
  expect_equal(result$covariates_used, c("age", "sex"))

  # Check samples used
  expect_equal(length(result$samples_used), 20)
})

test_that("adjust_olink_covariates handles adjusted method", {
  set.seed(456)
  count_matrix <- matrix(rnorm(100, mean = 50, sd = 10),
                         nrow = 5, ncol = 20)
  rownames(count_matrix) <- paste0("Protein_", 1:5)
  colnames(count_matrix) <- paste0("Sample_", 1:20)

  metadata <- data.frame(
    sample_id = paste0("Sample_", 1:20),
    age = rnorm(20, mean = 40, sd = 8),
    stringsAsFactors = FALSE
  )

  result <- adjust_olink_covariates(
    count_matrix = count_matrix,
    metadata = metadata,
    covariates = "age",
    sample_id_col = "sample_id",
    method = "adjusted",
    verbose = FALSE
  )

  expect_type(result, "list")
  expect_equal(nrow(result$adjusted_matrix), 5)
  expect_equal(ncol(result$adjusted_matrix), 20)
})

test_that("adjust_olink_covariates handles log transformation", {
  set.seed(789)
  count_matrix <- matrix(abs(rnorm(100, mean = 100, sd = 30)),
                         nrow = 5, ncol = 20)
  rownames(count_matrix) <- paste0("Protein_", 1:5)
  colnames(count_matrix) <- paste0("Sample_", 1:20)

  metadata <- data.frame(
    sample_id = paste0("Sample_", 1:20),
    batch = sample(c("A", "B", "C"), 20, replace = TRUE),
    stringsAsFactors = FALSE
  )

  # With log transformation, return log
  result_log <- adjust_olink_covariates(
    count_matrix = count_matrix,
    metadata = metadata,
    covariates = "batch",
    sample_id_col = "sample_id",
    log_transform = TRUE,
    return_log = TRUE,
    verbose = FALSE
  )

  # With log transformation, return original scale
  result_orig <- adjust_olink_covariates(
    count_matrix = count_matrix,
    metadata = metadata,
    covariates = "batch",
    sample_id_col = "sample_id",
    log_transform = TRUE,
    return_log = FALSE,
    verbose = FALSE
  )

  expect_type(result_log, "list")
  expect_type(result_orig, "list")

  # Log values should be smaller in magnitude
  expect_true(mean(abs(result_log$adjusted_matrix), na.rm = TRUE) <
              mean(abs(result_orig$adjusted_matrix), na.rm = TRUE))
})

test_that("adjust_olink_covariates handles missing data", {
  set.seed(111)
  count_matrix <- matrix(rnorm(100, mean = 80, sd = 15),
                         nrow = 5, ncol = 20)
  rownames(count_matrix) <- paste0("Protein_", 1:5)
  colnames(count_matrix) <- paste0("Sample_", 1:20)

  # Add some NA values
  count_matrix[1, 1:3] <- NA
  count_matrix[2, 5:7] <- NA

  metadata <- data.frame(
    sample_id = paste0("Sample_", 1:20),
    age = rnorm(20, mean = 50, sd = 12),
    stringsAsFactors = FALSE
  )

  # Add some NA in metadata
  metadata$age[c(2, 4)] <- NA

  result <- adjust_olink_covariates(
    count_matrix = count_matrix,
    metadata = metadata,
    covariates = "age",
    sample_id_col = "sample_id",
    verbose = FALSE
  )

  expect_type(result, "list")
  # Should still produce results for valid data
  expect_true(sum(!is.na(result$adjusted_matrix)) > 0)
})

test_that("adjust_olink_covariates handles multiple covariates", {
  set.seed(222)
  count_matrix <- matrix(rnorm(150, mean = 100, sd = 25),
                         nrow = 10, ncol = 15)
  rownames(count_matrix) <- paste0("Protein_", 1:10)
  colnames(count_matrix) <- paste0("Sample_", 1:15)

  metadata <- data.frame(
    sample_id = paste0("Sample_", 1:15),
    age = rnorm(15, mean = 45, sd = 10),
    sex = sample(c("M", "F"), 15, replace = TRUE),
    batch = sample(c("B1", "B2", "B3"), 15, replace = TRUE),
    bmi = rnorm(15, mean = 25, sd = 4),
    stringsAsFactors = FALSE
  )

  result <- adjust_olink_covariates(
    count_matrix = count_matrix,
    metadata = metadata,
    covariates = c("age", "sex", "batch", "bmi"),
    sample_id_col = "sample_id",
    verbose = FALSE
  )

  expect_type(result, "list")
  expect_equal(length(result$covariates_used), 4)
  expect_true(all(c("age", "sex", "batch", "bmi") %in% result$covariates_used))

  # Check model info contains R-squared values
  expect_true(length(result$model_info) > 0)
  expect_true("r_squared" %in% names(result$model_info[[1]]))
})

test_that("adjust_olink_covariates validates inputs correctly", {
  count_matrix <- matrix(1:20, nrow = 4, ncol = 5)
  metadata <- data.frame(id = 1:5, age = rnorm(5))

  # Should error on invalid count_matrix type
  expect_error(
    adjust_olink_covariates(
      count_matrix = "not a matrix",
      metadata = metadata,
      covariates = "age"
    ),
    "count_matrix must be a matrix or data.frame"
  )

  # Should error on invalid metadata type
  expect_error(
    adjust_olink_covariates(
      count_matrix = count_matrix,
      metadata = "not a dataframe",
      covariates = "age"
    ),
    "metadata must be a data.frame"
  )

  # Should error on empty covariates
  expect_error(
    adjust_olink_covariates(
      count_matrix = count_matrix,
      metadata = metadata,
      covariates = character(0)
    ),
    "At least one covariate must be specified"
  )
})

test_that("adjust_olink_covariates handles mismatched samples", {
  count_matrix <- matrix(rnorm(40), nrow = 4, ncol = 10)
  rownames(count_matrix) <- paste0("Protein_", 1:4)
  colnames(count_matrix) <- paste0("Sample_", 1:10)

  # Metadata with different sample names
  metadata <- data.frame(
    sample_id = paste0("DifferentSample_", 1:10),
    age = rnorm(10),
    stringsAsFactors = FALSE
  )

  # Should error on no common samples
  expect_error(
    adjust_olink_covariates(
      count_matrix = count_matrix,
      metadata = metadata,
      covariates = "age",
      sample_id_col = "sample_id",
      verbose = FALSE
    ),
    "No common samples found"
  )
})

test_that("adjust_olink_covariates handles partial sample overlap", {
  set.seed(333)
  count_matrix <- matrix(rnorm(60), nrow = 6, ncol = 10)
  rownames(count_matrix) <- paste0("Protein_", 1:6)
  colnames(count_matrix) <- paste0("Sample_", 1:10)

  # Metadata with only partial overlap
  metadata <- data.frame(
    sample_id = paste0("Sample_", 5:15),  # Only 5-10 overlap
    age = rnorm(11),
    stringsAsFactors = FALSE
  )

  result <- adjust_olink_covariates(
    count_matrix = count_matrix,
    metadata = metadata,
    covariates = "age",
    sample_id_col = "sample_id",
    verbose = FALSE
  )

  # Should only use overlapping samples
  expect_equal(length(result$samples_used), 6)
  expect_true(all(result$samples_used %in% paste0("Sample_", 5:10)))
})

test_that("adjust_olink_covariates uses rownames when sample_id_col is NULL", {
  set.seed(444)
  count_matrix <- matrix(rnorm(40), nrow = 4, ncol = 10)
  rownames(count_matrix) <- paste0("Protein_", 1:4)
  colnames(count_matrix) <- paste0("Sample_", 1:10)

  metadata <- data.frame(
    age = rnorm(10),
    sex = sample(c("M", "F"), 10, replace = TRUE),
    stringsAsFactors = FALSE
  )
  rownames(metadata) <- paste0("Sample_", 1:10)

  result <- adjust_olink_covariates(
    count_matrix = count_matrix,
    metadata = metadata,
    covariates = "age",
    sample_id_col = NULL,
    verbose = FALSE
  )

  expect_type(result, "list")
  expect_equal(length(result$samples_used), 10)
})

test_that("adjust_olink_covariates centers numeric covariates", {
  set.seed(555)
  count_matrix <- matrix(rnorm(50), nrow = 5, ncol = 10)
  rownames(count_matrix) <- paste0("Protein_", 1:5)
  colnames(count_matrix) <- paste0("Sample_", 1:10)

  metadata <- data.frame(
    sample_id = paste0("Sample_", 1:10),
    age = seq(20, 65, length.out = 10),
    stringsAsFactors = FALSE
  )

  result_centered <- adjust_olink_covariates(
    count_matrix = count_matrix,
    metadata = metadata,
    covariates = "age",
    sample_id_col = "sample_id",
    center_covariates = TRUE,
    verbose = FALSE
  )

  result_not_centered <- adjust_olink_covariates(
    count_matrix = count_matrix,
    metadata = metadata,
    covariates = "age",
    sample_id_col = "sample_id",
    center_covariates = FALSE,
    verbose = FALSE
  )

  expect_type(result_centered, "list")
  expect_type(result_not_centered, "list")

  # Both should produce valid results
  expect_true(all(!is.na(result_centered$adjusted_matrix)))
  expect_true(all(!is.na(result_not_centered$adjusted_matrix)))
})

test_that("model_info contains expected components", {
  set.seed(666)
  count_matrix <- matrix(rnorm(60), nrow = 3, ncol = 20)
  rownames(count_matrix) <- paste0("Protein_", 1:3)
  colnames(count_matrix) <- paste0("Sample_", 1:20)

  metadata <- data.frame(
    sample_id = paste0("Sample_", 1:20),
    age = rnorm(20, 45, 10),
    sex = sample(c("M", "F"), 20, replace = TRUE),
    stringsAsFactors = FALSE
  )

  result <- adjust_olink_covariates(
    count_matrix = count_matrix,
    metadata = metadata,
    covariates = c("age", "sex"),
    sample_id_col = "sample_id",
    verbose = FALSE
  )

  # Check model info structure
  expect_true(length(result$model_info) > 0)

  for (protein in names(result$model_info)) {
    info <- result$model_info[[protein]]
    expect_true("r_squared" %in% names(info))
    expect_true("adj_r_squared" %in% names(info))
    expect_true("n_samples" %in% names(info))
    expect_true("coefficients" %in% names(info))

    expect_true(is.numeric(info$r_squared))
    expect_true(info$r_squared >= 0 && info$r_squared <= 1)
    expect_true(is.numeric(info$n_samples))
    expect_true(info$n_samples > 0)
  }
})
