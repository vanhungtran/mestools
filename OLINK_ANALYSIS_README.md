# Olink Proteomics Data Analysis in mestools

Complete toolkit for analyzing Olink protein expression data with covariate adjustment and comprehensive regression analyses.

## 📋 Overview

This package provides three main functions for Olink data analysis:

1. **`adjust_olink_covariates()`** - Remove batch effects and confounders
2. **`olink_regression_analysis()`** - Univariate & multivariate association testing
3. **`plot_forest_results()`** - Visualize odds ratios and beta coefficients

## 🚀 Quick Start

### Installation

```r
# Install from GitHub
devtools::install_github("vanhungtran/mestools")

# Load the package
library(mestools)
```

### Basic Usage

```r
# 1. Adjust for batch effects and covariates
adjusted_data <- adjust_olink_covariates(
  count_matrix = your_protein_matrix,
  metadata = your_metadata,
  covariates = c("age", "sex", "batch"),
  sample_id_col = "sample_id",
  method = "adjusted"
)

# 2. Test associations with disease
results <- olink_regression_analysis(
  count_matrix = adjusted_data$adjusted_matrix,
  metadata = your_metadata,
  outcome = "disease_status",
  outcome_type = "binary",
  covariates = c("age", "sex", "bmi"),
  sample_id_col = "sample_id"
)

# 3. Visualize results
plot_forest_results(
  results = results,
  analysis_type = "multivariate",
  top_n = 20
)
```

## 📊 Features

### Covariate Adjustment

- **Remove confounding effects** from protein expression data
- **Batch effect correction** for multi-batch experiments
- **Flexible adjustment methods**:
  - Residuals method (for downstream analysis)
  - Adjusted values method (preserves scale)
- **Log transformation** support with back-transformation
- **Automatic centering** of numeric covariates
- **Model diagnostics** (R², coefficients, sample counts)

### Regression Analysis

#### Binary Outcomes (Logistic Regression)
- **Odds ratios with 95% CI**
- Case-control studies
- Disease association studies
- Treatment response prediction

#### Continuous Outcomes (Linear Regression)
- **Beta coefficients with 95% CI**
- Severity scores
- Biomarker levels
- Quantitative traits

#### Both Include:
- **Univariate analysis** (one protein at a time)
- **Multivariate analysis** (adjusted for covariates)
- **Multiple testing correction** (BH, Bonferroni, Holm, etc.)
- **Missing data handling**
- **Comprehensive output** with p-values, effect sizes, sample counts

### Visualization

- **Forest plots** for odds ratios or beta coefficients
- **Volcano plots** for binary outcomes
- **Manhattan plots** for genome-wide style views
- **Before/after plots** for covariate adjustment
- **Customizable aesthetics**

## 📖 Examples with Figures

See the **[olink_examples/](olink_examples/)** directory for:
- 📈 **10+ example visualizations**
- 📄 **Detailed README with interpretations**
- 📊 **Sample result tables**
- 💻 **Reproducible code examples**

### Example Visualizations

| Analysis Type | Visualization |
|--------------|---------------|
| Batch Effect Removal | ![](olink_examples/batch_effect_removal.png) |
| Forest Plot (Binary) | ![](olink_examples/forest_plot_univariate_binary.png) |
| Volcano Plot | ![](olink_examples/volcano_plot_binary.png) |
| Beta Coefficients | ![](olink_examples/forest_plot_continuous.png) |

*Full-size images available in [olink_examples/README.md](olink_examples/README.md)*

## 📚 Documentation

### Function Reference

#### `adjust_olink_covariates()`

```r
adjust_olink_covariates(
  count_matrix,           # Protein x Sample matrix
  metadata,               # Sample metadata data.frame
  covariates,             # Vector of covariate names
  sample_id_col = NULL,   # Sample ID column name
  method = "residuals",   # "residuals" or "adjusted"
  log_transform = FALSE,  # Log-transform before adjustment
  return_log = TRUE,      # Return log scale or back-transform
  center_covariates = TRUE,  # Center numeric covariates
  verbose = TRUE          # Print progress
)
```

**Returns:**
- `adjusted_matrix` - Covariate-adjusted expression matrix
- `model_info` - Model diagnostics for each protein
- `covariates_used` - Covariates that were adjusted
- `samples_used` - Samples included in analysis
- `proteins_failed` - Proteins where adjustment failed

#### `olink_regression_analysis()`

```r
olink_regression_analysis(
  count_matrix,           # Protein x Sample matrix
  metadata,               # Sample metadata data.frame
  outcome,                # Outcome variable name
  outcome_type = NULL,    # "binary" or "continuous" (auto-detect)
  covariates = NULL,      # Covariates for multivariate model
  sample_id_col = NULL,   # Sample ID column name
  log_transform = FALSE,  # Log-transform proteins
  conf_level = 0.95,      # Confidence level for CI
  p_adjust_method = "BH", # Multiple testing correction
  min_samples = 10,       # Minimum samples per protein
  verbose = TRUE          # Print progress
)
```

**Returns:**
- `univariate_results` - Univariate associations (data.frame)
- `multivariate_results` - Multivariate associations (data.frame)
- `summary_stats` - Analysis summary (list)
- `model_type` - "logistic" or "linear"

#### Result Data Frames

**Binary Outcome (Logistic):**
| Column | Description |
|--------|-------------|
| `protein` | Protein identifier |
| `OR` | Odds ratio |
| `OR_lower` | Lower 95% CI |
| `OR_upper` | Upper 95% CI |
| `log_OR` | Log odds ratio (coefficient) |
| `SE` | Standard error |
| `p_value` | Raw p-value |
| `p_adjusted` | FDR-adjusted p-value |
| `n_samples` | Number of samples used |

**Continuous Outcome (Linear):**
| Column | Description |
|--------|-------------|
| `protein` | Protein identifier |
| `beta` | Beta coefficient |
| `beta_lower` | Lower 95% CI |
| `beta_upper` | Upper 95% CI |
| `SE` | Standard error |
| `p_value` | Raw p-value |
| `p_adjusted` | FDR-adjusted p-value |
| `n_samples` | Number of samples used |

#### `plot_forest_results()`

```r
plot_forest_results(
  results,                      # Results from olink_regression_analysis()
  analysis_type = "univariate", # "univariate" or "multivariate"
  top_n = 20,                   # Number of top proteins to show
  sort_by = "p_value",          # "p_value" or "effect"
  show_pval = TRUE,             # Show p-values on plot
  title = NULL                  # Custom plot title
)
```

**Returns:** ggplot2 object

## 🔬 Use Cases

### 1. Batch Effect Correction

```r
# Multi-center study with batch effects
corrected <- adjust_olink_covariates(
  count_matrix = raw_data,
  metadata = sample_info,
  covariates = c("batch", "plate", "run_date"),
  method = "adjusted"
)

# Visualize correction
plot_adjustment_results(
  original_matrix = raw_data,
  adjusted_matrix = corrected$adjusted_matrix,
  metadata = sample_info,
  covariate = "batch"
)
```

### 2. Case-Control Association Study

```r
# Identify disease-associated proteins
results <- olink_regression_analysis(
  count_matrix = protein_data,
  metadata = clinical_data,
  outcome = "case_control",  # 0 = control, 1 = case
  outcome_type = "binary",
  covariates = c("age", "sex", "bmi", "smoking"),
  p_adjust_method = "BH"
)

# Get significant proteins (FDR < 0.05)
sig_proteins <- results$multivariate_results[
  results$multivariate_results$p_adjusted < 0.05,
]

# Create visualization
plot_forest_results(results, "multivariate", top_n = 30)
```

### 3. Biomarker Discovery for Continuous Trait

```r
# Find proteins associated with disease severity
results <- olink_regression_analysis(
  count_matrix = protein_data,
  metadata = clinical_data,
  outcome = "severity_score",  # Continuous 0-100
  outcome_type = "continuous",
  covariates = c("age", "sex", "disease_duration"),
  log_transform = TRUE
)

# Top positive associations (higher protein → higher severity)
top_positive <- results$univariate_results[
  results$univariate_results$beta > 0,
][1:10, ]

# Top negative associations (higher protein → lower severity)
top_negative <- results$univariate_results[
  results$univariate_results$beta < 0,
][1:10, ]
```

### 4. Multi-Stage Analysis Workflow

```r
# Stage 1: Correct for technical factors
adjusted <- adjust_olink_covariates(
  count_matrix = raw_data,
  metadata = metadata,
  covariates = c("batch", "plate"),
  method = "adjusted"
)

# Stage 2: Test biological associations
results_bio <- olink_regression_analysis(
  count_matrix = adjusted$adjusted_matrix,
  metadata = metadata,
  outcome = "disease_status",
  outcome_type = "binary",
  covariates = c("age", "sex", "ethnicity"),
  log_transform = FALSE  # Already adjusted
)

# Stage 3: Visualize and export
plot_forest_results(results_bio, "multivariate", top_n = 25)

write.csv(
  results_bio$multivariate_results,
  "protein_associations.csv",
  row.names = FALSE
)
```

## 🧪 Testing

The package includes comprehensive test suites:

```r
# Run all tests
devtools::test()

# Run specific test files
testthat::test_file("tests/testthat/test-olink-covariate-adjustment.R")
testthat::test_file("tests/testthat/test-olink-regression-analysis.R")
```

**Test coverage:**
- ✅ 15 tests for covariate adjustment
- ✅ 20 tests for regression analysis
- ✅ Input validation
- ✅ Missing data handling
- ✅ Edge cases
- ✅ Output structure verification

## 📝 Data Format Requirements

### Count Matrix

```r
# Protein × Sample matrix
#               Sample_1  Sample_2  Sample_3  ...
# Protein_A        150.2     142.8     165.3
# Protein_B         85.1      91.2      78.9
# Protein_C        210.5     198.7     225.1
# ...
```

- Rows = Proteins
- Columns = Samples
- Values = Expression levels (NPX, counts, or intensity)
- Row names = Protein IDs
- Column names = Sample IDs

### Metadata

```r
# Sample metadata data.frame
#   sample_id  age  sex  batch  outcome
#   Sample_1    45    M  Batch1       0
#   Sample_2    52    F  Batch1       1
#   Sample_3    38    M  Batch2       0
#   ...
```

- Rows = Samples
- Must include outcome variable
- Must include covariate columns
- Sample IDs should match count matrix column names

## ⚙️ Advanced Options

### Custom Confidence Levels

```r
# Use 99% confidence intervals
results <- olink_regression_analysis(
  count_matrix = data,
  metadata = metadata,
  outcome = "outcome",
  conf_level = 0.99  # Default is 0.95
)
```

### Different Multiple Testing Corrections

```r
# Bonferroni correction (most conservative)
results_bonf <- olink_regression_analysis(
  count_matrix = data,
  metadata = metadata,
  outcome = "outcome",
  p_adjust_method = "bonferroni"
)

# No correction (not recommended)
results_none <- olink_regression_analysis(
  count_matrix = data,
  metadata = metadata,
  outcome = "outcome",
  p_adjust_method = "none"
)
```

### Log Transformation Options

```r
# Transform during adjustment, keep log scale
adj1 <- adjust_olink_covariates(
  count_matrix = data,
  metadata = metadata,
  covariates = c("batch"),
  log_transform = TRUE,
  return_log = TRUE
)

# Transform and back-transform to original scale
adj2 <- adjust_olink_covariates(
  count_matrix = data,
  metadata = metadata,
  covariates = c("batch"),
  log_transform = TRUE,
  return_log = FALSE
)
```

## 🐛 Troubleshooting

### Issue: No proteins successfully analyzed

**Cause:** Insufficient complete cases after removing missing data

**Solution:**
```r
# Check for missing data
sum(is.na(count_matrix))
sum(is.na(metadata$outcome))

# Lower minimum sample requirement
results <- olink_regression_analysis(
  ...,
  min_samples = 5  # Default is 10
)
```

### Issue: Separation or convergence warnings

**Cause:** Perfect or near-perfect separation in logistic regression

**Solution:**
- Check for proteins with extreme values
- Consider log transformation
- Add more samples if possible
- Review covariate balance

### Issue: Sample ID mismatch

**Cause:** Sample IDs don't match between matrix and metadata

**Solution:**
```r
# Check overlap
common <- intersect(colnames(count_matrix), metadata$sample_id)
length(common)

# Ensure IDs match exactly (case-sensitive)
colnames(count_matrix)[1:5]
metadata$sample_id[1:5]
```

## 📚 Additional Resources

- **Full Examples:** [olink_examples/README.md](olink_examples/README.md)
- **Test Files:** `tests/testthat/test-olink-*.R`
- **Source Code:**
  - `R/olink-covariate-adjustment.R`
  - `R/olink-regression-analysis.R`
- **Generate Examples:** Run `source("generate_olink_examples.R")`

## 📄 Citation

If you use these functions in your research, please cite:

```bibtex
@software{mestools2026,
  title = {mestools: Microbiome and Epidemiology Statistical Tools},
  author = {Tran, Van Hung},
  year = {2026},
  url = {https://github.com/vanhungtran/mestools}
}
```

## 📧 Support

For questions, issues, or feature requests:
- **GitHub Issues:** https://github.com/vanhungtran/mestools/issues
- **Email:** [your email]

---

**Version:** 1.0.0
**Last Updated:** 2026-01-14
**R Version Required:** ≥ 4.0.0
**Dependencies:** ggplot2, tidyr (for visualization)
