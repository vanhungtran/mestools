# mestools <img src="man/figures/logo.png" align="right" height="139" alt="" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/vanhungtran/mestools/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/vanhungtran/mestools/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/vanhungtran/mestools/branch/main/graph/badge.svg)](https://codecov.io/gh/vanhungtran/mestools?branch=main)
[![CRAN status](https://www.r-pkg.org/badges/version/mestools)](https://CRAN.R-project.org/package=mestools)
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
<!-- badges: end -->

**mestools** is a comprehensive R package that provides utility functions for data manipulation, package development, and deployment automation. It's designed to streamline common analytical tasks and enhance developer productivity.

## Features

### 🔧 Data Analysis Tools
- **Quick data summaries** with comprehensive statistics
- **Safe file reading** with automatic format detection
- **Batch processing** with progress reporting and error handling
- **Random data generation** for testing and demonstrations

### 🧬 Genomics & Microbiome Analysis
- **GEO dataset downloading** - Batch download and process Gene Expression Omnibus datasets
- **16S rRNA primer tools** - Design and evaluate primers for microbiome sequencing
- **Primer coverage analysis** - Assess taxonomic coverage and specificity
- **In silico PCR** - Test primer pairs against reference databases

### 📦 Package Development
- **Automated GitHub deployment** with testing and documentation
- **Dependency checking** and installation management
- **Repository validation** using GitHub API

### 🏗️ Project Management
- **Standardized project structure** creation
- **Multi-source package installation** (CRAN, Bioconductor, GitHub)
- **Development workflow automation**

## Installation

You can install the development version of mestools from [GitHub](https://github.com/vanhungtran/mestools) with:

```r
# install.packages("pak")
pak::pak("vanhungtran/mestools")
```

Or using devtools:

```r
# install.packages("devtools")
devtools::install_github("vanhungtran/mestools")
```

## Quick Start

```r
library(mestools)

# Create a project structure
create_project_structure("my_analysis")

# Generate test data
df <- generate_random_df(1000, 8)

# Get quick summary
summary_info <- quick_summary(df)
print(summary_info)

# Check package dependencies
check_dependencies(c("ggplot2", "dplyr", "tidyr"))

# Download GEO datasets
geo_data <- read_geo_dataset("GSE102628")

# Get 16S rRNA primers
primers <- get_16s_primers()
v4_primers <- get_16s_primers(region = "V4")
```

## Key Features

### 🧬 Genomics & Microbiome Analysis (NEW!)

**GEO Dataset Tools:**
- Download and process Gene Expression Omnibus (GEO) datasets
- Batch processing of multiple datasets with error handling
- Extract expression matrices, phenotype data, and feature annotations
- 139 pre-curated GEO dataset IDs included
- Automatic retry and rate limiting for reliable downloads

**16S rRNA Primer Tools:**
- Access database of validated primer pairs for all 16S regions (V1-V9)
- Get primers optimized for specific applications (gut, soil, marine, etc.)
- Design custom primers for your sequences
- Evaluate primer coverage across bacterial/archaeal taxa
- Perform in silico PCR to test primer specificity
- Compare multiple primer pairs side-by-side
- Get recommendations based on sample type and sequencing platform

**Supported 16S Regions:**
- V1-V2, V1-V3, V3-V4, V4, V4-V5, V6-V8, V7-V9
- Earth Microbiome Project standard primers included
- Human Microbiome Project primers included

## Core Functions

### Data Analysis

#### `quick_summary()`
Get comprehensive data frame statistics:

```r
# Basic usage
data_summary <- quick_summary(mtcars)
print(data_summary$dimensions)  # [1] 32 11
print(data_summary$missing_values)  # Named vector of NA counts
```

#### `generate_random_df()`
Create random data for testing:

```r
# Generate test data
test_data <- generate_random_df(n_rows = 500, n_cols = 6, seed = 42)

# Creates mixed data types:
# - numeric1, numeric4: normal distributions  
# - character2, character5: random strings
# - factor3, factor6: factor variables
```

#### `read_file_safe()`
Safely read various file formats:

```r
# Automatic format detection
data <- read_file_safe("data.csv")
text <- read_file_safe("notes.txt") 
saved_object <- read_file_safe("results.rds")

# Manual type specification
data <- read_file_safe("file.dat", type = "csv")
```

#### `batch_apply()`
Process multiple objects with progress tracking:

```r
# Process list of data frames
file_list <- list("file1.csv", "file2.csv", "file3.csv")
results <- batch_apply(file_list, read_file_safe, .progress = TRUE)

# Apply function with error handling
results <- batch_apply(1:100, function(x) {
  if (x %% 10 == 0) stop("Divisible by 10!")
  return(x^2)
})
```

### Package Development

#### `deploy_package()`
Automated package deployment to GitHub:

```r
# Deploy with default settings
deploy_package()

# Custom deployment
deploy_package(
  repo_url = "https://github.com/username/mypackage.git",
  commit_message = "Major update with new features",
  run_tests = TRUE,
  run_checks = TRUE
)
```

#### `validate_github_repo()`
Check if a GitHub repository exists:

```r
# Validate repository
if (validate_github_repo("https://github.com/r-lib/usethis.git")) {
  message("Repository is accessible!")
}
```

#### `check_dependencies()`
Verify and install package dependencies:

```r
# Check if packages are available
check_dependencies(c("ggplot2", "dplyr", "tidyr"))

# Auto-install missing packages
check_dependencies(c("ggplot2", "dplyr", "tidyr"), auto_install = TRUE)
```

#### `install_and_load()`
Smart package installation from multiple sources:

```r
# Install from CRAN, Bioconductor, or GitHub
install_and_load(c(
  "ggplot2",                    # CRAN
  "DESeq2",                     # Bioconductor  
  "tidyverse/dplyr",           # Direct GitHub
  "some_package"               # GitHub search
))
```

### Project Management

#### `create_project_structure()`
Set up standardized project directories:

```r
# Default structure
create_project_structure("my_project")
# Creates: inputs/, scripts/, output/, docs/, ref/, reports/

# Custom directories
create_project_structure(
  "custom_project",
  directories = c("data", "analysis", "results", "figures")
)
```

### Genomics & Microbiome Tools

#### GEO Dataset Functions

Download and process Gene Expression Omnibus datasets:

```r
# Download a single GEO dataset
result <- read_geo_dataset("GSE102628")
expression_matrix <- result$expression_matrix
metadata <- result$phenotype_data

# Batch download multiple datasets
gse_list <- c("GSE102628", "GSE102641", "GSE102725")
results <- read_multiple_geo_datasets(gse_list)

# Get dataset summaries without downloading
summary_info <- get_geo_summary(c("GSE102628", "GSE102641"))

# Get predefined list of 139 datasets
all_datasets <- get_default_gse_list()

# Process all default datasets (or subset for testing)
test_results <- process_all_geo_datasets(subset_size = 5)
```

#### 16S rRNA Primer Tools

Design and evaluate primers for microbiome sequencing:

```r
# Get all available 16S primer pairs
all_primers <- get_16s_primers()

# Get primers for specific region
v4_primers <- get_16s_primers(region = "V4")
v3v4_primers <- get_16s_primers(region = "V3-V4")

# Get primers for specific application
gut_primers <- get_16s_primers(application = "gut")
environmental_primers <- get_16s_primers(application = "environmental")

# Design custom primers for your sequences
custom_primers <- design_16s_primers(
  sequences = my_16s_sequences,
  target_region = "V4",
  min_coverage = 0.95
)

# Evaluate primer coverage
evaluation <- evaluate_primer_coverage(
  forward_primer = "GTGCCAGCMGCCGCGGTAA",
  reverse_primer = "GGACTACHVGGGTWTCTAAT",
  reference_db = silva_database
)

# Test primers in silico
pcr_results <- test_primers_insilico(
  forward = "GTGCCAGCMGCCGCGGTAA",
  reverse = "GGACTACHVGGGTWTCTAAT",
  database = "silva",
  max_mismatches = 2
)

# Get primer recommendations
recommendations <- recommend_primers(
  sample_type = "human_gut",
  sequencing_platform = "illumina",
  read_length = 250
)
```

#### Differential Expression Visualization

Create publication-ready heatmaps with fold change and volcano plots:

```r
# Basic usage with significance annotations
plot_heatmap_with_fc(
  expression_matrix = expr_matrix,
  log2fc = fold_changes,
  pvalues = pvalues,
  sample_groups = groups,
  title = "Differential Expression"
)

# Advanced: Filter top genes with volcano plot
plot_heatmap_with_fc(
  expression_matrix = expr_matrix,
  log2fc = fold_changes,
  pvalues = pvalues,
  padj = adjusted_pvalues,
  top_n = 30,                          # Show top 30 genes
  highlight_genes = c("TP53", "BRCA1"), # Highlight specific genes
  color_theme = "viridis",              # Professional color theme
  plot_layout = "heatmap_fc_volcano",   # 3-panel layout
  title = "Top 30 Differentially Expressed Genes"
)

# Publication-ready figure
plot_heatmap_with_fc(
  expression_matrix = expr_matrix,
  log2fc = fold_changes,
  pvalues = pvalues,
  padj = adjusted_pvalues,
  top_n = 25,
  color_theme = "publication",  # Grayscale for journals
  gene_categories = pathways,   # Add pathway annotations
  title = "Differentially Expressed Genes"
)
```

**Key Features:**

- Statistical significance stars (*, **, ***)
- Gene filtering by fold change or top N
- Gene highlighting with custom colors
- 5 color themes (default, viridis, RdBu, publication, colorblind)
- Volcano plot integration
- Multi-factor annotations
- Synchronized 3-panel layouts

**Color Themes:**

- `"default"` - Blue-white-red (classic)
- `"viridis"` - Perceptually uniform (best for data visualization)
- `"RdBu"` - Red-Blue diverging (ColorBrewer standard)
- `"publication"` - Grayscale (journal-ready)
- `"colorblind"` - Accessible palette (inclusive design)

![Heatmap Example](man/figures/heatmap.png)

## Advanced Usage

### Workflow Integration

Combine functions for complete data analysis workflows:

```r
# 1. Set up project
create_project_structure("analysis_2024")

# 2. Check dependencies
required_pkgs <- c("tidyverse", "ggplot2", "plotly")
check_dependencies(required_pkgs, auto_install = TRUE)

# 3. Process multiple files
data_files <- list.files("data/", pattern = "\\.csv$", full.names = TRUE)
datasets <- batch_apply(data_files, read_file_safe)

# 4. Generate summaries
summaries <- batch_apply(datasets, quick_summary)

# 5. Create test data if needed
test_data <- generate_random_df(1000, 10)
```

### Microbiome Study Workflow

Complete workflow for 16S rRNA microbiome analysis:

```r
# 1. Set up microbiome project
create_project_structure("microbiome_study_2024")

# 2. Get recommended primers for your sample type
primers <- recommend_primers(
  sample_type = "human_gut",
  sequencing_platform = "illumina",
  read_length = 250
)

# 3. Evaluate primer coverage
coverage <- evaluate_primer_coverage(
  forward_primer = primers$forward_seq,
  reverse_primer = primers$reverse_seq,
  reference_db = "silva"
)

# 4. Test primers in silico
pcr_test <- test_primers_insilico(
  forward = primers$forward_seq,
  reverse = primers$reverse_seq,
  database = "silva",
  max_mismatches = 2
)

# 5. Compare multiple primer pairs
comparison <- compare_primer_pairs(
  primer_pairs = list(
    V3V4 = c("CCTACGGGNGGCWGCAG", "GACTACHVGGGTATCTAATCC"),
    V4 = c("GTGCCAGCMGCCGCGGTAA", "GGACTACHVGGGTWTCTAAT")
  ),
  reference_db = "silva"
)
```

### GEO Data Analysis Workflow

Download and analyze multiple GEO datasets:

```r
# 1. Get list of relevant datasets
gse_list <- get_default_gse_list()
gut_datasets <- gse_list[grep("gut|intestin", gse_list, ignore.case = TRUE)]

# 2. Download datasets
datasets <- read_multiple_geo_datasets(
  gse_ids = gut_datasets,
  destdir = "geo_data",
  sleep_between = 2
)

# 3. Process expression data
for (gse_id in names(datasets)) {
  expr_matrix <- datasets[[gse_id]]$expression_matrix
  phenotype <- datasets[[gse_id]]$phenotype_data
  
  # Your analysis here
  cat("Processing", gse_id, ":", 
      nrow(expr_matrix), "features x", 
      ncol(expr_matrix), "samples\n")
}

# 4. Save processed results
saveRDS(datasets, "processed_geo_datasets.rds")
```

### Package Development Workflow

Streamline your R package development:

```r
# 1. Check dependencies
check_dependencies(c("devtools", "usethis", "testthat"))

# 2. Validate target repository
validate_github_repo("https://github.com/username/mypackage.git")

# 3. Deploy package
deploy_package(
  commit_message = "Release version 1.0.0",
  run_tests = TRUE,
  run_checks = TRUE
)
```

## Error Handling

All functions include comprehensive error handling:

```r
# Functions fail gracefully with informative messages
try({
  read_file_safe("nonexistent.csv")
})
# Error: File does not exist: nonexistent.csv

# Batch operations continue despite individual failures
results <- batch_apply(c("good.csv", "bad.csv", "good2.csv"), read_file_safe)
# Warnings for failed operations, but processing continues
```

## Configuration

### Options

Set package-specific options:

```r
# Enable test mode (useful for package development)
options(mestools.test_mode = TRUE)

# Disable progress bars globally  
options(mestools.progress = FALSE)
```

### Dependencies

**Required packages:**
- `devtools`, `usethis`, `gert`, `here` (for deployment functions)
- `utils`, `tools`, `stats` (base functionality)

**Suggested packages:**
- `gh` (for GitHub repository validation)
- `remotes` (for GitHub package installation)
- `BiocManager` (for Bioconductor packages)
- `GEOquery` (for GEO dataset downloading)
- `Biostrings` (for 16S sequence analysis)
- `DECIPHER` (for primer design and evaluation)

## Contributing

Contributions are welcome! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for details.

### Development Setup

```r
# Clone repository
git clone https://github.com/vanhungtran/mestools.git

# Install development dependencies
devtools::install_dev_deps()

# Run tests
devtools::test()

# Check package
devtools::check()
```

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details.

## Citation

```r
citation("mestools")
```

## Changelog

See [NEWS.md](NEWS.md) for details about changes in each version.

## Getting Help

- **Documentation**: Use `?function_name` for help on specific functions
- **Issues**: Report bugs at [GitHub Issues](https://github.com/vanhungtran/mestools/issues)
- **Discussions**: Ask questions at [GitHub Discussions](https://github.com/vanhungtran/mestools/discussions)

## Related Packages

- [`usethis`](https://usethis.r-lib.org/) - Workflow functions for package development
- [`devtools`](https://devtools.r-lib.org/) - Package development tools
- [`here`](https://here.r-lib.org/) - Reproducible file paths
- [`gert`](https://docs.ropensci.org/gert/) - Simple git client

---

*Developed by [Lucas VHH TRAN](https://github.com/vanhungtran) • 2025*