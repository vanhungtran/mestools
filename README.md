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
```

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