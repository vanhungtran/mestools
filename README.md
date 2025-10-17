# mestools <img src="man/figures/logo.png" align="right" height="139" alt="mestools logo" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/vanhungtran/mestools/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/vanhungtran/mestools/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/vanhungtran/mestools/branch/main/graph/badge.svg)](https://codecov.io/gh/vanhungtran/mestools?branch=main)
[![CRAN status](https://www.r-pkg.org/badges/version/mestools)](https://CRAN.R-project.org/package=mestools)
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
<!-- badges: end -->

**mestools** is a practical collection of helpers that we keep reaching for when shipping data products. The package currently focuses on three themes that match the functions exported in [`NAMESPACE`](NAMESPACE):

| Category | When to reach for it | Key helpers |
| --- | --- | --- |
| **Data utilities** | Understand or fabricate tabular data quickly | `quick_summary()`, `read_file_safe()`, `batch_apply()`, `generate_random_df()` |
| **Project & package workflows** | Bootstrap an analysis repo or automate checks before pushing | `create_project_structure()`, `check_dependencies()`, `install_and_load()`, `validate_github_repo()`, `deploy_package()` |
| **GEO dataset integrations** | Work with Gene Expression Omnibus (GEO) accessions from R | `read_geo_dataset()`, `read_multiple_geo_datasets()`, `get_geo_summary()`, `get_default_gse_list()`, `process_all_geo_datasets()` |

Every helper is implemented under `R/`, documented in `man/`, and exercised in `tests/` so you can trust the examples below.

## Installation

Install the development version from GitHub with your preferred tooling:

```r
# install.packages("pak")
pak::pak("vanhungtran/mestools")
```

or

```r
# install.packages("devtools")
devtools::install_github("vanhungtran/mestools")
```

> **Heads-up**
> Some features rely on optional packages and external services:
>
> - GEO helpers call [`GEOquery`](https://bioconductor.org/packages/release/bioc/html/GEOquery.html) (via Bioconductor) and need an internet connection.
> - Deployment automation layers on `devtools`, `usethis`, `gert`, `here`, and `gh` for git operations.
> - `install_and_load()` will install from CRAN, Bioconductor (`BiocManager`), or GitHub (`remotes`, `gh`).

## Quick start

```r
library(mestools)

# Create a reproducible project scaffold
create_project_structure("my_analysis")

# Simulate some data and explore it
df <- generate_random_df(n_rows = 250, n_cols = 6, seed = 42)
quick_summary(df)

# Read a batch of files with consistent error handling
inputs <- list.files("my_analysis/inputs", full.names = TRUE)
results <- batch_apply(inputs, read_file_safe)

# Check that key packages are ready before you deploy
check_dependencies(c("ggplot2", "dplyr", "purrr"))
```

## Data utilities

### `quick_summary()` — instant diagnostics
Summarise a data frame without pulling in heavier dependencies. The helper returns dimensions, column classes, missing value counts, and approximate size.

```r
summary_info <- quick_summary(mtcars)
summary_info$dimensions
summary_info$missing_values
```

### `read_file_safe()` — format-aware loading
Load CSV, TSV, RDS, or plain-text files with informative errors. Unsupported types throw early so you can handle them upstream.

```r
csv_data <- read_file_safe("data/metrics.csv")
notes <- read_file_safe("docs/notes.txt")
```

### `batch_apply()` — resilient batch processing
Apply a function over a vector or list, collect errors instead of crashing, and (optionally) show a progress bar.

```r
files <- c("sample1.csv", "sample2.csv", "sample3.csv")
parsed <- batch_apply(files, read_file_safe, type = "csv", .progress = TRUE)
```

### `generate_random_df()` — reproducible demo data
Spin up a mixed-type data frame (numeric, character, factor columns) that is handy for tests or walkthroughs.

```r
test_data <- generate_random_df(n_rows = 500, n_cols = 9, seed = 123)
```

## Project and package workflows

### `create_project_structure()` — consistent scaffolds
Start new analysis work with the same folder layout every time.

```r
create_project_structure(
  path = "client-deliverable",
  directories = c("data", "analysis", "reports", "figures")
)
```

### Dependency helpers — `check_dependencies()` and `install_and_load()`
Verify that required packages are installed (optionally installing them), or iterate through a mix of CRAN, Bioconductor, and GitHub packages with one call.

```r
check_dependencies(c("usethis", "gert", "here"), auto_install = FALSE)
install_and_load(c("ggplot2", "BiocManager", "tidyverse/dplyr"))
```

### Deployment helpers — `validate_github_repo()` and `deploy_package()`
Validate that a target GitHub repository exists before running the automated deployment workflow. `deploy_package()` rebuilds documentation, runs tests/checks (when requested), commits via `gert`, and pushes to the configured remote.

```r
if (validate_github_repo("https://github.com/vanhungtran/mestools.git")) {
  deploy_package(run_tests = TRUE, run_checks = FALSE)
}
```

## Working with GEO datasets

The GEO helpers wrap `GEOquery` to simplify downloading, summarising, and batching Gene Expression Omnibus accessions. Network calls can take a while; use the summary tools or subsets when exploring.

```r
# Ensure GEOquery is available first
if (!requireNamespace("GEOquery", quietly = TRUE)) {
  BiocManager::install("GEOquery")
}

# Download and process a single series
geo <- read_geo_dataset("GSE102628")
str(geo$expression_matrix)

# Summarise multiple series without full downloads
summary_tbl <- get_geo_summary(c("GSE102628", "GSE102641"))

# Batch download using the curated defaults (with polite pauses)
subset <- get_default_gse_list()[1:3]
geo_list <- read_multiple_geo_datasets(subset)

# Fire-and-forget processing of the full catalogue (can take hours!)
# process_all_geo_datasets(destdir = "geo_data")
```

See [GEO_FUNCTIONS_README.md](GEO_FUNCTIONS_README.md) and [GEO_IMPLEMENTATION_SUMMARY.md](GEO_IMPLEMENTATION_SUMMARY.md) for more detailed walkthroughs, and `test_geo_functions*.R` for reproducible scripts.

## Additional documentation & support

- Browse NEWS in [NEWS.md](NEWS.md) for release notes.
- Function reference pages (in `man/`) are generated with roxygen2.
- Example GEO scripts live in [`test_geo_functions.R`](test_geo_functions.R) and [`test_geo_functions_simple.R`](test_geo_functions_simple.R).

## Contributing

We welcome improvements! See [CONTRIBUTING.md](CONTRIBUTING.md) for the contribution guide, coding standards, and release checklist.

### Local development tips

```r
# Install development dependencies
devtools::install_dev_deps()

# Run automated checks
devtools::test()
devtools::check()
```

## License & citation

This project is released under the MIT License – see [LICENSE.md](LICENSE.md) for details.

To cite the package in publications:

```r
citation("mestools")
```

## Getting help

- Documentation: use `?function_name` once the package is loaded.
- Issues: report bugs or feature requests at [GitHub Issues](https://github.com/vanhungtran/mestools/issues).
- Questions: start a thread in [GitHub Discussions](https://github.com/vanhungtran/mestools/discussions).
