# mestools (development version)

## mestools 0.1.0

### New Features

* **Data Analysis Functions**
  * `quick_summary()` - Comprehensive data frame analysis with dimensions, types, missing values, and memory usage
  * `generate_random_df()` - Create random test data with mixed column types (numeric, character, factor)
  * `read_file_safe()` - Safe file reading with automatic format detection (CSV, TSV, RDS, TXT)
  * `batch_apply()` - Batch processing with progress reporting and error handling

* **Package Development Tools**
  * `deploy_package()` - Complete automated deployment workflow to GitHub with testing and documentation
  * `validate_github_repo()` - GitHub repository existence validation using GitHub API
  * `check_dependencies()` - Package dependency verification with optional auto-installation

* **Project Management**
  * `create_project_structure()` - Standardized project directory structure creation
  * `install_and_load()` - Multi-source package installation (CRAN, Bioconductor, GitHub) with interactive GitHub search

### Package Infrastructure

* Comprehensive test suite with 90+ tests covering:
  * Core functionality testing
  * Edge case handling
  * Integration workflows
  * Performance validation
  * Error handling verification

* Complete documentation with:
  * Function help pages with examples
  * Package vignettes
  * README with usage examples
  * Development guidelines

### Dependencies

* **Imports**: `devtools`, `usethis`, `gert`, `here`, `utils`, `tools`, `stats`, `desc`, `knitr`, `rmarkdown`, `testthat`, `pkgdown`
* **Suggests**: `gh`, `remotes`, `BiocManager`, `covr`, `curl`, `pkgload`, `spelling`

### Technical Details

* **R Version**: Requires R >= 4.0.0
* **License**: MIT + file LICENSE
* **Encoding**: UTF-8
* **Roxygen**: 7.3.2
* **Test Framework**: testthat edition 3