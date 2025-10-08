# Contributing to mestools

Thank you for your interest in contributing to mestools! This document provides guidelines for contributing to the project.

## Types of Contributions

### 🐛 Bug Reports
- Use the [GitHub issue tracker](https://github.com/vanhungtran/mestools/issues)
- Include a minimal reproducible example
- Specify your R version and operating system
- Include relevant error messages and warnings

### 💡 Feature Requests
- Open a [GitHub issue](https://github.com/vanhungtran/mestools/issues) with the "enhancement" label
- Clearly describe the proposed feature and its use case
- Explain how it fits with the package's goals

### 🔧 Code Contributions
- Fork the repository
- Create a feature branch
- Make your changes with tests
- Submit a pull request

## Development Setup

### Prerequisites
- R >= 4.0.0
- RStudio (recommended)
- Git

### Getting Started

```r
# 1. Fork and clone the repository
git clone https://github.com/yourusername/mestools.git
cd mestools

# 2. Install development dependencies
R -e "devtools::install_dev_deps()"

# 3. Load the package for development
R -e "devtools::load_all()"
```

### Development Workflow

1. **Create a branch** for your feature or bug fix:
   ```bash
   git checkout -b feature/your-feature-name
   ```

2. **Make your changes** following the coding standards below

3. **Write tests** for new functionality:
   ```r
   # Add tests to tests/testthat/test-*.R
   devtools::test()
   ```

4. **Update documentation**:
   ```r
   # Update roxygen comments and rebuild docs
   devtools::document()
   ```

5. **Check the package**:
   ```r
   devtools::check()
   ```

6. **Commit and push**:
   ```bash
   git add .
   git commit -m "feat: add new feature description"
   git push origin feature/your-feature-name
   ```

7. **Create a pull request**

## Coding Standards

### R Code Style
Follow the [tidyverse style guide](https://style.tidyverse.org/):

```r
# Good
create_project_structure <- function(path = ".", 
                                     directories = c("inputs", "scripts")) {
  # Function body
}

# Bad
create_project_structure<-function(path=".",directories=c("inputs","scripts")){
  # Function body
}
```

### Function Documentation
Use roxygen2 for documentation:

```r
#' Brief Function Description
#'
#' Longer description of what the function does, including details
#' about the algorithm or approach if relevant.
#'
#' @param param1 Description of parameter 1
#' @param param2 Description of parameter 2
#' @return Description of return value
#' @export
#' @examples
#' \dontrun{
#' result <- your_function(param1 = "value")
#' }
your_function <- function(param1, param2 = NULL) {
  # Function implementation
}
```

### Error Handling
Use informative error messages:

```r
# Good
if (!is.data.frame(data)) {
  stop("Input must be a data frame, got: ", class(data)[1])
}

# Bad  
if (!is.data.frame(data)) {
  stop("Bad input")
}
```

### Progress Reporting
Use consistent progress reporting:

```r
if (.progress) {
  message("Processing ", length(objects), " items...")
  pb <- utils::txtProgressBar(min = 0, max = length(objects), style = 3)
}
```

## Testing Guidelines

### Test Structure
- Place tests in `tests/testthat/test-*.R` files
- Use descriptive test names
- Group related tests together

### Test Examples

```r
test_that("function_name handles normal cases", {
  result <- function_name(input)
  expect_equal(result, expected_output)
  expect_type(result, "list")
})

test_that("function_name handles edge cases", {
  expect_error(function_name(NULL), "Input cannot be NULL")
  expect_warning(function_name(problematic_input))
})

test_that("function_name handles empty inputs", {
  result <- function_name(character(0))
  expect_equal(length(result), 0)
})
```

### Test Coverage
- Aim for >90% test coverage
- Test both success and failure cases
- Include edge cases and boundary conditions
- Test error messages and warnings

### Running Tests

```r
# Run all tests
devtools::test()

# Run specific test file
devtools::test_file("tests/testthat/test-utils.R")

# Check test coverage
covr::package_coverage()
```

## Documentation

### Function Documentation
- Every exported function must have complete roxygen2 documentation
- Include working examples (use `\dontrun{}` for functions requiring external resources)
- Document all parameters and return values
- Provide meaningful descriptions

### Package Documentation
- Update `README.md` for user-facing changes
- Update `NEWS.md` for all changes
- Keep vignettes current with new features

### Examples
- Provide realistic, working examples
- Use built-in datasets when possible
- Show common use cases
- Include error handling examples

## Pull Request Process

1. **Pre-submission checklist**:
   - [ ] All tests pass (`devtools::test()`)
   - [ ] Package check passes (`devtools::check()`)
   - [ ] Documentation is updated
   - [ ] NEWS.md is updated
   - [ ] Code follows style guidelines

2. **Pull Request Description**:
   - Clearly describe what changes were made
   - Reference any related issues
   - Include examples of new functionality
   - Mention any breaking changes

3. **Review Process**:
   - Maintainers will review your code
   - Address any requested changes
   - Keep discussions constructive and professional

## Issue Reporting

### Bug Reports
Include:
- **Brief description** of the problem
- **Reproducible example** (minimal code that demonstrates the issue)
- **Expected behavior** vs **actual behavior**
- **Environment details**: R version, OS, package version
- **Error messages** (full output, not screenshots)

### Feature Requests
Include:
- **Clear description** of the desired feature
- **Use case** or problem it solves
- **Proposed API** (if you have ideas)
- **Examples** of how it would be used

## Code of Conduct

### Our Standards
- Be respectful and inclusive
- Focus on constructive feedback
- Help create a welcoming environment
- Be patient with newcomers

### Enforcement
Unacceptable behavior may result in temporary or permanent ban from the project.

## Getting Help

- **Documentation**: Check function help with `?function_name`
- **Issues**: Search existing issues before creating new ones
- **Discussions**: Use GitHub Discussions for questions
- **Email**: Contact maintainers for sensitive issues

## Recognition

Contributors will be recognized in:
- Package `DESCRIPTION` file
- `AUTHORS.md` file
- Release notes
- Package documentation

Thank you for contributing to mestools! 🚀