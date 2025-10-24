#!/usr/bin/env Rscript
# Test script for radar plot with confidence intervals

library(mestools)

# Create test data with clear group differences
set.seed(123)
n_per_group <- 20

test_data <- data.frame(
  group = rep(c("Control", "Treatment"), each = n_per_group),
  Quality_Score = c(rnorm(n_per_group, 60, 10), rnorm(n_per_group, 75, 12)),
  Accuracy = c(rnorm(n_per_group, 70, 8), rnorm(n_per_group, 80, 10)),
  Speed = c(rnorm(n_per_group, 65, 15), rnorm(n_per_group, 70, 12)),
  Reliability = c(rnorm(n_per_group, 75, 9), rnorm(n_per_group, 85, 11)),
  Efficiency = c(rnorm(n_per_group, 68, 11), rnorm(n_per_group, 78, 9))
)

cat("\n=== Test 1: ggplot engine with 95% CI (default) ===\n")
result1 <- radar_with_p(
  data = test_data,
  value_vars = c("Quality_Score", "Accuracy", "Speed", "Reliability", "Efficiency"),
  group_var = "group",
  mode = "groups",
  engine = "ggplot",
  p_method = "t.test",
  show_ci = TRUE,
  title = "Group Comparison with 95% CI"
)

cat("\nCI data (first few rows):\n")
print(head(result1$ci))
cat("\nCI level:", result1$ci_level, "\n")

cat("\n=== Test 2: ggplot engine with 99% CI ===\n")
result2 <- radar_with_p(
  data = test_data,
  value_vars = c("Quality_Score", "Accuracy", "Speed", "Reliability", "Efficiency"),
  group_var = "group",
  mode = "groups",
  engine = "ggplot",
  p_method = "t.test",
  show_ci = TRUE,
  ci_level = 0.99,
  ci_alpha = 0.2,
  title = "Group Comparison with 99% CI"
)

cat("\n=== Test 3: fmsb engine with 95% CI ===\n")
result3 <- radar_with_p(
  data = test_data,
  value_vars = c("Quality_Score", "Accuracy", "Speed", "Reliability", "Efficiency"),
  group_var = "group",
  mode = "groups",
  engine = "fmsb",
  p_method = "t.test",
  show_ci = TRUE,
  title = "Group Comparison (fmsb) with 95% CI"
)

cat("\n=== Test 4: Three groups with ANOVA and 95% CI ===\n")
test_data_3groups <- data.frame(
  group = rep(c("A", "B", "C"), each = 15),
  var1 = c(rnorm(15, 50, 10), rnorm(15, 60, 12), rnorm(15, 70, 11)),
  var2 = c(rnorm(15, 60, 9), rnorm(15, 65, 10), rnorm(15, 75, 13)),
  var3 = c(rnorm(15, 55, 11), rnorm(15, 70, 12), rnorm(15, 65, 10)),
  var4 = c(rnorm(15, 65, 13), rnorm(15, 75, 11), rnorm(15, 80, 14))
)

result4 <- radar_with_p(
  data = test_data_3groups,
  value_vars = c("var1", "var2", "var3", "var4"),
  group_var = "group",
  mode = "groups",
  engine = "ggplot",
  p_method = "anova",
  show_ci = TRUE,
  title = "Three Group Comparison with 95% CI"
)

cat("\n=== Test 5: Without confidence intervals (comparison) ===\n")
result5 <- radar_with_p(
  data = test_data,
  value_vars = c("Quality_Score", "Accuracy", "Speed", "Reliability", "Efficiency"),
  group_var = "group",
  mode = "groups",
  engine = "ggplot",
  p_method = "t.test",
  show_ci = FALSE,
  title = "Group Comparison without CI"
)

cat("\nCI data:", result5$ci, "\n")
cat("CI level:", result5$ci_level, "\n")

cat("\n=== Summary ===\n")
cat("All tests completed successfully!\n")
cat("- Test 1: ggplot with default 95% CI\n")
cat("- Test 2: ggplot with 99% CI and custom alpha\n")
cat("- Test 3: fmsb with 95% CI\n")
cat("- Test 4: Three groups with ANOVA and 95% CI\n")
cat("- Test 5: Without CI for comparison\n")
