# Test script for expanded 16S primer database
# Run with: Rscript test_expanded_primers.R

devtools::load_all()

cat("\n=== Testing Expanded 16S Primer Database ===\n\n")

# Test 1: Get all primers
cat("Test 1: Loading full primer database\n")
primers <- get_16s_primers()
cat("Total primers loaded:", nrow(primers), "\n")
cat("Columns:", paste(names(primers), collapse = ", "), "\n\n")

# Test 2: Print statistics
cat("Test 2: Primer Database Statistics\n")
print_primer_stats()

# Test 3: Search by region
cat("\n\nTest 3: Search for V4 primers\n")
v4_primers <- search_16s_primers(region = "V4")
if (nrow(v4_primers) > 0) {
  cat("Found V4 primers:\n")
  print(v4_primers[, c("forward_name", "reverse_name", "application", "reference")])
}

# Test 4: Search by amplicon size
cat("\n\nTest 4: Search for short amplicons (<300 bp)\n")
short_primers <- search_16s_primers(amplicon_size_max = 300)
if (nrow(short_primers) > 0) {
  cat("Short amplicon primers:\n")
  print(short_primers[, c("region", "amplicon_size", "application")])
}

# Test 5: Search by platform
cat("\n\nTest 5: Search for Illumina NovaSeq optimized primers\n")
illumina_primers <- search_16s_primers(platform = "NovaSeq")
if (nrow(illumina_primers) > 0) {
  cat("NovaSeq primers:\n")
  print(illumina_primers[, c("region", "forward_name", "reverse_name", "platform_optimized")])
}

# Test 6: Search by keyword
cat("\n\nTest 6: Search for marine-related primers\n")
marine_primers <- search_16s_primers(keyword = "marine")
if (nrow(marine_primers) > 0) {
  cat("Marine primers:\n")
  print(marine_primers[, c("region", "application", "notes")])
}

# Test 7: Search by reference
cat("\n\nTest 7: Search for Earth Microbiome Project primers\n")
emp_primers <- search_16s_primers(reference = "Caporaso")
if (nrow(emp_primers) > 0) {
  cat("EMP primers (Caporaso):\n")
  print(emp_primers[, c("forward_name", "reverse_name", "application", "reference")])
}

# Test 8: Complex search
cat("\n\nTest 8: Complex search (Very Low bias + Very High coverage)\n")
optimal_primers <- search_16s_primers(
  bias_level = "Very Low",
  coverage = "Very High"
)
if (nrow(optimal_primers) > 0) {
  cat("Optimal primers:\n")
  print(optimal_primers[, c("region", "forward_name", "reverse_name", "application")])
}

# Test 9: Get unique references
cat("\n\nTest 9: All references in database\n")
all_refs <- unique(primers$reference)
cat(paste(" -", all_refs), sep = "\n")

# Test 10: Summary by coverage
cat("\n\nTest 10: Primers by coverage rating\n")
coverage_summary <- table(primers$coverage)
print(coverage_summary)

cat("\n\n=== All Tests Complete ===\n")
cat("Primer database successfully expanded with", nrow(primers), "primer pairs\n")
cat("Spanning", length(unique(primers$region)), "unique regions\n")
cat("From", length(unique(primers$reference)), "different publications\n")
