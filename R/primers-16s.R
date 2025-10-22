# 16S rRNA Primer Functions for Microbiome Sequencing
# Author: mestools package
# Date: 2025-10-22

# --- 16S rRNA Primer Database and Functions ---

#' Get Standard 16S rRNA Primer Pairs
#'
#' Returns a comprehensive database of commonly used 16S rRNA primer pairs
#' for microbiome sequencing studies.
#'
#' @return A data.frame containing primer information with columns:
#'   region, name, forward_primer, reverse_primer, forward_name, reverse_name, 
#'   application, amplicon_size
#' @export
#' @examples
#' # Get all standard primer pairs
#' primers <- get_16s_primers()
#' print(primers)
#' 
#' # Filter for V4 region
#' v4_primers <- primers[primers$region == "V4", ]
get_16s_primers <- function() {
  primers_df <- data.frame(
    region = c("V1-V2", "V1-V3", "V3-V4", "V4", "V4-V5", "V6-V8", "V7-V9"),
    forward_name = c("27F", "27F", "341F", "515F", "515F", "939F", "1115F"),
    forward_primer = c(
      "AGAGTTTGATCMTGGCTCAG",
      "AGAGTTTGATCMTGGCTCAG",
      "CCTACGGGNGGCWGCAG",
      "GTGCCAGCMGCCGCGGTAA",
      "GTGCCAGCMGCCGCGGTAA",
      "TTGACGGGGGCCCGCAC",
      "GCAACGAGCGCAACCC"
    ),
    reverse_name = c("338R", "534R", "785R", "806R", "944R", "1378R", "1492R"),
    reverse_primer = c(
      "TGCTGCCTCCCGTAGGAGT",
      "ATTACCGCGGCTGCTGG",
      "GACTACHVGGGTATCTAATCC",
      "GGACTACHVGGGTWTCTAAT",
      "GAATTGGCACCTTGAGGAC",
      "CGGTGTGTACAAGGCCCGGGAACG",
      "GGTTACCTTGTTACGACTT"
    ),
    application = c(
      "Broad bacterial profiling",
      "Environmental microbiomes",
      "Gut microbiota (widely used)",
      "Universal (Earth Microbiome Project)",
      "Marine/environmental studies",
      "Archaeal and bacterial detection",
      "Full-length near-universal"
    ),
    amplicon_size = c(311, 507, 444, 291, 429, 439, 377),
    coverage = c("High", "High", "Very High", "Very High", "High", "Medium-High", "High"),
    bias_level = c("Low", "Low", "Very Low", "Very Low", "Low", "Medium", "Low"),
    stringsAsFactors = FALSE
  )
  
  return(primers_df)
}

#' Select Optimal 16S Primer Pair
#'
#' Recommends the best 16S rRNA primer pair based on study requirements.
#'
#' @param application Character string describing the study type: 
#'   "gut", "soil", "marine", "environmental", "universal", "archaeal"
#' @param read_length Integer. Sequencing read length (e.g., 150, 250, 300)
#' @param platform Character. Sequencing platform: "illumina", "pacbio", "nanopore"
#' @param region Character. Preferred hypervariable region (optional). E.g., "V3-V4", "V4"
#' @return A list containing recommended primer information and rationale
#' @export
#' @examples
#' \dontrun{
#' # Get recommendation for gut microbiome study
#' recommendation <- select_16s_primers(
#'   application = "gut",
#'   read_length = 250,
#'   platform = "illumina"
#' )
#' print(recommendation)
#' 
#' # For marine environment
#' marine_primers <- select_16s_primers(
#'   application = "marine",
#'   read_length = 300,
#'   platform = "illumina"
#' )
#' }
select_16s_primers <- function(application = "universal", 
                               read_length = 250,
                               platform = "illumina",
                               region = NULL) {
  
  primers_db <- get_16s_primers()
  
  # Application-specific recommendations
  recommendations <- list(
    gut = list(regions = c("V3-V4", "V4"), rationale = "Best resolution for gut bacteria"),
    soil = list(regions = c("V4", "V3-V4"), rationale = "High bacterial diversity coverage"),
    marine = list(regions = c("V4-V5", "V4"), rationale = "Good for marine bacteria and archaea"),
    environmental = list(regions = c("V4", "V1-V3"), rationale = "Broad environmental coverage"),
    universal = list(regions = c("V4", "V3-V4"), rationale = "Most widely validated"),
    archaeal = list(regions = c("V6-V8", "V4"), rationale = "Enhanced archaeal detection"),
    full_length = list(regions = c("V7-V9", "V1-V2"), rationale = "Near full-length coverage")
  )
  
  # Get application recommendation
  app_lower <- tolower(application)
  if (!app_lower %in% names(recommendations)) {
    message("Unknown application. Defaulting to 'universal'")
    app_lower <- "universal"
  }
  
  rec <- recommendations[[app_lower]]
  
  # If region is specified, use it; otherwise use recommendation
  if (!is.null(region)) {
    selected_region <- region
  } else {
    selected_region <- rec$regions[1]
  }
  
  # Get primer pair for selected region
  primer_info <- primers_db[primers_db$region == selected_region, ]
  
  if (nrow(primer_info) == 0) {
    # Fall back to V4 (most universal)
    message("Region not found. Defaulting to V4 region.")
    primer_info <- primers_db[primers_db$region == "V4", ]
  }
  
  # Platform-specific considerations
  platform_notes <- switch(tolower(platform),
    illumina = "Paired-end sequencing recommended. Consider index adapters.",
    pacbio = "Suitable for long-read sequencing. Full-length 16S possible.",
    nanopore = "Long-read capable. Consider V7-V9 for full coverage.",
    "Standard short-read platform assumed."
  )
  
  # Check if amplicon size is suitable for read length
  amplicon_suitable <- primer_info$amplicon_size <= (read_length * 2 - 50)
  coverage_warning <- if (!amplicon_suitable) {
    paste0("⚠️ WARNING: Amplicon size (", primer_info$amplicon_size, 
           "bp) may not be fully covered by ", read_length, "bp reads")
  } else {
    paste0("✅ Amplicon size (", primer_info$amplicon_size, 
           "bp) suitable for ", read_length, "bp paired-end reads")
  }
  
  result <- list(
    recommendation = selected_region,
    application = application,
    forward_primer = primer_info$forward_primer,
    forward_name = primer_info$forward_name,
    reverse_primer = primer_info$reverse_primer,
    reverse_name = primer_info$reverse_name,
    amplicon_size = primer_info$amplicon_size,
    rationale = rec$rationale,
    coverage = primer_info$coverage,
    bias_level = primer_info$bias_level,
    platform_notes = platform_notes,
    coverage_warning = coverage_warning,
    sequencing_setup = list(
      read_length = read_length,
      platform = platform,
      recommended_overlap = max(0, (read_length * 2) - primer_info$amplicon_size)
    )
  )
  
  class(result) <- c("primer_recommendation", "list")
  return(result)
}

#' Print Primer Recommendation
#'
#' @param x A primer_recommendation object
#' @param ... Additional arguments (unused)
#' @export
print.primer_recommendation <- function(x, ...) {
  cat("\n🧬 16S rRNA Primer Recommendation\n")
  cat("═══════════════════════════════════════════════\n\n")
  
  cat("📍 Target Region: ", x$recommendation, "\n")
  cat("🎯 Application: ", x$application, "\n")
  cat("💡 Rationale: ", x$rationale, "\n\n")
  
  cat("🔬 PRIMER DETAILS:\n")
  cat("─────────────────────────────────────────────\n")
  cat("Forward Primer (", x$forward_name, "):\n")
  cat("  5'-", x$forward_primer, "-3'\n\n")
  cat("Reverse Primer (", x$reverse_name, "):\n")
  cat("  5'-", x$reverse_primer, "-3'\n\n")
  
  cat("📏 Amplicon Size: ", x$amplicon_size, " bp\n")
  cat("📊 Coverage: ", x$coverage, "\n")
  cat("⚖️  Bias Level: ", x$bias_level, "\n\n")
  
  cat("🔧 SEQUENCING SETUP:\n")
  cat("─────────────────────────────────────────────\n")
  cat("Platform: ", x$sequencing_setup$platform, "\n")
  cat("Read Length: ", x$sequencing_setup$read_length, " bp\n")
  cat(x$coverage_warning, "\n")
  cat("\n💭 Note: ", x$platform_notes, "\n\n")
  
  invisible(x)
}

#' Design Custom 16S Primers
#'
#' Provides guidance for designing custom 16S rRNA primers based on sequence input.
#' This function helps evaluate primer properties and suggests improvements.
#'
#' @param forward_primer Character. Forward primer sequence (5' to 3')
#' @param reverse_primer Character. Reverse primer sequence (5' to 3')
#' @param check_degeneracy Logical. Whether to check for degenerate bases
#' @return A list containing primer analysis results
#' @export
#' @examples
#' # Analyze custom primer pair
#' analysis <- design_custom_16s_primers(
#'   forward_primer = "GTGCCAGCMGCCGCGGTAA",
#'   reverse_primer = "GGACTACHVGGGTWTCTAAT"
#' )
#' print(analysis)
design_custom_16s_primers <- function(forward_primer, 
                                     reverse_primer,
                                     check_degeneracy = TRUE) {
  
  # Convert to uppercase
  forward_primer <- toupper(forward_primer)
  reverse_primer <- toupper(reverse_primer)
  
  # Basic primer statistics
  analyze_primer <- function(seq) {
    length <- nchar(seq)
    
    # Count nucleotides
    bases <- strsplit(seq, "")[[1]]
    gc_count <- sum(bases %in% c("G", "C"))
    gc_content <- (gc_count / length) * 100
    
    # Check for degenerate bases
    degenerate_bases <- c("R", "Y", "S", "W", "K", "M", "B", "D", "H", "V", "N")
    deg_count <- sum(bases %in% degenerate_bases)
    has_degeneracy <- deg_count > 0
    
    # Calculate melting temperature (basic formula)
    # Tm = 4(G+C) + 2(A+T) for short primers (<14bp)
    # More accurate: 64.9 + 41*(G+C-16.4)/length for longer primers
    if (length < 14) {
      at_count <- sum(bases %in% c("A", "T"))
      tm <- 4 * gc_count + 2 * at_count
    } else {
      tm <- 64.9 + 41 * (gc_count - 16.4) / length
    }
    
    # Check for runs (homopolymers)
    has_poly_run <- grepl("(A{4,}|T{4,}|G{4,}|C{4,})", seq)
    
    # Check for 3' stability (last 5 bases)
    three_prime <- substr(seq, nchar(seq) - 4, nchar(seq))
    three_prime_gc <- sum(strsplit(three_prime, "")[[1]] %in% c("G", "C"))
    three_prime_stable <- three_prime_gc >= 2 && three_prime_gc <= 3
    
    return(list(
      length = length,
      gc_content = round(gc_content, 1),
      tm = round(tm, 1),
      has_degeneracy = has_degeneracy,
      degenerate_count = deg_count,
      has_poly_run = has_poly_run,
      three_prime_stable = three_prime_stable,
      sequence = seq
    ))
  }
  
  forward_stats <- analyze_primer(forward_primer)
  reverse_stats <- analyze_primer(reverse_primer)
  
  # Quality checks and recommendations
  issues <- character()
  recommendations <- character()
  
  # Length check (optimal 18-25 bp)
  if (forward_stats$length < 18 || forward_stats$length > 30) {
    issues <- c(issues, paste0("Forward primer length (", forward_stats$length, 
                               "bp) outside optimal range (18-25bp)"))
  }
  if (reverse_stats$length < 18 || reverse_stats$length > 30) {
    issues <- c(issues, paste0("Reverse primer length (", reverse_stats$length, 
                               "bp) outside optimal range (18-25bp)"))
  }
  
  # GC content check (optimal 40-60%)
  if (forward_stats$gc_content < 40 || forward_stats$gc_content > 60) {
    issues <- c(issues, paste0("Forward primer GC content (", forward_stats$gc_content, 
                               "%) outside optimal range (40-60%)"))
    recommendations <- c(recommendations, "Consider adjusting primer to achieve 40-60% GC content")
  }
  if (reverse_stats$gc_content < 40 || reverse_stats$gc_content > 60) {
    issues <- c(issues, paste0("Reverse primer GC content (", reverse_stats$gc_content, 
                               "%) outside optimal range (40-60%)"))
  }
  
  # Tm difference check (should be within 5°C)
  tm_diff <- abs(forward_stats$tm - reverse_stats$tm)
  if (tm_diff > 5) {
    issues <- c(issues, paste0("Tm difference (", round(tm_diff, 1), 
                               "°C) exceeds recommended 5°C"))
    recommendations <- c(recommendations, "Adjust primer lengths to match Tm values")
  }
  
  # Homopolymer runs
  if (forward_stats$has_poly_run) {
    issues <- c(issues, "Forward primer contains homopolymer runs (4+ same bases)")
    recommendations <- c(recommendations, "Avoid homopolymer runs to prevent sequencing errors")
  }
  if (reverse_stats$has_poly_run) {
    issues <- c(issues, "Reverse primer contains homopolymer runs")
  }
  
  # 3' stability
  if (!forward_stats$three_prime_stable) {
    issues <- c(issues, "Forward primer 3' end may lack optimal stability")
    recommendations <- c(recommendations, "Ensure 2-3 G/C bases in last 5 positions of 3' end")
  }
  if (!reverse_stats$three_prime_stable) {
    issues <- c(issues, "Reverse primer 3' end may lack optimal stability")
  }
  
  # Overall quality score
  quality_score <- 100
  quality_score <- quality_score - length(issues) * 10
  quality_score <- max(0, min(100, quality_score))
  
  quality_rating <- if (quality_score >= 80) {
    "Excellent"
  } else if (quality_score >= 60) {
    "Good"
  } else if (quality_score >= 40) {
    "Acceptable"
  } else {
    "Needs Improvement"
  }
  
  result <- list(
    forward = forward_stats,
    reverse = reverse_stats,
    tm_difference = round(tm_diff, 1),
    quality_score = quality_score,
    quality_rating = quality_rating,
    issues = if (length(issues) > 0) issues else "No issues detected",
    recommendations = if (length(recommendations) > 0) recommendations else "Primers look good!",
    validation_tools = c(
      "In silico PCR: NCBI e-PCR or ecoPCR",
      "Coverage testing: SILVA TestPrime (https://www.arb-silva.de/search/testprime/)",
      "Specificity: NCBI Primer-BLAST",
      "R package: DECIPHER for primer design"
    )
  )
  
  class(result) <- c("primer_analysis", "list")
  return(result)
}

#' Print Primer Analysis
#'
#' @param x A primer_analysis object
#' @param ... Additional arguments (unused)
#' @export
print.primer_analysis <- function(x, ...) {
  cat("\n🔬 Custom 16S Primer Analysis\n")
  cat("═══════════════════════════════════════════════\n\n")
  
  cat("FORWARD PRIMER:\n")
  cat("Sequence: 5'-", x$forward$sequence, "-3'\n")
  cat("Length: ", x$forward$length, " bp\n")
  cat("GC Content: ", x$forward$gc_content, "%\n")
  cat("Tm: ", x$forward$tm, "°C\n")
  if (x$forward$has_degeneracy) {
    cat("Degenerate bases: ", x$forward$degenerate_count, "\n")
  }
  cat("\n")
  
  cat("REVERSE PRIMER:\n")
  cat("Sequence: 5'-", x$reverse$sequence, "-3'\n")
  cat("Length: ", x$reverse$length, " bp\n")
  cat("GC Content: ", x$reverse$gc_content, "%\n")
  cat("Tm: ", x$reverse$tm, "°C\n")
  if (x$reverse$has_degeneracy) {
    cat("Degenerate bases: ", x$reverse$degenerate_count, "\n")
  }
  cat("\n")
  
  cat("─────────────────────────────────────────────\n")
  cat("Tm Difference: ", x$tm_difference, "°C\n")
  cat("Quality Score: ", x$quality_score, "/100 (", x$quality_rating, ")\n\n")
  
  if (!identical(x$issues, "No issues detected")) {
    cat("⚠️  ISSUES DETECTED:\n")
    for (i in seq_along(x$issues)) {
      cat("  ", i, ". ", x$issues[i], "\n", sep = "")
    }
    cat("\n")
  }
  
  if (!identical(x$recommendations, "Primers look good!")) {
    cat("💡 RECOMMENDATIONS:\n")
    for (i in seq_along(x$recommendations)) {
      cat("  ", i, ". ", x$recommendations[i], "\n", sep = "")
    }
    cat("\n")
  }
  
  cat("🔧 VALIDATION TOOLS:\n")
  for (tool in x$validation_tools) {
    cat("  • ", tool, "\n", sep = "")
  }
  cat("\n")
  
  invisible(x)
}

#' Compare Multiple 16S Primer Pairs
#'
#' Compares multiple 16S primer pairs side-by-side for selection.
#'
#' @param regions Character vector of regions to compare (e.g., c("V3-V4", "V4"))
#' @param criteria Character vector of comparison criteria: 
#'   "coverage", "bias", "size", "all"
#' @return A comparison data.frame
#' @export
#' @examples
#' # Compare V3-V4 and V4 regions
#' comparison <- compare_16s_primers(c("V3-V4", "V4"))
#' print(comparison)
compare_16s_primers <- function(regions = NULL, criteria = "all") {
  
  primers_db <- get_16s_primers()
  
  if (is.null(regions)) {
    # Compare all regions
    comparison <- primers_db
  } else {
    # Filter for specified regions
    comparison <- primers_db[primers_db$region %in% regions, ]
  }
  
  if (nrow(comparison) == 0) {
    stop("No primers found for specified regions")
  }
  
  # Add ranking based on criteria
  if ("coverage" %in% criteria || criteria == "all") {
    coverage_rank <- rank(-match(comparison$coverage, 
                                 c("Very High", "High", "Medium-High", "Medium", "Low")))
    comparison$coverage_rank <- coverage_rank
  }
  
  if ("bias" %in% criteria || criteria == "all") {
    bias_rank <- rank(match(comparison$bias_level, 
                           c("Very Low", "Low", "Medium", "Medium-High", "High")))
    comparison$bias_rank <- bias_rank
  }
  
  # Add summary score
  if (criteria == "all") {
    comparison$overall_score <- round(
      (comparison$coverage_rank + comparison$bias_rank) / 2, 2
    )
    comparison <- comparison[order(comparison$overall_score), ]
  }
  
  message("📊 Primer Pair Comparison")
  message("Top recommendation: ", comparison$region[1], " (", 
          comparison$forward_name[1], "/", comparison$reverse_name[1], ")")
  
  return(comparison)
}

#' Generate Primer Order Sheet
#'
#' Creates a formatted primer order sheet for laboratory use.
#'
#' @param region Character. Target region (e.g., "V4", "V3-V4")
#' @param include_adaptors Logical. Whether to include Illumina adaptor sequences
#' @param concentration Character. Desired primer concentration
#' @param scale Character. Synthesis scale (e.g., "25nm", "100nm")
#' @return A formatted text output for primer ordering
#' @export
#' @examples
#' # Generate order sheet for V4 primers
#' order_sheet <- generate_primer_order(region = "V4")
#' cat(order_sheet)
generate_primer_order <- function(region = "V4",
                                  include_adaptors = TRUE,
                                  concentration = "100 μM",
                                  scale = "25 nmol") {
  
  primers_db <- get_16s_primers()
  primer_info <- primers_db[primers_db$region == region, ]
  
  if (nrow(primer_info) == 0) {
    stop("Region not found in database")
  }
  
  # Illumina adaptors (if requested)
  illumina_forward <- "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"
  illumina_reverse <- "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG"
  
  output <- paste0(
    "\n═══════════════════════════════════════════════\n",
    "  16S rRNA PRIMER ORDER SHEET\n",
    "═══════════════════════════════════════════════\n\n",
    "Target Region: ", primer_info$region, "\n",
    "Application: ", primer_info$application, "\n",
    "Expected Amplicon Size: ", primer_info$amplicon_size, " bp\n\n",
    "─────────────────────────────────────────────\n",
    "FORWARD PRIMER (", primer_info$forward_name, ")\n",
    "─────────────────────────────────────────────\n",
    "Primer Name: ", primer_info$forward_name, "_16S_", primer_info$region, "\n"
  )
  
  if (include_adaptors) {
    output <- paste0(output,
      "Sequence (5'→3'):\n",
      illumina_forward, "\n",
      primer_info$forward_primer, "\n\n",
      "Full Length: ", nchar(illumina_forward) + nchar(primer_info$forward_primer), " bp\n"
    )
  } else {
    output <- paste0(output,
      "Sequence (5'→3'): ", primer_info$forward_primer, "\n",
      "Length: ", nchar(primer_info$forward_primer), " bp\n"
    )
  }
  
  output <- paste0(output,
    "Purification: HPLC or PAGE\n",
    "Scale: ", scale, "\n",
    "Concentration: ", concentration, "\n\n",
    "─────────────────────────────────────────────\n",
    "REVERSE PRIMER (", primer_info$reverse_name, ")\n",
    "─────────────────────────────────────────────\n",
    "Primer Name: ", primer_info$reverse_name, "_16S_", primer_info$region, "\n"
  )
  
  if (include_adaptors) {
    output <- paste0(output,
      "Sequence (5'→3'):\n",
      illumina_reverse, "\n",
      primer_info$reverse_primer, "\n\n",
      "Full Length: ", nchar(illumina_reverse) + nchar(primer_info$reverse_primer), " bp\n"
    )
  } else {
    output <- paste0(output,
      "Sequence (5'→3'): ", primer_info$reverse_primer, "\n",
      "Length: ", nchar(primer_info$reverse_primer), " bp\n"
    )
  }
  
  output <- paste0(output,
    "Purification: HPLC or PAGE\n",
    "Scale: ", scale, "\n",
    "Concentration: ", concentration, "\n\n",
    "─────────────────────────────────────────────\n",
    "NOTES:\n",
    "• Store at -20°C\n",
    "• Make working stocks at 10 μM\n",
    "• Avoid freeze-thaw cycles\n",
    "• For Illumina sequencing, indices will be added separately\n\n",
    "Date: ", Sys.Date(), "\n",
    "═══════════════════════════════════════════════\n"
  )
  
  class(output) <- c("primer_order", "character")
  return(output)
}

#' Print Primer Order Sheet
#'
#' @param x A primer_order object
#' @param ... Additional arguments (unused)
#' @export
print.primer_order <- function(x, ...) {
  cat(x)
  invisible(x)
}