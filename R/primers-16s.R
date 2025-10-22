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


# --- Sanger Sequencing Analysis Functions ---

#' Load and Process Sanger Sequencing Files
#'
#' Loads ABI files, renames them by removing date/time stamps, and creates 
#' groups of forward and reverse reads. Assumes file names are in the 
#' convention: YYYY_MM_DD_16S_ACC#_FWD_WELL_DATE_STAMP.ab1
#'
#' @param path Character. Path to directory containing .ab1 files
#' @return Character vector of unique group names
#' @export
#' @examples
#' \dontrun{
#' groups <- load.files("path/to/abi/files")
#' }
load.files <- function(path = path) {
  
  if (!requireNamespace("stringr", quietly = TRUE)) {
    stop("Package 'stringr' is required. Install with: install.packages('stringr')")
  }
  
  # Get the file names and the list of file names you want to replace them with
  # Assumes file names are in the following convention: YYYY_MM_DD_16S_ACC#_FWD_WELL_DATE_STAMP.ab1
  old_file_names <- dir(path, pattern = ".ab1", full.names = TRUE)
  new_file_names <- gsub("FWD_\\w+\\d+_\\d+_\\d+.*", "FWD.ab1", old_file_names)
  new_file_names <- gsub("REV_\\w+\\d+_\\d+_\\d+.*", "REV.ab1", new_file_names)
  
  file.rename(from = old_file_names, to = new_file_names)
  
  # Index forward and reverse reads
  files_cleaned = new_file_names
  f_matches = stringr::str_match(files_cleaned, "FWD")
  f_indices = which(!is.na(f_matches))
  r_matches = stringr::str_match(files_cleaned, "REV")
  r_indices = which(!is.na(r_matches))
  
  # ONLY keep files that match either FWD or REV suffix
  keep = c(f_indices, r_indices)
  files_cleaned = files_cleaned[keep]
  
  # Remove the suffixes to create a groupname
  files_cleaned = gsub("_FWD.*", "", files_cleaned)
  files_cleaned = gsub("_REV.*", "", files_cleaned)
  
  # Create groups
  group_dataframe = data.frame("file.path" = new_file_names[keep], "group" = files_cleaned)
  groups = unique(group_dataframe$group)
  
  return(groups)
}

#' Generate Sanger Contig from Forward and Reverse Reads
#'
#' Creates a consensus sequence from forward and reverse Sanger reads using
#' the sangeranalyseR package. CB.Contig function generates a contig from a 
#' fwd and reverse read using sangercontig function. FASTA consensus is output 
#' into fasta file. Quality data for fwd and rev read are subsetted and returned 
#' as a summary object.
#'
#' @param path Character. Path to directory containing .ab1 files
#' @param contigName Character. Name for the contig
#' @param suffixForwardRegExp Character. Regular expression for forward read suffix
#' @param suffixReverseRegExp Character. Regular expression for reverse read suffix
#' @param file_name_fwd Character. Forward read filename
#' @param file_name_rev Character. Reverse read filename
#' @return List containing summary statistics and contig object
#' @export
#' @examples
#' \dontrun{
#' contig <- CB.Contig(
#'   path = "path/to/files",
#'   contigName = "Sample01",
#'   suffixForwardRegExp = "_FWD",
#'   suffixReverseRegExp = "_REV",
#'   file_name_fwd = "Sample01_FWD.ab1",
#'   file_name_rev = "Sample01_REV.ab1"
#' )
#' }
CB.Contig <- function(path, contigName, suffixForwardRegExp, 
                                 suffixReverseRegExp, file_name_fwd, file_name_rev) {
  
  if (!requireNamespace("sangeranalyseR", quietly = TRUE)) {
    stop("Package 'sangeranalyseR' is required. Install with: BiocManager::install('sangeranalyseR')")
  }
  
  message("Reading forward and reverse reads and generating contig")
  
  sangerContig <- sangeranalyseR::SangerContig(
    inputSource           = "ABIF",
    ABIF_Directory        = path,
    contigName            = contigName,
    REGEX_SuffixForward   = suffixForwardRegExp,
    REGEX_SuffixReverse   = suffixReverseRegExp,
    TrimmingMethod        = "M2",
    M1TrimmingCutoff      = NULL,
    M2CutoffQualityScore  = 40,
    M2SlidingWindowSize   = 10,
    minReadLength         = 0,
    signalRatioCutoff     = 0.33,
    showTrimmed           = TRUE,
    geneticCode           = sangeranalyseR::GENETIC_CODE
  )
  
  message("Exporting fasta sequence")
  
  # Create output directories if they don't exist
  if (!dir.exists("../Fasta_Sequences/")) {
    dir.create("../Fasta_Sequences/", recursive = TRUE)
  }
  if (!dir.exists("../Results/")) {
    dir.create("../Results/", recursive = TRUE)
  }
  
  sangeranalyseR::writeFasta(sangerContig, outputDir = "../Fasta_Sequences/", selection = "contig")
  
  message("Exporting contig alignment")
  
  alignment <- sangerContig@alignment
  alignment <- Biostrings::DNAMultipleAlignment(alignment)
  filepath <- file.path("../Results", paste0("contigAlign_", contigName, ".txt"))
  DECIPHER::WriteAlignments(alignment, filepath)
  
  message("Subsetting quality data")
  
  QualityFWD <- sangerContig@forwardReadList[[file_name_fwd]]@QualityReport
  QualityREV <- sangerContig@reverseReadList[[file_name_rev]]@QualityReport
  
  message("Generating read summary")
  
  read.summary <- c(
    "consensus.length"       = sangerContig@contigSeq@length,
    "trimmed.seq.length.FWD" = QualityFWD@trimmedSeqLength,
    "trimmed.seq.length.REV" = QualityREV@trimmedSeqLength,
    "trimmed.Mean.qual.FWD"  = QualityFWD@trimmedMeanQualityScore,
    "trimmed.Mean.qual.REV"  = QualityREV@trimmedMeanQualityScore
  )
  
  return(list("summary" = read.summary, "contig" = sangerContig))
}

#' Summarize Sanger Sequencing Data
#'
#' Runs the contig generation function and organizes quality data into a dataframe.
#' Summarize.Sanger function runs sanger contig function and puts subsetted quality 
#' data into a dataframe.
#'
#' @param group Character. Group identifier
#' @param path Character. Path to directory containing .ab1 files
#' @param summarylist List. List to store summary data
#' @return Updated summarylist with quality metrics
#' @export
#' @examples
#' \dontrun{
#' summary <- Summarize.Sanger(
#'   group = "Sample01",
#'   path = "path/to/files",
#'   summarylist = list()
#' )
#' }
Summarize.Sanger <- function(group, path = path, summarylist = summarylist) {
  
  file_name_fwd <- paste0(group, "_FWD.ab1")
  file_name_rev <- paste0(group, "_REV.ab1")
  contigName <- basename(group)
  
  col_names <- c(
    "Consensus length",
    "trim length FWD",
    "trim length REV",
    "trim MeanQual FWD",
    "trim MeanQual REV"
  )
  
  consensus_sequence <- CB.Contig(
    path                = path,
    contigName          = contigName,
    suffixForwardRegExp = "_FWD",
    suffixReverseRegExp = "_REV",
    file_name_fwd       = file_name_fwd,
    file_name_rev       = file_name_rev
  )
  
  # Separate summary data from the consensus sequence object 
  # Turn it into a dataframe with row names from above
  summary <- consensus_sequence$summary
  summary <- t(summary)
  summary <- as.data.frame(summary)
  colnames(summary) <- col_names
  
  # Add a column called sample which is equal to the contig name
  summary$sample <- contigName
  
  summarylist[[group]] <- summary
  
  return(summarylist)
}

#' Analyze Sanger Sequences in Batch
#'
#' Processes all Sanger sequence files in a directory and generates quality reports.
#' Function to run summarize.sanger on groups of files and output summary of quality 
#' results for all samples.
#'
#' @param path Character. Path to directory containing .ab1 files
#' @return Invisibly returns the summary data.frame
#' @export
#' @examples
#' \dontrun{
#' # Analyze all sequences in directory
#' results <- analyze.sequences("path/to/abi/files")
#' }
analyze.sequences <- function(path) {
  
  # Load the files and group into groups based on accession #
  groups <- load.files(path) 
  
  # Generate an empty summary list to put the summary data in
  summarylist <- list()
  
  # Run Summarize.Sanger function on all files in the path to generate fasta files
  summarylist <- lapply(groups, FUN = Summarize.Sanger, path = path, summarylist = summarylist)
  
  # Then concatenate the summary data into a single df
  summary_data <- do.call(rbind, summarylist)
  
  # Change the order of the columns so the contig name is the first column of the df
  summary_data <- summary_data[, c(6, 1:5)]
  
  # Export the summary data into a csv in the Results folder
  resultpath <- file.path("../Results", paste0("Quality_Report_", basename(path), ".csv"))
  write.csv(summary_data, file = resultpath, row.names = FALSE)
  
  message("Analysis complete! Results saved to: ", resultpath)
  
  invisible(summary_data)
}

#' Process Single Sanger Read
#'
#' Analyzes a single Sanger sequencing read (forward or reverse only).
#'
#' @param readFileName Character. Path to the .ab1 file
#' @param readFeature Character. Read feature identifier
#' @return List containing summary statistics and read object
#' @export
#' @examples
#' \dontrun{
#' single_result <- single.read(
#'   readFileName = "path/to/Sample01_FWD.ab1",
#'   readFeature = "Forward"
#' )
#' }
single.read <- function(readFileName, readFeature) {
  
  if (!requireNamespace("sangeranalyseR", quietly = TRUE)) {
    stop("Package 'sangeranalyseR' is required. Install with: BiocManager::install('sangeranalyseR')")
  }
  
  sangerRead <- sangeranalyseR::SangerRead(
    inputSource          = "ABIF",
    readFeature          = readFeature,
    readFileName         = readFileName,
    geneticCode          = sangeranalyseR::GENETIC_CODE,
    TrimmingMethod       = "M2",
    M1TrimmingCutoff     = NULL,
    M2CutoffQualityScore = 15,
    M2SlidingWindowSize  = 10, 
    baseNumPerRow        = 100,
    heightPerRow         = 200,
    signalRatioCutoff    = 0.33,
    showTrimmed          = TRUE
  )
  
  # Create output directory if it doesn't exist
  if (!dir.exists("../Fasta_Sequences/")) {
    dir.create("../Fasta_Sequences/", recursive = TRUE)
  }
  
  sangeranalyseR::writeFasta(sangerRead, outputDir = "../Fasta_Sequences", compress = FALSE)
  
  Quality <- sangerRead@QualityReport
  
  message("Generating read summary")
  
  read.summary <- c(
    "trimmed.seq.length" = Quality@trimmedSeqLength,
    "trimmed.Mean.qual"  = Quality@trimmedMeanQualityScore
  )
  
  return(list("summary" = read.summary, "Read" = sangerRead))
}

#' Summarize Single Sanger Read
#'
#' Creates a summary dataframe for a single Sanger read. Summarize.Single function 
#' runs sanger contig function and puts subsetted quality data into a dataframe.
#'
#' @param readFileName Character. Path to the .ab1 file
#' @param readFeature Character. Read feature identifier
#' @param summarylist List. List to store summary data
#' @return Updated summarylist with quality metrics
#' @export
#' @examples
#' \dontrun{
#' summary <- Summarize.Single(
#'   readFileName = "path/to/Sample01_FWD.ab1",
#'   readFeature = "Forward",
#'   summarylist = list()
#' )
#' }
Summarize.Single <- function(readFileName, readFeature, summarylist = summarylist) {
  
  singleName = basename(readFileName)
  
  col_names = c("trim length", "trim MeanQual")
  
  single_sequence <- single.read(
    readFileName = readFileName,
    readFeature  = readFeature
  )
  
  # Separate summary data from the read sequence object 
  # Turn it into a dataframe with row names from above
  summary <- single_sequence$summary
  summary <- t(summary)
  summary <- as.data.frame(summary)
  colnames(summary) <- col_names
  
  # Add a column called sample which is equal to the read name
  summary$sample <- singleName
  
  summarylist[[readFileName]] <- summary
  
  return(summarylist)
}

#' Analyze Single Sanger Sequence
#'
#' Processes a single Sanger read and generates a quality report. Function to run 
#' summarize.sanger on single reads and output summary of quality results for all samples.
#'
#' @param readFileName Character. Path to the .ab1 file
#' @param readFeature Character. Read feature identifier
#' @return Invisibly returns the summary data
#' @export
#' @examples
#' \dontrun{
#' # Analyze single forward read
#' result <- analyze.single.sequence(
#'   readFileName = "path/to/Sample01_FWD.ab1",
#'   readFeature = "Forward"
#' )
#' }
analyze.single.sequence <- function(readFileName, readFeature) {
  
  # Generate an empty summary list to put the summary data in
  summarylist = list()
  
  # Run Summarize.Single function on all files in the path to generate fasta files
  summarylist <- Summarize.Single(
    readFileName = readFileName, 
    readFeature  = readFeature, 
    summarylist  = summarylist
  )
  
  # Extract the summary data
  summary_data <- summarylist[[1]]
  
  # Change the order of the columns so the read name is the first column
  summary_data <- summary_data[, c(3, 1, 2)]
  
  # Create output directory if it doesn't exist
  if (!dir.exists("../Results/")) {
    dir.create("../Results/", recursive = TRUE)
  }
  
  # Export the summary data into a csv in the Results folder
  resultpath <- file.path("../Results", paste0("Quality_Report_", basename(readFileName), ".csv"))
  write.csv(summary_data, file = resultpath, row.names = FALSE)
  
  message("Analysis complete! Results saved to: ", resultpath)
  
  invisible(summary_data)
}


# --- BLAST Functions for 16S and ITS Database Searches ---

#' Edit FASTA File to Remove Spaces
#'
#' Removes spaces from FASTA file headers and replaces them with underscores.
#' This is necessary for proper BLAST processing.
#'
#' @param x Character. Path to FASTA file
#' @return NULL. File is modified in place
#' @export
#' @examples
#' \dontrun{
#' edit.fasta("path/to/sequences.fasta")
#' }
edit.fasta <- function(x) {
  
  if (!requireNamespace("seqinr", quietly = TRUE)) {
    stop("Package 'seqinr' is required. Install with: install.packages('seqinr')")
  }
  
  fasta <- seqinr::read.fasta(x, whole.header = TRUE, as.string = TRUE)
  noSpace <- lapply(names(fasta), function(y) gsub('\\s+', "_", y))
  seqinr::write.fasta(sequences = fasta, names = noSpace, file.out = x)
}

#' Run Local BLAST Search
#'
#' Executes a local BLAST search using blastn command and returns results as a tibble.
#' Requires NCBI BLAST+ to be installed and accessible in the system PATH.
#'
#' @param x Character. Path to query FASTA file
#' @param blastn Character. Path to blastn executable (default: "blastn")
#' @param blast_db Character. Path to BLAST database
#' @param input Character. Input file path (default: x)
#' @param evalue Numeric. E-value threshold (default: 0.05)
#' @param format Character. Output format specification
#' @param max_target_seqs Integer. Maximum number of aligned sequences to keep (default: 100)
#' @param gapopen Integer. Cost to open a gap (default: 0)
#' @param word_size Integer. Word size for wordfinder algorithm (default: 28)
#' @param qcov_hsp_perc Numeric. Minimum query coverage per HSP (default: 90)
#' @return Tibble containing BLAST results
#' @export
#' @examples
#' \dontrun{
#' blast_results <- Blast.CB(
#'   x = "sequences.fasta",
#'   blast_db = "path/to/16S_ribosomal_RNA"
#' )
#' }
Blast.CB <- function(x, 
                     blastn = "blastn",
                     blast_db = blast_db,
                     input = x,
                     evalue = 0.05,
                     format = '"6 sscinames pident length qcovhsp mismatch gapopen qstart qend sstart send evalue bitscore"',
                     max_target_seqs = 100,
                     gapopen = 0,
                     word_size = 28,
                     qcov_hsp_perc = 90) {
  
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required. Install with: install.packages('dplyr')")
  }
  
  if (!requireNamespace("tidyr", quietly = TRUE)) {
    stop("Package 'tidyr' is required. Install with: install.packages('tidyr')")
  }
  
  colnames <- c("Name",
                "pident",
                "length",
                "query_coverage",
                "mismatch",
                "gapopen",
                "qstart",
                "qend",
                "sstart",
                "send",
                "evalue",
                "bitscore")
  
  blast_out <- system2(
    command = blastn, 
    args = c("-db", blast_db, 
             "-query", input, 
             "-outfmt", format, 
             "-evalue", evalue, 
             "-gapopen", gapopen,
             "-max_target_seqs", max_target_seqs,
             "-word_size", word_size,
             "-qcov_hsp_perc", qcov_hsp_perc),
    wait = TRUE,
    stdout = TRUE
  )
  
  blast_out <- dplyr::as_tibble(blast_out)
  blast_out <- tidyr::separate(blast_out, col = value, into = colnames, sep = "\t", convert = TRUE)
  
  return(blast_out)
}

#' BLAST Single FASTA File
#'
#' Processes a single FASTA file: removes spaces from headers, runs BLAST,
#' and saves results to CSV.
#'
#' @param file_name Character. Path to FASTA file
#' @param blast_db Character. Path to BLAST database
#' @param DBname Character. Database name for output file naming
#' @return NULL. Results are saved to CSV file
#' @export
#' @examples
#' \dontrun{
#' Blast.all(
#'   file_name = "sample_16S.fasta",
#'   blast_db = "path/to/16S_ribosomal_RNA",
#'   DBname = "16S"
#' )
#' }
Blast.all <- function(file_name, blast_db, DBname) {
  
  name = basename(file_name)
  DB = DBname
  
  edit.fasta(file_name)
  Blast_output <- Blast.CB(file_name, blast_db = blast_db) 
  
  # Create output directory if it doesn't exist
  if (!dir.exists("../Results/")) {
    dir.create("../Results/", recursive = TRUE)
  }
  
  mypath <- file.path("../Results/", paste0("Result_", DB, "_", name, ".csv"))
  write.csv(Blast_output, file = mypath, row.names = FALSE)
  
  message("BLAST results saved to: ", mypath)
}

#' BLAST Multiple Files Against 16S and ITS Databases
#'
#' Processes all FASTA files in a directory, automatically detecting whether
#' they contain 16S or ITS sequences based on filename, and BLASTs against
#' the appropriate database.
#'
#' @param Blastpath Character. Path to directory containing FASTA files
#' @param blast16Sdb Character. Path to 16S ribosomal RNA BLAST database
#' @param blastITSdb Character. Path to ITS RefSeq Fungi BLAST database
#' @param DBname Character. Database name for output file naming
#' @return NULL. Results are saved to CSV files
#' @export
#' @examples
#' \dontrun{
#' # Set BLAST environment (do this once per session)
#' Sys.setenv(PATH = paste(Sys.getenv("PATH"), 
#'                         "path/to/ncbi-blast/bin", 
#'                         sep = .Platform$path.sep))
#' Sys.setenv(BLASTDB = "path/to/ncbi-blast/db")
#' 
#' # BLAST all files
#' Blast.Files(
#'   Blastpath = "path/to/fasta/files",
#'   blast16Sdb = "path/to/db/16S_ribosomal_RNA",
#'   blastITSdb = "path/to/db/ITS_RefSeq_Fungi",
#'   DBname = "MyProject"
#' )
#' }
Blast.Files <- function(Blastpath, 
                        blast16Sdb = "../ncbi-blast-2.13.0+/db/16S_ribosomal_RNA", 
                        blastITSdb = "../ncbi-blast-2.13.0+/db/ITS_RefSeq_Fungi", 
                        DBname = DBname) {
  
  if (!requireNamespace("stringr", quietly = TRUE)) {
    stop("Package 'stringr' is required. Install with: install.packages('stringr')")
  }
  
  file_names <- dir(Blastpath, pattern = ".fa|.fasta", full.names = TRUE)
  
  message("Generating 16S sequence list")
  
  matches_16S = stringr::str_match(file_names, "16S")
  indices_16S = which(!is.na(matches_16S))
  file_names_16S = file_names[indices_16S]
  
  if (length(file_names_16S) > 0) {
    message("BLASTing ", length(file_names_16S), " 16S sequences")
    lapply(file_names_16S, FUN = Blast.all, blast_db = blast16Sdb, DBname = DBname)
  } else {
    message("No 16S sequences found")
  }
  
  message("Generating ITS sequence list")
  
  matches_ITS = stringr::str_match(file_names, "ITS")
  indices_ITS = which(!is.na(matches_ITS))
  file_names_ITS = file_names[indices_ITS]
  
  if (length(file_names_ITS) > 0) {
    message("BLASTing ", length(file_names_ITS), " ITS sequences")
    lapply(file_names_ITS, FUN = Blast.all, blast_db = blastITSdb, DBname = DBname)
  } else {
    message("No ITS sequences found")
  }
  
  message("BLAST analysis complete!")
}

#' Set BLAST Environment Variables
#'
#' Helper function to set up BLAST+ environment variables for the current R session.
#' Call this before using any BLAST functions.
#'
#' @param blast_path Character. Path to NCBI BLAST+ installation directory
#' @return NULL. Environment variables are set
#' @export
#' @examples
#' \dontrun{
#' setup_blast_env("path/to/ncbi-blast-2.13.0+")
#' }
setup_blast_env <- function(blast_path = "../ncbi-blast-2.13.0+") {
  
  bin_path <- file.path(blast_path, "bin")
  db_path <- file.path(blast_path, "db")
  
  Sys.setenv(PATH = paste(Sys.getenv("PATH"), bin_path, sep = .Platform$path.sep))
  Sys.setenv(PATH = paste(Sys.getenv("PATH"), db_path, sep = .Platform$path.sep))
  Sys.setenv(BLASTDB = db_path)
  
  message("BLAST environment configured:")
  message("  Binary path: ", bin_path)
  message("  Database path: ", db_path)
  message("  BLASTDB: ", Sys.getenv("BLASTDB"))
}