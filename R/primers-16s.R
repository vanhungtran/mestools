# 16S rRNA Primer Functions for Microbiome Sequencing
# Author: mestools package
# Date: 2025-10-22

# --- 16S rRNA Primer Database and Functions ---

#' Get Standard 16S rRNA Primer Pairs
#'
#' Returns a comprehensive database of commonly used 16S rRNA primer pairs
#' for microbiome sequencing studies. Includes 20+ validated primer sets
#' from major publications and sequencing projects.
#'
#' @return A data.frame containing primer information with columns:
#'   \itemize{
#'     \item region: Hypervariable region(s) targeted
#'     \item forward_name: Forward primer identifier
#'     \item forward_primer: Forward primer sequence (5' to 3')
#'     \item reverse_name: Reverse primer identifier
#'     \item reverse_primer: Reverse primer sequence (5' to 3')
#'     \item application: Recommended use case
#'     \item amplicon_size: Expected amplicon length (bp)
#'     \item coverage: Taxonomic coverage rating
#'     \item bias_level: PCR bias assessment
#'     \item reference: Original publication
#'     \item platform_optimized: Recommended sequencing platform
#'     \item notes: Additional usage information
#'   }
#' @export
#' @examples
#' # Get all standard primer pairs
#' primers <- get_16s_primers()
#' print(primers)
#' 
#' # Filter for V4 region
#' v4_primers <- primers[primers$region == "V4", ]
#' 
#' # Find primers with very low bias
#' low_bias <- primers[primers$bias_level == "Very Low", ]
#' 
#' # Get Earth Microbiome Project primers
#' emp_primers <- primers[grepl("EMP", primers$application), ]
#' 
#' # View all available references
#' unique(primers$reference)
#' 
#' @references
#' \itemize{
#'   \item Lane (1991) 16S/23S rRNA sequencing. In Nucleic Acid Techniques.
#'   \item Caporaso et al. (2011) Global patterns of 16S rRNA diversity. PNAS.
#'   \item Klindworth et al. (2013) Evaluation of primers for amplicon sequencing. NAR.
#'   \item Apprill et al. (2015) Minor revision to V4 primers. Aquat Microb Ecol.
#'   \item Parada et al. (2016) Every base matters. Environ Microbiol.
#'   \item Walters et al. (2016) Improved primers for Illumina sequencing. ISME J.
#' }
get_16s_primers <- function() {
  primers_df <- data.frame(
    region = c(
      # Original primers
      "V1-V2", "V1-V3", "V3-V4", "V4", "V4-V5", "V6-V8", "V7-V9",
      # Additional widely-used primers
      "V1-V3", "V3-V4", "V4", "V4-V5", "V3-V5", "V5-V6", "V5-V7",
      "V6-V8", "V1-V9", "V2-V3", "V4-V6", "V3-V4", "V4"
    ),
    forward_name = c(
      # Original
      "27F", "27F", "341F", "515F", "515F", "939F", "1115F",
      # Additional
      "8F", "357F", "518F", "563F", "515F-Y", "784F", "799F",
      "926F", "27F-Full", "104F", "515F-Parada", "347F", "519F"
    ),
    forward_primer = c(
      # Original
      "AGAGTTTGATCMTGGCTCAG",
      "AGAGTTTGATCMTGGCTCAG",
      "CCTACGGGNGGCWGCAG",
      "GTGCCAGCMGCCGCGGTAA",
      "GTGCCAGCMGCCGCGGTAA",
      "TTGACGGGGGCCCGCAC",
      "GCAACGAGCGCAACCC",
      # Additional
      "AGAGTTTGATCCTGGCTCAG",           # 8F - Bacterial universal
      "CCTACGGGAGGCAGCAG",               # 357F - Common alternative
      "CCAGCAGCYGCGGTAAN",               # 518F - Alternative V4
      "GCCAGCMGCCGCGGTAA",               # 563F - Modified 515F
      "GTGYCAGCMGCCGCGGTAA",             # 515F-Y - Apprill modification
      "AGGATTAGATACCCTGGTA",             # 784F - V5-V6
      "AACMGGATTAGATACCCKG",             # 799F - V5-V7
      "AAACTYAAAKGAATTGACGG",            # 926F - V6-V8 alternative
      "AGAGTTTGATYMTGGCTCAG",            # 27F-Full - Universal with Y
      "GGCGVACGGGTGAGTAA",               # 104F - V2-V3
      "GTGYCAGCMGCCGCGGTAA",             # 515F-Parada - Improved V4
      "GGAGGCAGCAGTRRGGAAT",             # 347F - Alternative V3-V4
      "CAGCMGCCGCGGTAATWC"               # 519F - EMP modified
    ),
    reverse_name = c(
      # Original
      "338R", "534R", "785R", "806R", "944R", "1378R", "1492R",
      # Additional
      "519R", "926R", "786R", "926R", "806R-Y", "1061R", "1193R",
      "1392R", "1492R-Full", "338R", "907R", "803R", "785R"
    ),
    reverse_primer = c(
      # Original
      "TGCTGCCTCCCGTAGGAGT",
      "ATTACCGCGGCTGCTGG",
      "GACTACHVGGGTATCTAATCC",
      "GGACTACHVGGGTWTCTAAT",
      "GAATTGGCACCTTGAGGAC",
      "CGGTGTGTACAAGGCCCGGGAACG",
      "GGTTACCTTGTTACGACTT",
      # Additional
      "GWATTACCGCGGCKGCTG",             # 519R - Broad coverage
      "CCGYCAATTYMTTTRAGTTT",           # 926R - V3-V4 alternative
      "GGACTACCAGGGTATCTAAT",           # 786R - Alternative V4
      "CCGTCAATTCMTTTGAGTTT",           # 926R - V4-V5
      "GGACTACNVGGGTWTCTAAT",           # 806R-Y - Apprill modification
      "CRRCACGAGCTGACGAC",               # 1061R - V5-V6
      "ACGTCATCCCCACCTTCC",              # 1193R - V5-V7
      "ACGGGCGGTGTGTRC",                 # 1392R - V6-V8
      "RGYTACCTTGTTACGACTT",             # 1492R-Full - Full length
      "ACTCCTACGGGAGGCAGC",              # 338R - V2-V3
      "CCGTCAATTCMTTTRAGTTT",            # 907R - V4-V6
      "CTACCRGGGTATCTAATCC",             # 803R - V3-V4 alt
      "TACNVGGGTATCTAATCC"               # 785R - Alternative
    ),
    application = c(
      # Original
      "Broad bacterial profiling",
      "Environmental microbiomes",
      "Gut microbiota (widely used)",
      "Universal (Earth Microbiome Project)",
      "Marine/environmental studies",
      "Archaeal and bacterial detection",
      "Full-length near-universal",
      # Additional
      "Bacterial universal primer",
      "Human microbiome, high specificity",
      "Alternative Earth Microbiome Project",
      "Modified EMP for marine samples",
      "Cyanobacteria/chloroplast improved",
      "Deep sequencing applications",
      "Soil/environmental diversity",
      "Archaeal diversity studies",
      "Full-length 16S for long reads",
      "Human-associated microbiomes",
      "Complex environmental samples",
      "Klindworth optimized primers",
      "Modified EMP with improved coverage"
    ),
    amplicon_size = c(
      # Original
      311, 507, 444, 291, 429, 439, 377,
      # Additional
      511, 569, 268, 363, 291, 277, 394,
      466, 1465, 234, 392, 456, 438
    ),
    coverage = c(
      # Original
      "High", "High", "Very High", "Very High", "High", "Medium-High", "High",
      # Additional
      "Very High", "High", "Very High", "High", "Very High", "Medium", "High",
      "Medium-High", "Very High", "High", "High", "High", "Very High"
    ),
    bias_level = c(
      # Original
      "Low", "Low", "Very Low", "Very Low", "Low", "Medium", "Low",
      # Additional
      "Low", "Very Low", "Very Low", "Low", "Very Low", "Medium", "Low",
      "Medium", "Low", "Low", "Low", "Low", "Very Low"
    ),
    reference = c(
      # Original
      "Lane 1991",
      "Lane 1991 / Muyzer 1993",
      "Klindworth 2013",
      "Caporaso 2011 (EMP)",
      "Parada 2016",
      "PMHz 2007",
      "Turner 1999",
      # Additional
      "Edwards 1989",
      "Nadkarni 2002",
      "Caporaso 2012",
      "Parada 2016",
      "Apprill 2015",
      "Anderson 2008",
      "Chelius 2001",
      "DeLong 1992",
      "Weisburg 1991",
      "Sundquist 2007",
      "Tamaki 2011",
      "Klindworth 2013",
      "Walters 2016"
    ),
    platform_optimized = c(
      # Original
      "Illumina MiSeq", "Illumina MiSeq", "Illumina MiSeq/NovaSeq", 
      "Illumina MiSeq/NovaSeq", "Illumina MiSeq", "454/Illumina", "PacBio/Nanopore",
      # Additional
      "Sanger/PacBio", "Illumina MiSeq", "Illumina MiSeq/NovaSeq",
      "Illumina MiSeq", "Illumina MiSeq/NovaSeq", "Illumina MiSeq",
      "Illumina HiSeq", "454/Ion Torrent", "PacBio/Nanopore",
      "Illumina MiSeq", "Illumina MiSeq", "Illumina MiSeq/NovaSeq",
      "Illumina MiSeq/NovaSeq"
    ),
    notes = c(
      # Original
      "Classical bacterial primer, broad coverage",
      "Good for diverse environmental samples",
      "Most widely used, excellent for gut microbiome",
      "Earth Microbiome Project standard",
      "Improved coverage for marine samples",
      "Includes archaea, longer amplicon",
      "Nearly complete 16S, requires long reads",
      # Additional
      "Original bacterial universal primer",
      "High specificity, low bias",
      "Alternative to 515F/806R",
      "Better for cyanobacteria-rich samples",
      "Reduced bias against SAR11 clade and Thaumarchaeota",
      "Good taxonomic resolution, moderate length",
      "Excellent for soil microbial diversity",
      "Archaeal primer, use with caution for bacteria",
      "Best for species-level identification",
      "Short amplicon, good for degraded DNA",
      "Longer amplicon, good taxonomic resolution",
      "Optimized for maximum coverage",
      "Improved version of EMP primers"
    ),
    stringsAsFactors = FALSE
  )
  
  return(primers_df)
}

#' Search and Filter 16S Primer Database
#'
#' Advanced search function to filter primers by multiple criteria including
#' region, amplicon size, platform, coverage, and reference.
#'
#' @param region Character. Filter by hypervariable region (e.g., "V4", "V3-V4")
#' @param amplicon_size_min Integer. Minimum amplicon size (bp)
#' @param amplicon_size_max Integer. Maximum amplicon size (bp)
#' @param platform Character. Filter by sequencing platform
#' @param coverage Character. Filter by coverage rating ("High", "Very High", etc.)
#' @param bias_level Character. Filter by bias level ("Low", "Very Low", etc.)
#' @param reference Character. Filter by publication/reference
#' @param keyword Character. Search in application or notes fields
#' @return Filtered data.frame of primers matching criteria
#' @export
#' @examples
#' \dontrun{
#' # Find all V4 primers
#' v4_primers <- search_16s_primers(region = "V4")
#' 
#' # Find short amplicons for degraded DNA
#' short_primers <- search_16s_primers(amplicon_size_max = 300)
#' 
#' # Find primers optimized for Illumina with very high coverage
#' illumina_primers <- search_16s_primers(
#'   platform = "Illumina", 
#'   coverage = "Very High"
#' )
#' 
#' # Find primers for marine studies
#' marine_primers <- search_16s_primers(keyword = "marine")
#' 
#' # Find Earth Microbiome Project primers
#' emp_primers <- search_16s_primers(reference = "Caporaso")
#' 
#' # Complex search: V3-V4 region, low bias, for Illumina
#' optimal <- search_16s_primers(
#'   region = "V3-V4",
#'   bias_level = "Very Low",
#'   platform = "Illumina"
#' )
#' }
search_16s_primers <- function(region = NULL,
                               amplicon_size_min = NULL,
                               amplicon_size_max = NULL,
                               platform = NULL,
                               coverage = NULL,
                               bias_level = NULL,
                               reference = NULL,
                               keyword = NULL) {
  
  primers <- get_16s_primers()
  
  # Apply filters
  if (!is.null(region)) {
    primers <- primers[primers$region == region, ]
  }
  
  if (!is.null(amplicon_size_min)) {
    primers <- primers[primers$amplicon_size >= amplicon_size_min, ]
  }
  
  if (!is.null(amplicon_size_max)) {
    primers <- primers[primers$amplicon_size <= amplicon_size_max, ]
  }
  
  if (!is.null(platform)) {
    primers <- primers[grepl(platform, primers$platform_optimized, ignore.case = TRUE), ]
  }
  
  if (!is.null(coverage)) {
    primers <- primers[primers$coverage == coverage, ]
  }
  
  if (!is.null(bias_level)) {
    primers <- primers[primers$bias_level == bias_level, ]
  }
  
  if (!is.null(reference)) {
    primers <- primers[grepl(reference, primers$reference, ignore.case = TRUE), ]
  }
  
  if (!is.null(keyword)) {
    keyword_match <- grepl(keyword, primers$application, ignore.case = TRUE) |
                     grepl(keyword, primers$notes, ignore.case = TRUE)
    primers <- primers[keyword_match, ]
  }
  
  if (nrow(primers) == 0) {
    message("No primers found matching the specified criteria")
    return(data.frame())
  }
  
  message("Found ", nrow(primers), " primer pair(s) matching criteria")
  
  return(primers)
}

#' Get Primer Statistics and Summary
#'
#' Provides summary statistics about the 16S primer database including
#' coverage by region, platform distribution, and reference counts.
#'
#' @return A list containing summary statistics
#' @export
#' @examples
#' \dontrun{
#' stats <- get_primer_stats()
#' print(stats)
#' }
get_primer_stats <- function() {
  
  primers <- get_16s_primers()
  
  stats <- list(
    total_primers = nrow(primers),
    
    by_region = table(primers$region),
    
    by_platform = table(primers$platform_optimized),
    
    by_coverage = table(primers$coverage),
    
    by_bias = table(primers$bias_level),
    
    amplicon_range = range(primers$amplicon_size),
    
    mean_amplicon = mean(primers$amplicon_size),
    
    references = unique(primers$reference),
    
    regions_available = unique(primers$region),
    
    platforms_supported = unique(primers$platform_optimized)
  )
  
  return(stats)
}

#' Print Primer Statistics
#'
#' Pretty-print summary of primer database statistics
#'
#' @export
#' @examples
#' \dontrun{
#' print_primer_stats()
#' }
print_primer_stats <- function() {
  
  stats <- get_primer_stats()
  
  cat("\n=== 16S Primer Database Statistics ===\n\n")
  
  cat("Total primer pairs:", stats$total_primers, "\n\n")
  
  cat("Coverage by Region:\n")
  print(stats$by_region)
  
  cat("\nCoverage Rating Distribution:\n")
  print(stats$by_coverage)
  
  cat("\nBias Level Distribution:\n")
  print(stats$by_bias)
  
  cat("\nAmplicon Size Range:", stats$amplicon_range[1], "-", stats$amplicon_range[2], "bp\n")
  cat("Mean Amplicon Size:", round(stats$mean_amplicon, 0), "bp\n\n")
  
  cat("Regions Available:\n")
  cat(paste(" -", stats$regions_available), sep = "\n")
  
  cat("\nPlatforms Supported:\n")
  cat(paste(" -", stats$platforms_supported), sep = "\n")
  
  cat("\nReferences Included:\n")
  cat(paste(" -", stats$references), sep = "\n")
  
  cat("\n")
  
  invisible(stats)
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
  
  # Input validation
  if (missing(path) || is.null(path) || path == "") {
    stop("Please provide a valid path to the directory containing .ab1 files")
  }
  
  if (!dir.exists(path)) {
    stop("Directory does not exist: ", path)
  }
  
  if (!requireNamespace("stringr", quietly = TRUE)) {
    stop("Package 'stringr' is required. Install with: install.packages('stringr')")
  }
  
  # Get the file names and the list of file names you want to replace them with
  # Assumes file names are in the following convention: YYYY_MM_DD_16S_ACC#_FWD_WELL_DATE_STAMP.ab1
  old_file_names <- dir(path, pattern = ".ab1", full.names = TRUE)
  
  if (length(old_file_names) == 0) {
    stop("No .ab1 files found in directory: ", path)
  }
  
  message("Found ", length(old_file_names), " .ab1 files")
  
  new_file_names <- gsub("FWD_\\w+\\d+_\\d+_\\d+.*", "FWD.ab1", old_file_names)
  new_file_names <- gsub("REV_\\w+\\d+_\\d+_\\d+.*", "REV.ab1", new_file_names)
  
  # Only rename if names have changed
  files_to_rename <- old_file_names != new_file_names
  if (any(files_to_rename)) {
    rename_result <- file.rename(from = old_file_names[files_to_rename], 
                                  to = new_file_names[files_to_rename])
    if (!all(rename_result)) {
      warning("Some files could not be renamed")
    }
  }
  
  # Index forward and reverse reads
  files_cleaned <- new_file_names
  f_matches <- stringr::str_match(files_cleaned, "FWD")
  f_indices <- which(!is.na(f_matches))
  r_matches <- stringr::str_match(files_cleaned, "REV")
  r_indices <- which(!is.na(r_matches))
  
  # ONLY keep files that match either FWD or REV suffix
  keep <- c(f_indices, r_indices)
  
  if (length(keep) == 0) {
    stop("No files with 'FWD' or 'REV' suffix found in directory: ", path)
  }
  
  files_cleaned <- files_cleaned[keep]
  
  # Remove the suffixes to create a groupname
  files_cleaned <- gsub("_FWD.*", "", files_cleaned)
  files_cleaned <- gsub("_REV.*", "", files_cleaned)
  
  # Create groups
  group_dataframe <- data.frame("file.path" = new_file_names[keep], 
                                "group" = files_cleaned,
                                stringsAsFactors = FALSE)
  groups <- unique(group_dataframe$group)
  
  message("Created ", length(groups), " groups from ", length(keep), " files")
  
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
  
  # Input validation
  if (missing(path) || !dir.exists(path)) {
    stop("Please provide a valid path to the directory containing .ab1 files")
  }
  
  if (missing(contigName) || is.null(contigName) || contigName == "") {
    stop("Please provide a valid contig name")
  }
  
  if (!requireNamespace("sangeranalyseR", quietly = TRUE)) {
    stop("Package 'sangeranalyseR' is required. Install with: BiocManager::install('sangeranalyseR')")
  }
  
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop("Package 'Biostrings' is required. Install with: BiocManager::install('Biostrings')")
  }
  
  if (!requireNamespace("DECIPHER", quietly = TRUE)) {
    stop("Package 'DECIPHER' is required. Install with: BiocManager::install('DECIPHER')")
  }
  
  message("Reading forward and reverse reads and generating contig for: ", contigName)
  
  # Wrap in tryCatch for better error handling
  sangerContig <- tryCatch({
    sangeranalyseR::SangerContig(
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
  }, error = function(e) {
    stop("Failed to create contig for ", contigName, ": ", conditionMessage(e))
  })
  
  message("Exporting fasta sequence")
  
  # Create output directories if they don't exist
  fasta_dir <- "../Fasta_Sequences/"
  results_dir <- "../Results/"
  
  if (!dir.exists(fasta_dir)) {
    dir.create(fasta_dir, recursive = TRUE)
    message("Created directory: ", fasta_dir)
  }
  if (!dir.exists(results_dir)) {
    dir.create(results_dir, recursive = TRUE)
    message("Created directory: ", results_dir)
  }
  
  tryCatch({
    sangeranalyseR::writeFasta(sangerContig, outputDir = fasta_dir, selection = "contig")
    message("FASTA file saved to: ", fasta_dir)
  }, error = function(e) {
    warning("Failed to write FASTA file: ", conditionMessage(e))
  })
  
  message("Exporting contig alignment")
  
  tryCatch({
    alignment <- sangerContig@alignment
    alignment <- Biostrings::DNAMultipleAlignment(alignment)
    filepath <- file.path(results_dir, paste0("contigAlign_", contigName, ".txt"))
    DECIPHER::WriteAlignments(alignment, filepath)
    message("Alignment saved to: ", filepath)
  }, error = function(e) {
    warning("Failed to write alignment file: ", conditionMessage(e))
  })
  
  message("Extracting quality data")
  
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
  
  message("Contig processing complete for: ", contigName)
  
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
  
  # Input validation
  if (missing(group) || is.null(group)) {
    stop("Please provide a valid group identifier")
  }
  
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
  
  message("Processing group: ", contigName)
  
  consensus_sequence <- tryCatch({
    CB.Contig(
      path                = path,
      contigName          = contigName,
      suffixForwardRegExp = "_FWD",
      suffixReverseRegExp = "_REV",
      file_name_fwd       = file_name_fwd,
      file_name_rev       = file_name_rev
    )
  }, error = function(e) {
    warning("Failed to process group ", contigName, ": ", conditionMessage(e))
    return(NULL)
  })
  
  if (is.null(consensus_sequence)) {
    return(summarylist)
  }
  
  # Separate summary data from the consensus sequence object 
  # Turn it into a dataframe with row names from above
  summary <- consensus_sequence$summary
  summary <- t(summary)
  summary <- as.data.frame(summary, stringsAsFactors = FALSE)
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
  
  # Input validation
  if (missing(path) || !dir.exists(path)) {
    stop("Please provide a valid path to the directory containing .ab1 files")
  }
  
  message("Starting batch analysis of Sanger sequences...")
  
  # Load the files and group into groups based on accession #
  groups <- load.files(path) 
  
  if (length(groups) == 0) {
    stop("No valid groups found in directory: ", path)
  }
  
  # Generate an empty summary list to put the summary data in
  summarylist <- list()
  
  # Run Summarize.Sanger function on all files in the path to generate fasta files
  message("Processing ", length(groups), " groups...")
  summarylist <- lapply(groups, FUN = Summarize.Sanger, path = path, summarylist = summarylist)
  
  # Remove NULL entries (failed processing)
  summarylist <- Filter(Negate(is.null), summarylist)
  
  if (length(summarylist) == 0) {
    stop("No sequences were successfully processed")
  }
  
  # Then concatenate the summary data into a single df
  summary_data <- do.call(rbind, summarylist)
  
  # Change the order of the columns so the contig name is the first column of the df
  summary_data <- summary_data[, c(6, 1:5)]
  
  # Export the summary data into a csv in the Results folder
  results_dir <- "../Results/"
  if (!dir.exists(results_dir)) {
    dir.create(results_dir, recursive = TRUE)
  }
  
  resultpath <- file.path(results_dir, paste0("Quality_Report_", basename(path), ".csv"))
  write.csv(summary_data, file = resultpath, row.names = FALSE)
  
  message("Analysis complete! Results saved to: ", resultpath)
  message("Successfully processed ", nrow(summary_data), " sequences")
  
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
  
  # Input validation
  if (missing(readFileName) || !file.exists(readFileName)) {
    stop("Please provide a valid path to an .ab1 file")
  }
  
  if (missing(readFeature) || is.null(readFeature)) {
    stop("Please provide a valid read feature (e.g., 'Forward' or 'Reverse')")
  }
  
  if (!requireNamespace("sangeranalyseR", quietly = TRUE)) {
    stop("Package 'sangeranalyseR' is required. Install with: BiocManager::install('sangeranalyseR')")
  }
  
  message("Processing single read: ", basename(readFileName))
  
  sangerRead <- tryCatch({
    sangeranalyseR::SangerRead(
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
  }, error = function(e) {
    stop("Failed to process read ", basename(readFileName), ": ", conditionMessage(e))
  })
  
  # Create output directory if it doesn't exist
  fasta_dir <- "../Fasta_Sequences/"
  if (!dir.exists(fasta_dir)) {
    dir.create(fasta_dir, recursive = TRUE)
    message("Created directory: ", fasta_dir)
  }
  
  tryCatch({
    sangeranalyseR::writeFasta(sangerRead, outputDir = fasta_dir, compress = FALSE)
    message("FASTA file saved to: ", fasta_dir)
  }, error = function(e) {
    warning("Failed to write FASTA file: ", conditionMessage(e))
  })
  
  Quality <- sangerRead@QualityReport
  
  message("Generating read summary")
  
  read.summary <- c(
    "trimmed.seq.length" = Quality@trimmedSeqLength,
    "trimmed.Mean.qual"  = Quality@trimmedMeanQualityScore
  )
  
  message("Read processing complete")
  
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
  
  # Input validation
  if (missing(readFileName) || !file.exists(readFileName)) {
    stop("Please provide a valid path to an .ab1 file")
  }
  
  singleName <- basename(readFileName)
  
  col_names <- c("trim length", "trim MeanQual")
  
  message("Summarizing single read: ", singleName)
  
  single_sequence <- tryCatch({
    single.read(
      readFileName = readFileName,
      readFeature  = readFeature
    )
  }, error = function(e) {
    warning("Failed to process ", singleName, ": ", conditionMessage(e))
    return(NULL)
  })
  
  if (is.null(single_sequence)) {
    return(summarylist)
  }
  
  # Separate summary data from the read sequence object 
  # Turn it into a dataframe with row names from above
  summary <- single_sequence$summary
  summary <- t(summary)
  summary <- as.data.frame(summary, stringsAsFactors = FALSE)
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
  
  # Input validation
  if (missing(readFileName) || !file.exists(readFileName)) {
    stop("Please provide a valid path to an .ab1 file")
  }
  
  message("Starting single sequence analysis...")
  
  # Generate an empty summary list to put the summary data in
  summarylist <- list()
  
  # Run Summarize.Single function to generate fasta file
  summarylist <- Summarize.Single(
    readFileName = readFileName, 
    readFeature  = readFeature, 
    summarylist  = summarylist
  )
  
  if (length(summarylist) == 0) {
    stop("Failed to process the sequence")
  }
  
  # Extract the summary data
  summary_data <- summarylist[[1]]
  
  # Change the order of the columns so the read name is the first column
  summary_data <- summary_data[, c(3, 1, 2)]
  
  # Create output directory if it doesn't exist
  results_dir <- "../Results/"
  if (!dir.exists(results_dir)) {
    dir.create(results_dir, recursive = TRUE)
  }
  
  # Export the summary data into a csv in the Results folder
  resultpath <- file.path(results_dir, paste0("Quality_Report_", basename(readFileName), ".csv"))
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
  
  # Input validation
  if (missing(x) || !file.exists(x)) {
    stop("Please provide a valid path to a FASTA file")
  }
  
  if (!requireNamespace("seqinr", quietly = TRUE)) {
    stop("Package 'seqinr' is required. Install with: install.packages('seqinr')")
  }
  
  message("Editing FASTA file: ", basename(x))
  
  tryCatch({
    fasta <- seqinr::read.fasta(x, whole.header = TRUE, as.string = TRUE)
    
    if (length(fasta) == 0) {
      warning("No sequences found in FASTA file: ", x)
      return(invisible(NULL))
    }
    
    # Replace spaces with underscores in headers
    noSpace <- lapply(names(fasta), function(y) gsub('\\s+', "_", y))
    
    seqinr::write.fasta(sequences = fasta, names = noSpace, file.out = x)
    
    message("Successfully edited ", length(fasta), " sequences in ", basename(x))
    
  }, error = function(e) {
    stop("Failed to edit FASTA file: ", conditionMessage(e))
  })
  
  invisible(NULL)
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
  
  # Input validation
  if (missing(x) || !file.exists(x)) {
    stop("Please provide a valid path to a query FASTA file")
  }
  
  if (missing(blast_db)) {
    stop("Please provide a valid path to a BLAST database")
  }
  
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required. Install with: install.packages('dplyr')")
  }
  
  if (!requireNamespace("tidyr", quietly = TRUE)) {
    stop("Package 'tidyr' is required. Install with: install.packages('tidyr')")
  }
  
  # Check if blastn is available
  blastn_check <- suppressWarnings(system2("which", blastn, stdout = TRUE, stderr = TRUE))
  if (length(blastn_check) == 0 || attr(blastn_check, "status") == 1) {
    stop("blastn command not found. Please install NCBI BLAST+ and set up the environment with setup_blast_env()")
  }
  
  message("Running BLAST search for: ", basename(x))
  message("  Database: ", blast_db)
  message("  E-value threshold: ", evalue)
  
  col_names <- c("Name",
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
  
  blast_out <- tryCatch({
    system2(
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
      stdout = TRUE,
      stderr = TRUE
    )
  }, error = function(e) {
    stop("BLAST command failed: ", conditionMessage(e))
  })
  
  # Check if BLAST returned an error
  if (!is.null(attr(blast_out, "status")) && attr(blast_out, "status") != 0) {
    stop("BLAST search failed. Error: ", paste(blast_out, collapse = "\n"))
  }
  
  if (length(blast_out) == 0) {
    warning("No BLAST hits found for ", basename(x))
    return(dplyr::tibble())
  }
  
  message("  Found ", length(blast_out), " BLAST hits")
  
  # Convert to tibble and parse
  blast_out <- dplyr::tibble(value = blast_out)
  blast_out <- tidyr::separate(blast_out, col = "value", into = col_names, 
                                sep = "\t", convert = TRUE, fill = "warn")
  
  message("BLAST search complete")
  
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
  
  # Input validation
  if (missing(file_name) || !file.exists(file_name)) {
    stop("Please provide a valid path to a FASTA file")
  }
  
  if (missing(blast_db)) {
    stop("Please provide a valid path to a BLAST database")
  }
  
  if (missing(DBname) || is.null(DBname)) {
    DBname <- "BLAST"
  }
  
  name <- basename(file_name)
  
  message("Processing: ", name)
  
  # Edit the FASTA file to remove spaces
  tryCatch({
    edit.fasta(file_name)
  }, error = function(e) {
    warning("Failed to edit FASTA file: ", conditionMessage(e))
  })
  
  # Run BLAST
  Blast_output <- tryCatch({
    Blast.CB(file_name, blast_db = blast_db)
  }, error = function(e) {
    stop("BLAST failed for ", name, ": ", conditionMessage(e))
  })
  
  # Create output directory if it doesn't exist
  results_dir <- "../Results/"
  if (!dir.exists(results_dir)) {
    dir.create(results_dir, recursive = TRUE)
    message("Created directory: ", results_dir)
  }
  
  # Save results
  mypath <- file.path(results_dir, paste0("Result_", DBname, "_", name, ".csv"))
  write.csv(Blast_output, file = mypath, row.names = FALSE)
  
  message("BLAST results saved to: ", mypath)
  message("  Number of hits: ", nrow(Blast_output))
  
  invisible(Blast_output)
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
                        DBname = "Analysis") {
  
  # Input validation
  if (missing(Blastpath) || !dir.exists(Blastpath)) {
    stop("Please provide a valid path to the directory containing FASTA files")
  }
  
  if (!requireNamespace("stringr", quietly = TRUE)) {
    stop("Package 'stringr' is required. Install with: install.packages('stringr')")
  }
  
  message("Starting batch BLAST analysis...")
  message("  Directory: ", Blastpath)
  
  # Find FASTA files
  file_names <- dir(Blastpath, pattern = ".fa|.fasta", full.names = TRUE)
  
  if (length(file_names) == 0) {
    stop("No FASTA files (.fa or .fasta) found in directory: ", Blastpath)
  }
  
  message("  Found ", length(file_names), " FASTA files")
  
  # Process 16S sequences
  message("\n=== Processing 16S Sequences ===")
  matches_16S <- stringr::str_match(file_names, "16S")
  indices_16S <- which(!is.na(matches_16S))
  file_names_16S <- file_names[indices_16S]
  
  if (length(file_names_16S) > 0) {
    message("Found ", length(file_names_16S), " 16S sequences")
    message("BLASTing against database: ", blast16Sdb)
    
    results_16S <- lapply(file_names_16S, function(f) {
      tryCatch({
        Blast.all(f, blast_db = blast16Sdb, DBname = paste0(DBname, "_16S"))
      }, error = function(e) {
        warning("Failed to BLAST ", basename(f), ": ", conditionMessage(e))
        return(NULL)
      })
    })
    
    success_16S <- sum(sapply(results_16S, function(x) !is.null(x)))
    message("Successfully processed ", success_16S, " out of ", length(file_names_16S), " 16S sequences")
  } else {
    message("No 16S sequences found (no files containing '16S' in filename)")
  }
  
  # Process ITS sequences
  message("\n=== Processing ITS Sequences ===")
  matches_ITS <- stringr::str_match(file_names, "ITS")
  indices_ITS <- which(!is.na(matches_ITS))
  file_names_ITS <- file_names[indices_ITS]
  
  if (length(file_names_ITS) > 0) {
    message("Found ", length(file_names_ITS), " ITS sequences")
    message("BLASTing against database: ", blastITSdb)
    
    results_ITS <- lapply(file_names_ITS, function(f) {
      tryCatch({
        Blast.all(f, blast_db = blastITSdb, DBname = paste0(DBname, "_ITS"))
      }, error = function(e) {
        warning("Failed to BLAST ", basename(f), ": ", conditionMessage(e))
        return(NULL)
      })
    })
    
    success_ITS <- sum(sapply(results_ITS, function(x) !is.null(x)))
    message("Successfully processed ", success_ITS, " out of ", length(file_names_ITS), " ITS sequences")
  } else {
    message("No ITS sequences found (no files containing 'ITS' in filename)")
  }
  
  message("\n=== BLAST Analysis Complete! ===")
  message("Total files processed: ", length(file_names_16S) + length(file_names_ITS))
  message("Results saved to: ../Results/")
  
  invisible(NULL)
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