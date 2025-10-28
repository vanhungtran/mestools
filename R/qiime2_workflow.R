# QIIME2 Workflow in R
# Author: Generated from Qiime2_pipeline.ipynb
# Date: 2025-10-28

# ==============================================================================
# Setup: Helper Functions to Run QIIME2 Commands from R
# ==============================================================================

# Option 1: Simple system() approach
run_qiime <- function(command) {
  conda_env <- "qiime2-amplicon-2025.7"
  conda_path <- "/Users/vhtran/anaconda3"
  
  full_command <- sprintf(
    "source %s/etc/profile.d/conda.sh && conda activate %s && %s",
    conda_path, conda_env, command
  )
  system(full_command, intern = FALSE)
}

# Option 2: More robust system2() approach with error handling
run_qiime2 <- function(command, verbose = TRUE) {
  bash_script <- sprintf(
    "source /Users/vhtran/anaconda3/etc/profile.d/conda.sh && conda activate qiime2-amplicon-2025.7 && %s",
    command
  )
  
  result <- system2("bash", 
                    args = c("-c", bash_script),
                    stdout = TRUE, 
                    stderr = TRUE)
  
  if (verbose) {
    cat(result, sep = "\n")
  }
  
  invisible(result)
}

# ==============================================================================
# Import Data Files
# ==============================================================================

cat("Importing phylogenetic tree...\n")
run_qiime2("qiime tools import \\
  --input-path /Users/vhtran/Downloads/Workshop/FromQiita/pax/132071_insertion_tree.relabelled.tre \\
  --output-path pax_insertion_tree.qza \\
  --type 'Phylogeny[Rooted]'")

cat("Importing feature table...\n")
run_qiime2("qiime tools import \\
  --input-path /Users/vhtran/Downloads/Workshop/FromQiita/pax/132071_reference-hit.biom \\
  --type 'FeatureTable[Frequency]' \\
  --input-format BIOMV210Format \\
  --output-path pax_feature_table.qza")

cat("Importing representative sequences...\n")
run_qiime2("qiime tools import \\
  --input-path /Users/vhtran/Downloads/Workshop/FromQiita/pax/132071_reference-hit.seqs.fa \\
  --output-path representative-sequences.qza \\
  --type 'FeatureData[Sequence]'")

# ==============================================================================
# Feature Table Analysis
# ==============================================================================

cat("Summarizing feature table...\n")
run_qiime2("qiime feature-table summarize \\
  --i-table pax_feature_table.qza \\
  --o-visualization pax_feature_table.qzv")

cat("Filtering features (min frequency = 10)...\n")
run_qiime2("qiime feature-table filter-features \\
  --i-table pax_feature_table.qza \\
  --p-min-frequency 10 \\
  --o-filtered-table filtered_pax_feature_table.qza")

cat("Summarizing filtered feature table...\n")
run_qiime2("qiime feature-table summarize \\
  --i-table filtered_pax_feature_table.qza \\
  --o-visualization filtered_pax_feature_table.qzv")

# ==============================================================================
# Taxonomic Classification
# ==============================================================================

cat("Running taxonomic classification...\n")
run_qiime2("qiime feature-classifier classify-sklearn \\
  --i-reads representative-sequences.qza \\
  --i-classifier /Users/vhtran/Downloads/Workshop/2024.09.backbone.v4.nb.qza \\
  --o-classification classification.qza")

cat("Creating taxonomy barplot...\n")
run_qiime2("qiime taxa barplot \\
  --i-table pax_feature_table.qza \\
  --i-taxonomy classification.qza \\
  --m-metadata-file /Users/vhtran/Downloads/Workshop/sequences/pax/mapping_files/pax_metadata.txt \\
  --o-visualization taxa_boxplot.qzv")

cat("Filtering mitochondria and chloroplast...\n")
run_qiime2("qiime taxa filter-table \\
  --i-table filtered_pax_feature_table.qza \\
  --i-taxonomy classification.qza \\
  --p-exclude mitochondria,chloroplast \\
  --o-filtered-table double_filtered_pax_feature_table.qza")

cat("Summarizing double-filtered feature table...\n")
run_qiime2("qiime feature-table summarize \\
  --i-table double_filtered_pax_feature_table.qza \\
  --o-visualization double_filtered_pax_feature_table.qzv")

# ==============================================================================
# Diversity Analysis
# ==============================================================================

cat("Running alpha rarefaction...\n")
run_qiime2("qiime diversity alpha-rarefaction \\
  --i-table double_filtered_pax_feature_table.qza \\
  --i-phylogeny pax_insertion_tree.qza \\
  --m-metadata-file /Users/vhtran/Downloads/Workshop/sequences/pax/mapping_files/pax_metadata.txt \\
  --p-min-depth 30000 \\
  --p-max-depth 65000 \\
  --p-steps 10 \\
  --p-iterations 10 \\
  --o-visualization pax_rarefactioncurves.qzv")

cat("Rarefying feature table...\n")
run_qiime2("qiime feature-table rarefy \\
  --i-table double_filtered_pax_feature_table.qza \\
  --p-sampling-depth 46000 \\
  --o-rarefied-table pax_rarefied-feature-table.qza")

cat("Calculating beta diversity (UniFrac)...\n")
run_qiime2("qiime diversity beta-phylogenetic \\
  --i-table pax_rarefied-feature-table.qza \\
  --i-phylogeny pax_insertion_tree.qza \\
  --p-metric 'unweighted_unifrac' \\
  --o-distance-matrix pax_unweighted_distance_matrix.qza")

cat("Running PCoA...\n")
run_qiime2("qiime diversity pcoa \\
  --i-distance-matrix pax_unweighted_distance_matrix.qza \\
  --o-pcoa pax_unweighted_pcoa.qza")

cat("Creating Emperor plot...\n")
run_qiime2("qiime emperor plot \\
  --i-pcoa pax_unweighted_pcoa.qza \\
  --m-metadata-file /Users/vhtran/Downloads/Workshop/sequences/pax/mapping_files/pax_metadata.txt \\
  --o-visualization pax_emperor.qzv")

# ==============================================================================
# Option 3: Read QIIME2 Artifacts into R using qiime2R
# ==============================================================================

cat("\n\n=== Reading QIIME2 artifacts into R ===\n")

# Install qiime2R if needed
if (!require("qiime2R", quietly = TRUE)) {
  if (!require("devtools", quietly = TRUE)) install.packages("devtools")
  devtools::install_github("jbisanz/qiime2R")
}

# Load required libraries
library(qiime2R)
library(phyloseq)
library(ggplot2)
library(dplyr)

# Read QIIME2 artifacts
cat("Reading feature table...\n")
feature_table <- read_qza("pax_feature_table.qza")

cat("Reading taxonomy...\n")
taxonomy <- read_qza("classification.qza")

cat("Reading phylogenetic tree...\n")
tree <- read_qza("pax_insertion_tree.qza")

cat("Reading metadata...\n")
metadata <- read.table(
  "/Users/vhtran/Downloads/Workshop/sequences/pax/mapping_files/pax_metadata.txt",
  sep = "\t", 
  header = TRUE, 
  row.names = 1,
  comment.char = ""
)

# Create phyloseq object
cat("Creating phyloseq object...\n")
physeq <- phyloseq(
  otu_table(feature_table$data, taxa_are_rows = TRUE),
  tax_table(as.matrix(taxonomy$data)),
  phy_tree(tree$data),
  sample_data(metadata)
)

cat("\nPhyloseq object created successfully!\n")
print(physeq)

# ==============================================================================
# Example Visualizations in R
# ==============================================================================

cat("\n\n=== Creating visualizations ===\n")

# Alpha diversity plot
p1 <- plot_richness(physeq, measures = c("Shannon", "Simpson")) +
  theme_bw() +
  labs(title = "Alpha Diversity")

# Save plot
ggsave("alpha_diversity.pdf", p1, width = 10, height = 6)
cat("Saved: alpha_diversity.pdf\n")

# Ordination plot
ord <- ordinate(physeq, method = "PCoA", distance = "unifrac")
p2 <- plot_ordination(physeq, ord, type = "samples") +
  theme_bw() +
  labs(title = "PCoA - UniFrac Distance")

ggsave("pcoa_plot.pdf", p2, width = 10, height = 8)
cat("Saved: pcoa_plot.pdf\n")

# Taxonomy barplot
p3 <- plot_bar(physeq, fill = "Phylum") +
  theme_bw() +
  labs(title = "Taxonomic Composition")

ggsave("taxonomy_barplot.pdf", p3, width = 14, height = 8)
cat("Saved: taxonomy_barplot.pdf\n")

cat("\n=== Workflow complete! ===\n")
