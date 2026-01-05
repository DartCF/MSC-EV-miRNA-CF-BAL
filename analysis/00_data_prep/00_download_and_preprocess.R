# ==============================================================================
# Download and Preprocess miRNA Data from GEO
# ==============================================================================
#
# PURPOSE:
# Downloads raw count data from GEO (GSE282919), maps GEO sample positions to
# curated sample IDs, removes control probes and housekeeping genes, and 
# prepares data for downstream analysis.
#
# IMPORTANT NOTE ON SAMPLE MAPPING:
# The GEO file contains 81 samples numbered 1-81 by position. However, 3 samples
# (at GEO positions 27, 36, 46) were excluded during QC. The remaining 78 samples
# were renumbered as sample_1 through sample_78 in the curated analysis.
# 
# This script uses geo_to_curated_sample_mapping.csv to correctly map GEO 
# positions to curated sample IDs based on total read count matching.
#
# PREREQUISITES:
# Copy these files to data/metadata/ before running this script:
#   - sample_metadata.csv (treatment group assignments)
#   - geo_to_curated_sample_mapping.csv (GEO position to curated sample mapping)
#
# USAGE:
#   source(here("analysis", "00_data_prep", "00_download_and_preprocess.R"))
#
# OUTPUTS:
#   - data/processed/mirna_counts_filtered.csv   (2,083 miRNAs × 78 samples)
#
# GEO ACCESSION: GSE282919
# ==============================================================================

library(here)
library(readxl)
library(dplyr)
library(readr)

# ==============================================================================
# Configuration
# ==============================================================================

GEO_ACCESSION <- "GSE282919"
GEO_SUPP_FILE <- "GSE282919_VLP00436_miRNA_P33_21NOV2019_QC_and_Processed_Data.xlsx"

# Probes to remove (controls and housekeeping genes)
PROBES_TO_REMOVE <- c(
  # Control probes
  "CTRL_ANT1", "CTRL_ANT2", "CTRL_ANT3", "CTRL_ANT4", "CTRL_ANT5", 
  "CTRL_miR_POS",
  # Housekeeping genes
  "HK_ACTB", "HK_B2M", "HK_GAPDH", "HK_PPIA", "HK_RNU47", "HK_RNU75",
  "HK_RNY3", "HK_RPL19", "HK_RPL27", "HK_RPS12", "HK_RPS20", 
  "HK_SNORA66", "HK_YWHAZ"
)

# GEO positions excluded during QC (identified by total count matching)
# These are at positions 27, 36, 46 in the GEO file, NOT at the end
GEO_POSITIONS_EXCLUDED_QC <- c(27, 36, 46)

# ==============================================================================
# Create directories
# ==============================================================================

dir.create(here("data", "raw"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("data", "processed"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("data", "metadata"), recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# Step 1: Check for required metadata files
# ==============================================================================

metadata_file <- here("data", "metadata", "sample_metadata.csv")
mapping_file <- here("data", "metadata", "geo_to_curated_sample_mapping.csv")

if (!file.exists(metadata_file)) {
  stop(
    "\n",
    "================================================================================\n",
    "ERROR: sample_metadata.csv not found!\n",
    "================================================================================\n",
    "\n",
    "Please copy sample_metadata.csv to:\n",
    "  ", metadata_file, "\n",
    "\n",
    "This file contains the treatment group assignments and should have been\n",
    "provided with the repository setup materials.\n",
    "================================================================================\n"
  )
}

if (!file.exists(mapping_file)) {
  stop(
    "\n",
    "================================================================================\n",
    "ERROR: geo_to_curated_sample_mapping.csv not found!\n",
    "================================================================================\n",
    "\n",
    "Please copy geo_to_curated_sample_mapping.csv to:\n",
    "  ", mapping_file, "\n",
    "\n",
    "This file maps GEO sample positions (1-81) to curated sample IDs (sample_1\n",
    "through sample_78). It was created by matching total read counts between\n",
    "the GEO file and the original curated analysis.\n",
    "================================================================================\n"
  )
}

sample_metadata <- read_csv(metadata_file, show_col_types = FALSE)
cat("Loaded sample metadata:", nrow(sample_metadata), "samples\n")

geo_mapping <- read_csv(mapping_file, show_col_types = FALSE)
cat("Loaded GEO-to-curated sample mapping:", nrow(geo_mapping), "GEO positions\n")

# Verify mapping file structure
required_cols <- c("geo_position", "curated_sample", "in_curated_analysis")
missing_cols <- setdiff(required_cols, names(geo_mapping))
if (length(missing_cols) > 0) {
  stop("GEO mapping file is missing required columns: ", paste(missing_cols, collapse = ", "))
}

# Report on excluded samples
n_excluded <- sum(!geo_mapping$in_curated_analysis)
excluded_positions <- geo_mapping$geo_position[!geo_mapping$in_curated_analysis]
cat(sprintf("  GEO positions excluded from analysis: %s\n", 
            paste(excluded_positions, collapse = ", ")))

# ==============================================================================
# Step 2: Download from GEO (if not already present)
# ==============================================================================

geo_file <- here("data", "raw", GEO_SUPP_FILE)

if (!file.exists(geo_file)) {
  cat("Downloading supplementary file from GEO...\n")
  
  # Construct download URL
  geo_url <- paste0(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE282nnn/", GEO_ACCESSION, 
    "/suppl/", GEO_SUPP_FILE
  )
  
  download.file(geo_url, destfile = geo_file, mode = "wb")
  cat("Download complete:", geo_file, "\n\n")
} else {
  cat("GEO supplementary file already exists:", geo_file, "\n\n")
}

# ==============================================================================
# Step 3: Load and parse raw data
# ==============================================================================

cat("Loading raw data from GEO file...\n")

# Read the Raw sheet (has non-standard header structure)
raw_data <- read_excel(geo_file, sheet = "Raw", col_names = FALSE)

# File structure (1-indexed):
# Row 1: Investigator Name
# Row 7: Sample ID (1, 2, 3, ...)
# Row 8: Well
# Row 9: Sample Name
# Row 10: Total Counts
# Row 11+: miRNA count data

# Verify the number of samples matches the mapping
n_geo_samples <- nrow(geo_mapping)
n_data_cols <- ncol(raw_data) - 1  # subtract Gene_ID column

if (n_data_cols != n_geo_samples) {
  warning(sprintf("GEO file has %d sample columns but mapping has %d entries", 
                  n_data_cols, n_geo_samples))
}

# Extract count data (starts at row 11)
# Column 1 is probe names, columns 2+ are sample data
probe_names <- as.character(raw_data[[1]][11:nrow(raw_data)])
count_matrix <- as.matrix(raw_data[11:nrow(raw_data), 2:(n_geo_samples + 1)])

# Convert to numeric
mode(count_matrix) <- "numeric"
rownames(count_matrix) <- probe_names

# Assign column names using the mapping (geo_position -> curated_sample)
# Use EXCLUDED_geo_N as temporary name for excluded samples
col_names <- ifelse(
  is.na(geo_mapping$curated_sample),
  paste0("EXCLUDED_geo_", geo_mapping$geo_position),
  geo_mapping$curated_sample
)
colnames(count_matrix) <- col_names

cat(sprintf("  Raw data: %d probes × %d samples\n", nrow(count_matrix), ncol(count_matrix)))

# ==============================================================================
# Step 4: Remove QC-failed samples
# ==============================================================================

# Identify excluded samples (those without curated_sample mapping)
samples_to_remove <- col_names[grepl("^EXCLUDED_", col_names)]
n_removed <- length(samples_to_remove)

cat(sprintf("\nRemoving %d QC-failed samples at GEO positions: %s\n", 
            n_removed, 
            paste(GEO_POSITIONS_EXCLUDED_QC, collapse = ", ")))

count_matrix <- count_matrix[, !colnames(count_matrix) %in% samples_to_remove]
cat(sprintf("  After QC filtering: %d samples\n", ncol(count_matrix)))

# ==============================================================================
# Step 5: Remove Total Counts row, control probes, and housekeeping genes
# ==============================================================================

# Remove "Total Counts" row if present
if ("Total Counts" %in% rownames(count_matrix)) {
  count_matrix <- count_matrix[rownames(count_matrix) != "Total Counts", ]
  cat("\nRemoved 'Total Counts' row\n")
}

# Remove control probes and housekeeping genes
probes_present <- intersect(PROBES_TO_REMOVE, rownames(count_matrix))
cat(sprintf("Removing %d control/housekeeping probes:\n", length(probes_present)))
for (p in probes_present) {
  cat(sprintf("  - %s\n", p))
}

count_matrix <- count_matrix[!rownames(count_matrix) %in% PROBES_TO_REMOVE, ]
cat(sprintf("  After filtering: %d miRNAs\n", nrow(count_matrix)))

# ==============================================================================
# Step 6: Verify sample alignment with metadata
# ==============================================================================

cat("\nVerifying sample alignment with metadata...\n")

samples_in_counts <- colnames(count_matrix)
samples_in_metadata <- sample_metadata$Sample

# Check for mismatches
in_counts_not_meta <- setdiff(samples_in_counts, samples_in_metadata)
in_meta_not_counts <- setdiff(samples_in_metadata, samples_in_counts)

if (length(in_counts_not_meta) > 0) {
  warning("Samples in count matrix but not in metadata: ", 
          paste(in_counts_not_meta, collapse = ", "))
}

if (length(in_meta_not_counts) > 0) {
  cat("  Samples in metadata but not in counts (expected - QC failures): ", 
      paste(in_meta_not_counts, collapse = ", "), "\n")
}

# Verify all count samples have metadata
if (all(samples_in_counts %in% samples_in_metadata)) {
  cat("  ✓ All samples in count matrix have metadata\n")
} else {
  stop("Some samples in count matrix are missing from metadata!")
}

# ==============================================================================
# Step 7: Validate mapping by checking total counts
# ==============================================================================

cat("\nValidating GEO-to-curated mapping by total counts...\n")

# Calculate total counts for the filtered data
calculated_totals <- colSums(count_matrix)

# Get expected totals from mapping (for curated samples only)
mapping_curated <- geo_mapping[geo_mapping$in_curated_analysis, ]
expected_totals <- setNames(mapping_curated$curated_total, mapping_curated$curated_sample)

# Compare
mismatches <- 0
for (sample in names(calculated_totals)) {
  if (sample %in% names(expected_totals)) {
    if (abs(calculated_totals[sample] - expected_totals[sample]) > 1) {
      cat(sprintf("  WARNING: Total count mismatch for %s: calculated=%d, expected=%d\n",
                  sample, calculated_totals[sample], expected_totals[sample]))
      mismatches <- mismatches + 1
    }
  }
}

if (mismatches == 0) {
  cat("  ✓ All sample total counts match the mapping file\n")
} else {
  warning(sprintf("%d samples have total count mismatches", mismatches))
}

# ==============================================================================
# Step 8: Save processed data
# ==============================================================================

cat("\nSaving processed data...\n")

# Save filtered count matrix
count_df <- as.data.frame(count_matrix)
count_df$Gene_ID <- rownames(count_matrix)
count_df <- count_df[, c("Gene_ID", colnames(count_matrix))]

write_csv(count_df, here("data", "processed", "mirna_counts_filtered.csv"))
cat("  ✓ data/processed/mirna_counts_filtered.csv\n")

# ==============================================================================
# Step 9: Summary
# ==============================================================================

cat("\n")
cat("================================================================================\n")
cat("                         DATA PREPROCESSING COMPLETE                            \n")
cat("================================================================================\n\n")

cat("INPUT:\n")
cat(sprintf("  GEO Accession: %s\n", GEO_ACCESSION))
cat("  Raw samples in GEO: 81\n")
cat(sprintf("  QC-excluded samples (GEO positions): %s\n", 
            paste(GEO_POSITIONS_EXCLUDED_QC, collapse = ", ")))
cat("  Raw probes: ~2,103\n\n")

cat("OUTPUT:\n")
cat(sprintf("  Samples: %d (after removing %d QC failures)\n", 
            ncol(count_matrix), n_removed))
cat(sprintf("  miRNAs: %d (after removing %d controls/housekeeping)\n", 
            nrow(count_matrix), length(PROBES_TO_REMOVE)))
cat("\n")

cat("SAMPLE GROUPS (from metadata):\n")
metadata_subset <- sample_metadata %>% filter(Sample %in% samples_in_counts)
group_counts <- table(metadata_subset$Treatment_Group)
manuscript_groups <- c("CF_Asp_Neg", "HC", "UNSTIM")

for (grp in names(group_counts)) {
  n <- group_counts[grp]
  in_study <- ifelse(grp %in% manuscript_groups, " [IN CURRENT STUDY]", "")
  cat(sprintf("  %s: %d%s\n", grp, n, in_study))
}

n_in_study <- sum(metadata_subset$Treatment_Group %in% manuscript_groups)
cat(sprintf("\n  Total samples in current study: %d\n", n_in_study))

cat("\n")
cat("NEXT STEPS:\n")
cat("  1. Run analysis/01_compositional/CompositionalAnalysis.Rmd\n")
cat("\n")
cat("================================================================================\n")
