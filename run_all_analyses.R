# ==============================================================================
# MASTER ANALYSIS SCRIPT: MSC-EV-miRNA-CF-BAL
# ==============================================================================
# Reproduces the complete analysis pipeline from raw data to manuscript figures.
#
# USAGE:
#   source(here("run_all_analyses.R"))
#
# PREREQUISITES:
#   - Copy metadata files to data/metadata/:
#       - sample_metadata.csv
#       - geo_to_curated_sample_mapping.csv
#       - CF_miRNA_annotations.csv
#
# OUTPUTS:
#   - All figures (PNG/PDF) and tables (CSV) in results/ directories
#   - Console output summarizing each step
#
# GEO ACCESSION: GSE282919
# ==============================================================================

library(here)

cat("\n")
cat("================================================================================\n")
cat("     MSC-EV-miRNA-CF-BAL: FULL ANALYSIS PIPELINE\n")
cat("================================================================================\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# ==============================================================================
# STEP 0: Data Preprocessing
# ==============================================================================
# Downloads from GEO and applies sample mapping to fix QC-excluded sample positions

cat("\n")
cat("================================================================================\n")
cat("     STEP 0: Data Preprocessing\n")
cat("================================================================================\n\n")

source(here("analysis", "00_data_prep", "00_download_and_preprocess.R"))

cat("✓ Step 0 complete: Data preprocessed\n")

# ==============================================================================
# STEP 1: Compositional Analysis
# ==============================================================================
# Generates major miRNA list, fractional abundances, treatment comparisons
# Outputs: Figs 1, 2, 3A, 3B, 3C

cat("\n")
cat("================================================================================\n")
cat("     STEP 1: Compositional Analysis\n")
cat("================================================================================\n\n")

# Main compositional analysis (Figs 1, 2, 3A, 3B)
rmarkdown::render(
  here("analysis", "01_compositional", "CompositionalAnalysis.Rmd"),
  output_dir = here("results", "01_compositional"),
  quiet = FALSE
)

cat("✓ CompositionalAnalysis.Rmd rendered\n")

# Simpson diversity figure (Fig 3C)
source(here("analysis", "01_compositional", "simpson_diversity_figure.R"))

cat("✓ Step 1 complete: Compositional analysis finished\n")
cat("  Outputs: Figs 1, 2, 3A, 3B, 3C\n")
cat("  Critical output: results/01_compositional/major_miRNAs_1pct_summary.csv\n")

# ==============================================================================
# STEP 2: Minor miRNA Analysis
# ==============================================================================
# CF-relevant miRNAs, prevalence analysis, volcano plots
# Outputs: Fig 6, Table 4

cat("\n")
cat("================================================================================\n")
cat("     STEP 2: Minor miRNA Analysis\n")
cat("================================================================================\n\n")

rmarkdown::render(
  here("analysis", "02_minor_miRNA", "MinorMiRNA_Prevalence_Analysis.Rmd"),
  output_dir = here("results", "02_minor_miRNA"),
  quiet = FALSE
)

cat("✓ Step 2 complete: Minor miRNA analysis finished\n")
cat("  Outputs: Fig 6, Table 4 (CF_relevant_miRNAs_annotated.csv)\n")

# ==============================================================================
# STEP 3: DESeq2 Differential Expression Analysis
# ==============================================================================
# CF vs HC comparison, abundance tier classification
# Outputs: Fig 7

cat("\n")
cat("================================================================================\n")
cat("     STEP 3: DESeq2 Analysis\n")
cat("================================================================================\n\n")

source(here("analysis", "03_deseq2", "DESeq2_CF_vs_HC_analysis.R"))

cat("✓ Step 3 complete: DESeq2 analysis finished\n")
cat("  Outputs: Fig 7\n")

# ==============================================================================
# STEP 4: Pathway Analysis
# ==============================================================================
# IMPORTANT: Must run AFTER Step 1 generates major_miRNAs_1pct_summary.csv
# Outputs: Figs 4, 5A, 5B, 5C; Tables 1, 2, 3

cat("\n")
cat("================================================================================\n")
cat("     STEP 4: Pathway Analysis\n")
cat("================================================================================\n\n")

source(here("analysis", "04_pathway", "create_pathway_figures.R"))

cat("✓ Step 4 complete: Pathway analysis finished\n")
cat("  Outputs: Figs 4, 5A, 5B, 5C; Tables 1, 2, 3\n")

# ==============================================================================
# SUMMARY
# ==============================================================================

cat("\n")
cat("================================================================================\n")
cat("     PIPELINE COMPLETE\n")
cat("================================================================================\n\n")

cat("Finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("FIGURES GENERATED:\n")
cat("  ✓ Fig 1:  Major miRNA stacked bar by sample\n")
cat("  ✓ Fig 2:  Major miRNA mean fractional abundance\n")
cat("  ✓ Fig 3A: Major miRNAs with significant treatment differences\n")
cat("  ✓ Fig 3B: Density distribution by treatment\n")
cat("  ✓ Fig 3C: Simpson diversity by treatment\n")
cat("  ✓ Fig 4:  miRNA target overlap summary\n")
cat("  ✓ Fig 5A: Pathways enriched in major miRNA targets\n")
cat("  ✓ Fig 5B: GO Biological Process enrichment dotplot\n")
cat("  ✓ Fig 5C: CF pathway categories stacked barplot\n")
cat("  ✓ Fig 6:  Volcano plot CF-relevant miRNAs\n")
cat("  ✓ Fig 7:  DESeq2 volcano faceted by abundance tier\n")
cat("\n")

cat("TABLES GENERATED:\n")
cat("  ✓ Table 1: Pathway Category Summary\n")
cat("  ✓ Table 2: Convergent Genes - Inflammation (gene_pathway_counts.csv)\n")
cat("  ✓ Table 3: Convergent Genes - Tissue Remodeling (gene_pathway_counts.csv)\n")
cat("  ✓ Table 4: CF-Relevant miRNAs (CF_relevant_miRNAs_annotated.csv)\n")
cat("  ✓ Table 5: Context-Dependent Effects of miR-21 (manually curated)\n")
cat("\n")

cat("================================================================================\n")
cat("To archive outputs for submission, run:\n")
cat("  source(here('archive_for_submission.R'))\n")
cat("================================================================================\n\n")
