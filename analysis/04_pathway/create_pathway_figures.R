# ==============================================================================
# CREATE PATHWAY FIGURES FOR MANUSCRIPT
# ==============================================================================
# Generates all pathway analysis figures and tables for the manuscript.
# This script orchestrates the pathway visualization scripts and can be run
# independently or called from the root run_all_analyses.R.
#
# PREREQUISITES:
#   - Step 1 (Compositional Analysis) must be complete
#   - results/01_compositional/major_miRNAs_1pct_summary.csv must exist
#
# OUTPUTS:
#   - Fig 4:  miRNA target overlap summary (mirna_overlap_summary.png)
#   - Fig 5A: Pathway counts barplot (pathway_counts_barplot.png)
#   - Fig 5B: GO BP enrichment dotplot (CF_major_stringent_GO_BP_dotplot.png)
#   - Fig 5C: Pathway categories stacked barplot (pathway_categories_stacked_barplot.png)
#   - Table 1: Pathway Category Summary (pathway_categories_summary.csv)
#   - Tables 2, 3: Convergent Genes (gene_pathway_counts.csv)
#
# USAGE:
#   source(here("analysis", "04_pathway", "create_pathway_figures.R"))
#
# ==============================================================================

library(here)

cat("\n")
cat("================================================================================\n")
cat("     CREATE PATHWAY FIGURES FOR MANUSCRIPT\n")
cat("================================================================================\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

pathway_dir <- here("analysis", "04_pathway")

# ==============================================================================
# STEP 1: Run Pathway Enrichment (generates Fig 5B)
# ==============================================================================
# This sources the factorial analysis framework and runs Analysis C
# (Major miRNAs, Stringent filtering)

cat("Loading pathway analysis framework...\n")
source(file.path(pathway_dir, "run_factorial_analysis.R"))

cat("\nRunning pathway enrichment for Major miRNAs (Stringent)...\n")
cat("This generates Fig 5B: GO BP enrichment dotplot\n\n")

results_C <- run_analysis("C")

cat("\n✓ Pathway enrichment complete\n")
cat("  Output: Fig 5B (CF_major_stringent_GO_BP_dotplot.png)\n\n")

# ==============================================================================
# STEP 2: Generate Additional Pathway Visualizations
# ==============================================================================

cat("Running pathway visualization scripts...\n\n")

# ------------------------------------------------------------------------------
# Fig 5A: Pathway counts barplot
# Search term: pathway_counts_barplot
# ------------------------------------------------------------------------------
cat("--- Fig 5A: Pathway counts barplot ---\n")
source(file.path(pathway_dir, "plot_pathway_counts.R"))
cat("✓ plot_pathway_counts.R complete\n\n")

# ------------------------------------------------------------------------------
# Fig 5C, Table 1: Pathway categories
# Search terms: pathway_categories_stacked_barplot, pathway_categories_summary
# ------------------------------------------------------------------------------
cat("--- Fig 5C & Table 1: Pathway categories ---\n")
source(file.path(pathway_dir, "categorize_pathways.R"))
cat("✓ categorize_pathways.R complete\n\n")

# ------------------------------------------------------------------------------
# Fig 4: miRNA target overlap
# Search term: mirna_overlap_summary
# ------------------------------------------------------------------------------
cat("--- Fig 4: miRNA target overlap ---\n")
source(file.path(pathway_dir, "analyze_mirna_target_overlap.R"))
cat("✓ analyze_mirna_target_overlap.R complete\n\n")

# ------------------------------------------------------------------------------
# Tables 2, 3: Convergent genes (pathway similarity)
# Search term: gene_pathway_counts
# ------------------------------------------------------------------------------
cat("--- Tables 2 & 3: Convergent genes ---\n")
source(file.path(pathway_dir, "analyze_pathway_similarity.R"))
cat("✓ analyze_pathway_similarity.R complete\n\n")

# ==============================================================================
# SUMMARY
# ==============================================================================

cat("================================================================================\n")
cat("     PATHWAY FIGURES COMPLETE\n")
cat("================================================================================\n\n")

cat("Finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("FIGURES GENERATED:\n")
cat("  ✓ Fig 4:  mirna_overlap_summary.png/.pdf\n")
cat("  ✓ Fig 5A: pathway_counts_barplot.png/.pdf\n")
cat("  ✓ Fig 5B: CF_major_stringent_GO_BP_dotplot.png/.pdf\n")
cat("  ✓ Fig 5C: pathway_categories_stacked_barplot.png/.pdf\n")
cat("\n")

cat("TABLES GENERATED:\n")
cat("  ✓ Table 1: pathway_categories_summary.csv\n")
cat("  ✓ Tables 2, 3: gene_pathway_counts.csv (filter for inflammation/remodeling)\n")
cat("\n")

cat("================================================================================\n\n")
