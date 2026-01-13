# ==============================================================================
# Archive Outputs for Manuscript Submission
# ==============================================================================
# Collects all figures, tables, and supporting files into a single archive
# for sharing with collaborators or journal submission.
#
# This script reflects the final manuscript structure (5 figures, 4 tables)
# as agreed with co-authors in January 2026.
#
# USAGE:
#   source(here("archive_for_submission.R"))
#
# OUTPUT:
#   - submission_archive_YYYYMMDD.zip in project root
# ==============================================================================

library(here)

# Create timestamp for archive name
timestamp <- format(Sys.time(), "%Y%m%d")
archive_name <- paste0("submission_archive_", timestamp)
archive_dir <- here(archive_name)

cat("\n")
cat("================================================================================\n")
cat("     Creating Submission Archive\n")
cat("================================================================================\n\n")

# Create temporary directory for archive contents
dir.create(archive_dir, showWarnings = FALSE)
dir.create(file.path(archive_dir, "figures"), showWarnings = FALSE)
dir.create(file.path(archive_dir, "tables"), showWarnings = FALSE)

# ==============================================================================
# FIGURE MANIFEST
# ==============================================================================
#
# Fig 1: Major miRNA composition stacked by sample
# Fig 2: Major miRNA mean fractional abundance  
# Fig 3: Major miRNAs with significant treatment differences
# Fig 4: miRNA target overlap summary (4 panels: A, B, C, D)
# Fig 5: CF-relevant pathway categories (multi-database)
#
# ==============================================================================

cat("Collecting figures...\n")

figures <- list(
  # Compositional analysis (from results/01_compositional/)
  "Fig_1" = here("results", "01_compositional", "major_miRNAs_stacked_by_sample.png"),
  "Fig_2" = here("results", "01_compositional", "major_miRNAs_barplot.png"),
  "Fig_3" = here("results", "01_compositional", "major_miRNAs_significant_differences.png"),
  

  # Pathway analysis (from results/04_pathway/)
  "Fig_4" = here("results", "04_pathway", "major", "stringent", "mirna_overlap_summary.png"),
  "Fig_5" = here("results", "04_pathway", "pathway_categories_stacked_barplot.png")
)

for (fig_name in names(figures)) {
  src <- figures[[fig_name]]
  if (file.exists(src)) {
    file.copy(src, file.path(archive_dir, "figures", paste0(fig_name, ".png")))
    # Also copy PDF if exists
    pdf_src <- sub("\\.png$", ".pdf", src)
    if (file.exists(pdf_src)) {
      file.copy(pdf_src, file.path(archive_dir, "figures", paste0(fig_name, ".pdf")))
    }
    cat("
  ✓", fig_name, "\n")
  } else {
    cat("  ✗", fig_name, "(not found:", src, ")\n")
  }
}

# ==============================================================================
# TABLE MANIFEST
# ==============================================================================
#
# Table 1: Pathway category summary (counts by functional category)
# Table 2: Convergent genes - Inflammation pathways
# Table 3: Convergent genes - Tissue remodeling pathways
# Table 4: CF-relevant miRNAs in hMSC EV cargo (with literature references)
#
# Note: Tables 2 and 3 derive from the same source file (gene_pathway_counts.csv)
#       filtered by pathway category.
#
# ==============================================================================

cat("\nCollecting tables...\n")

tables <- list(
  "Table_1_pathway_categories" = here("results", "04_pathway", "pathway_categories_summary.csv"),
  "Table_2_3_convergent_genes" = here("results", "04_pathway", "major", "stringent", 
                                       "pathway_similarity_results", "gene_pathway_counts.csv"),
  "Table_4_CF_relevant_miRNAs" = here("results", "02_minor_miRNA", "CF_relevant_miRNAs_annotated.csv")
)

for (tbl_name in names(tables)) {
  src <- tables[[tbl_name]]
  if (file.exists(src)) {
    file.copy(src, file.path(archive_dir, "tables", paste0(tbl_name, ".csv")))
    cat("  ✓", tbl_name, "\n")
  } else {
    cat("  ✗", tbl_name, "(not found:", src, ")\n")
  }
}

# Also copy Word document for Table 4 if exists
docx_src <- here("results", "02_minor_miRNA", "CF_Relevant_miRNAs_Table.docx")
if (file.exists(docx_src)) {
  file.copy(docx_src, file.path(archive_dir, "tables", "Table_4_CF_relevant_miRNAs.docx"))
  cat("  ✓ Table_4 Word document\n")
}

# ==============================================================================
# Copy Manifest
# ==============================================================================

cat("\nAdding manifest...\n")

manifest_src <- here("manuscript_manifest.csv")
if (file.exists(manifest_src)) {
  file.copy(manifest_src, file.path(archive_dir, "manuscript_manifest.csv"))
  cat("  ✓ manuscript_manifest.csv\n")
} else {
  cat("  ⚠ manuscript_manifest.csv not found in project root\n")
}

# ==============================================================================
# Create ZIP Archive
# ==============================================================================

cat("\nCreating ZIP archive...\n")

zip_file <- here(paste0(archive_name, ".zip"))

# Use R's zip function
old_wd <- setwd(here())
zip(zip_file, files = archive_name, flags = "-r")
setwd(old_wd)

# Clean up temporary directory
unlink(archive_dir, recursive = TRUE)

cat("\n")
cat("================================================================================\n")
cat("     Archive Created\n")
cat("================================================================================\n\n")

cat("Archive location:", zip_file, "\n\n")

cat("Contents:\n")
cat("  figures/     - 5 figures (PNG + PDF)\n")
cat("  tables/      - 4 tables (CSV + Word)\n")
cat("  manuscript_manifest.csv - Figure/table to script mapping\n")

cat("\n")
cat("Ready for submission to Cytotherapy!\n")
cat("================================================================================\n\n")
