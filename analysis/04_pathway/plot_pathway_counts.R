# ==============================================================================
# Barplot of Significant Pathways by Database
# ==============================================================================
# Counts and visualizes enriched pathways (q < 0.05) from pathway enrichment analysis
#
# USAGE:
#   source(here("analysis", "04_pathway", "plot_pathway_counts.R"))
#
# INPUT:
#   results/04_pathway/major/stringent/CF_major_stringent_*_enrichment.csv
#
# OUTPUT:
#   results/04_pathway/pathway_counts_barplot.png (Fig 5A)
#   results/04_pathway/pathway_counts_barplot.pdf
#
# ==============================================================================

library(here)

# ------------------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------------------

input_dir <- here("results", "04_pathway", "major", "stringent")
output_dir <- here("results", "04_pathway")

files <- list(
  "GO_BP" = file.path(input_dir, "CF_major_stringent_GO_BP_enrichment.csv"),
  "GO_CC" = file.path(input_dir, "CF_major_stringent_GO_CC_enrichment.csv"),
  "GO_MF" = file.path(input_dir, "CF_major_stringent_GO_MF_enrichment.csv"),
  "KEGG" = file.path(input_dir, "CF_major_stringent_KEGG_enrichment.csv"),
  "Reactome" = file.path(input_dir, "CF_major_stringent_Reactome_enrichment.csv")
)

qvalue_threshold <- 0.05

# Verify input files exist
missing_files <- names(files)[!sapply(files, file.exists)]
if (length(missing_files) > 0) {
  stop("Missing input files for: ", paste(missing_files, collapse = ", "),
       "\nRun run_all_analyses.R first to generate enrichment results.")
}

# Create output directory if needed
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ------------------------------------------------------------------------------
# Count significant pathways
# ------------------------------------------------------------------------------

counts <- sapply(names(files), function(db) {
  df <- read.csv(files[[db]], stringsAsFactors = FALSE)
  sum(df$qvalue < qvalue_threshold, na.rm = TRUE)
})

# Create data frame for plotting
plot_data <- data.frame(
  Database = names(counts),
  Count = as.numeric(counts),
  stringsAsFactors = FALSE
)

# Order by count (descending)
plot_data <- plot_data[order(-plot_data$Count), ]
plot_data$Database <- factor(plot_data$Database, levels = plot_data$Database)

# Print summary
cat("\n")
cat("================================================================================\n")
cat("              PATHWAY ENRICHMENT COUNTS                                         \n")
cat("================================================================================\n\n")

cat("Significant pathways (q <", qvalue_threshold, ") per database:\n\n")
print(plot_data, row.names = FALSE)
cat("\nTotal:", sum(plot_data$Count), "\n")

# ------------------------------------------------------------------------------
# Create barplot
# ------------------------------------------------------------------------------

# Define colors
bar_colors <- c(
  "GO_BP" = "#3498db",    # Blue
  "GO_CC" = "#2ecc71",    # Green
  "GO_MF" = "#9b59b6",    # Purple
  "KEGG" = "#e74c3c",     # Red
  "Reactome" = "#f39c12"  # Orange
)

# Function to create the plot
create_plot <- function() {
  par(mar = c(5, 5, 4, 2))
  
  bp <- barplot(
    plot_data$Count,
    names.arg = plot_data$Database,
    col = bar_colors[as.character(plot_data$Database)],
    border = NA,
    las = 1,
    ylim = c(0, max(plot_data$Count) * 1.15),
    ylab = "Number of Significant Pathways",
    xlab = "Pathway Database",
    main = "Pathways Enriched in Major miRNA Targets",
    cex.names = 1.1,
    cex.lab = 1.2,
    cex.main = 1.3
  )
  
  # Add count labels above bars
  text(
    bp, 
    plot_data$Count + max(plot_data$Count) * 0.03, 
    labels = plot_data$Count,
    cex = 1.0,
    font = 2
  )
  
  # Add subtitle with threshold
  mtext(
    paste0("q < ", qvalue_threshold),
    side = 3, 
    line = 0.3, 
    cex = 0.9,
    col = "gray40"
  )
}

# ------------------------------------------------------------------------------
# Save outputs
# ------------------------------------------------------------------------------

# Save to PDF (publication quality)
pdf_file <- file.path(output_dir, "pathway_counts_barplot.pdf")
pdf(pdf_file, width = 8, height = 6)
create_plot()
dev.off()

# Save to PNG (for Google Docs / presentations)
png_file <- file.path(output_dir, "pathway_counts_barplot.png")
png(png_file, width = 8, height = 6, units = "in", res = 300)
create_plot()
dev.off()

cat("\nOutputs saved:\n")
cat("  ✓", pdf_file, "\n")
cat("  ✓", png_file, "\n")
cat("\n================================================================================\n")
