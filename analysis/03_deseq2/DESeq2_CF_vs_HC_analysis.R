# ==============================================================================
# DESeq2 Analysis: CF BAL vs HC BAL Exposure
# ==============================================================================
#
# PURPOSE:
# Perform conventional differential expression analysis to identify miRNAs
# significantly changed between CF and HC BAL-exposed BM-hMSC EVs, then
# contextualize these findings by abundance tier.
#
# INPUTS:
#   - data/raw/GSE282919_VLP00436_miRNA_P33_21NOV2019_QC_and_Processed_Data.xlsx
#   - data/metadata/sample_metadata.csv
#   - results/02_minor_miRNA/differential_expression_CF_vs_HC.csv (optional, for comparison)
#
# OUTPUTS:
#   - results/03_deseq2/DESeq2_results_CF_vs_HC.csv
#   - results/03_deseq2/DESeq2_DE_miRNAs_CF_vs_HC.csv
#   - results/03_deseq2/DE_rate_by_abundance_tier.csv
#   - results/03_deseq2/DESeq2_vs_compositional_comparison.csv
#   - results/03_deseq2/volcano_CF_vs_HC_by_abundance.png/.pdf
#   - results/03_deseq2/volcano_faceted_by_abundance.png/.pdf
#
# ==============================================================================

# Load required libraries
library(here)
library(readxl)
library(readr)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(ggrepel)

# Create output directory
output_dir <- here("results", "03_deseq2")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Helper function to save figures
save_figure <- function(plot, filename, width = 10, height = 8) {
  ggsave(file.path(output_dir, paste0(filename, ".png")), plot, 
         width = width, height = height, dpi = 300)
  ggsave(file.path(output_dir, paste0(filename, ".pdf")), plot, 
         width = width, height = height)
  cat("Saved:", filename, ".png/.pdf\n")
}

# ==============================================================================
# 1. DATA LOADING
# ==============================================================================

# Define input file paths
raw_data_file <- here("data", "raw", 
                      "GSE282919_VLP00436_miRNA_P33_21NOV2019_QC_and_Processed_Data.xlsx")
metadata_file <- here("data", "metadata", "sample_metadata.csv")

# Load sample metadata
sample_metadata <- read_csv(metadata_file, show_col_types = FALSE)

# Create sample mapping for compatibility
sample_mapping <- sample_metadata %>%
  dplyr::select(Sample, Treatment_Group) %>%
  as.data.frame()
rownames(sample_mapping) <- sample_mapping$Sample
names(sample_mapping)[names(sample_mapping) == "Treatment_Group"] <- "comparison_7"

# Read the data from GEO Excel file (skip first 10 rows, data starts at row 11)
raw_data <- read_excel(raw_data_file, sheet = "Raw", skip = 10, col_names = FALSE)

# Set column names: first column is Gene_ID, rest are sample columns
n_samples_file <- ncol(raw_data) - 1
colnames(raw_data) <- c("Gene_ID", paste0("sample_", 1:n_samples_file))

# Filter out control (CTRL) and housekeeping (HK) probes
raw_filtered <- raw_data %>%
  filter(!grepl("^CTRL|^HK", Gene_ID))

cat("Probes removed (CTRL/HK):", nrow(raw_data) - nrow(raw_filtered), "\n")
cat("Probes retained:", nrow(raw_filtered), "\n\n")

# ==============================================================================
# 2. SUBSET TO CF vs HC COMPARISON
# ==============================================================================

# Define groups for this comparison
comparison_groups <- c("CF_Asp_Neg", "HC")
samples_cf_hc <- sample_metadata %>%
  filter(Treatment_Group %in% comparison_groups, In_Current_Study == TRUE) %>%
  pull(Sample)

cat("=== Sample Summary ===\n")
cat("CF_Asp_Neg samples:", sum(sample_mapping[samples_cf_hc, "comparison_7"] == "CF_Asp_Neg"), "\n")
cat("HC samples:", sum(sample_mapping[samples_cf_hc, "comparison_7"] == "HC"), "\n")
cat("Total for comparison:", length(samples_cf_hc), "\n\n")

# Extract count data for CF vs HC
count_data <- raw_filtered[, c("Gene_ID", samples_cf_hc)]
rownames(count_data) <- count_data$Gene_ID
count_data <- count_data[, samples_cf_hc]

# Ensure integer counts for DESeq2
count_matrix <- as.matrix(count_data)
count_matrix <- apply(count_matrix, 2, as.numeric)
rownames(count_matrix) <- raw_filtered$Gene_ID
mode(count_matrix) <- "integer"

# Create colData for DESeq2
col_data <- data.frame(
  row.names = samples_cf_hc,
  condition = factor(sample_mapping[samples_cf_hc, "comparison_7"],
                     levels = c("HC", "CF_Asp_Neg"))  # HC is reference
)

cat("=== Condition Summary ===\n")
print(table(col_data$condition))
cat("\n")

# ==============================================================================
# 3. CALCULATE FRACTIONAL ABUNDANCE (for later categorization)
# ==============================================================================

# Column sums (total reads per sample)
col_totals <- colSums(count_matrix)

# Fractional abundance matrix
frac_matrix <- sweep(count_matrix, 2, col_totals, FUN = "/")

# Mean fractional abundance by condition
cf_samples <- rownames(col_data)[col_data$condition == "CF_Asp_Neg"]
hc_samples <- rownames(col_data)[col_data$condition == "HC"]

mean_frac_cf <- rowMeans(frac_matrix[, cf_samples])
mean_frac_hc <- rowMeans(frac_matrix[, hc_samples])
mean_frac_overall <- rowMeans(frac_matrix)

# Create abundance summary
abundance_summary <- data.frame(
  miRNA = rownames(count_matrix),
  mean_frac_CF = mean_frac_cf,
  mean_frac_HC = mean_frac_hc,
  mean_frac_overall = mean_frac_overall,
  mean_pct_overall = mean_frac_overall * 100
)

# Calculate RPM (reads per million)
rpm_matrix <- frac_matrix * 1e6
mean_rpm_cf <- rowMeans(rpm_matrix[, cf_samples])
mean_rpm_hc <- rowMeans(rpm_matrix[, hc_samples])
mean_rpm_overall <- rowMeans(rpm_matrix)

abundance_summary$mean_RPM_CF <- mean_rpm_cf
abundance_summary$mean_RPM_HC <- mean_rpm_hc
abundance_summary$mean_RPM_overall <- mean_rpm_overall

# Flag miRNAs below functional threshold (< 1 RPM in CF samples)
abundance_summary$below_1RPM_CF <- abundance_summary$mean_RPM_CF < 1

cat("=== RPM Distribution in CF Samples ===\n")
cat(sprintf("miRNAs with mean RPM < 1 in CF:    %d (%.1f%%)\n",
            sum(abundance_summary$below_1RPM_CF),
            100 * mean(abundance_summary$below_1RPM_CF)))
cat(sprintf("miRNAs with mean RPM >= 1 in CF:   %d (%.1f%%)\n",
            sum(!abundance_summary$below_1RPM_CF),
            100 * mean(!abundance_summary$below_1RPM_CF)))
cat("\n")

# Categorize by abundance tier
abundance_summary$abundance_tier <- cut(
  abundance_summary$mean_pct_overall,
  breaks = c(-Inf, 0.01, 0.1, 1, Inf),
  labels = c("Trace (<0.01%)", "Low (0.01-0.1%)", "Minor (0.1-1%)", "Major (≥1%)"),
  right = FALSE
)

cat("=== Abundance Tier Distribution ===\n")
print(table(abundance_summary$abundance_tier))
cat("\n")

# ==============================================================================
# 4. DESeq2 ANALYSIS
# ==============================================================================

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = col_data,
  design = ~ condition
)

# Filter low-count miRNAs (at least 10 reads total across samples)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
cat("miRNAs retained after filtering:", sum(keep), "of", length(keep), "\n\n")

# Run DESeq2
dds <- DESeq(dds)

# Extract results (CF vs HC, positive log2FC means higher in CF)
res <- results(dds, contrast = c("condition", "CF_Asp_Neg", "HC"))
res_df <- as.data.frame(res)
res_df$miRNA <- rownames(res_df)

# Merge with abundance data
res_df <- merge(res_df, abundance_summary, by = "miRNA", all.x = TRUE)

# ==============================================================================
# 5. RESULTS SUMMARY
# ==============================================================================

# Define significance threshold
fdr_threshold <- 0.05

# Count DE miRNAs
res_df$is_DE <- !is.na(res_df$padj) & res_df$padj < fdr_threshold
n_de <- sum(res_df$is_DE, na.rm = TRUE)

# Flag DE miRNAs below functional threshold
res_df$DE_below_1RPM_CF <- res_df$is_DE & res_df$below_1RPM_CF

cat("=======================================================\n")
cat("DESeq2 RESULTS SUMMARY\n")
cat("=======================================================\n\n")

cat("Total miRNAs tested:", nrow(res_df), "\n")
cat("Differentially expressed (FDR <", fdr_threshold, "):", n_de, "\n\n")

# DE by direction
de_mirnas <- res_df[res_df$is_DE, ]
n_up <- sum(de_mirnas$log2FoldChange > 0, na.rm = TRUE)
n_down <- sum(de_mirnas$log2FoldChange < 0, na.rm = TRUE)
cat("  - Upregulated in CF:", n_up, "\n")
cat("  - Downregulated in CF:", n_down, "\n\n")

# DE by abundance tier
cat("=== DE miRNAs by Abundance Tier ===\n")
de_by_tier <- table(de_mirnas$abundance_tier)
print(de_by_tier)
cat("\n")

# Calculate percentage of DE miRNAs in each tier
tier_totals <- table(res_df$abundance_tier)
cat("=== Percentage DE within each Abundance Tier ===\n")
for (tier in names(tier_totals)) {
  n_tier_total <- tier_totals[tier]
  n_tier_de <- ifelse(tier %in% names(de_by_tier), de_by_tier[tier], 0)
  pct <- round(100 * n_tier_de / n_tier_total, 1)
  cat(sprintf("  %s: %d of %d (%.1f%%)\n", tier, n_tier_de, n_tier_total, pct))
}
cat("\n")

# DE miRNAs below functional threshold
cat("=== DE miRNAs Below Functional Threshold ===\n")
n_de_below_1rpm <- sum(de_mirnas$below_1RPM_CF, na.rm = TRUE)
n_de_above_1rpm <- sum(!de_mirnas$below_1RPM_CF, na.rm = TRUE)
cat(sprintf("DE miRNAs with mean RPM < 1 in CF:  %d of %d (%.1f%%)\n",
            n_de_below_1rpm, n_de, 100 * n_de_below_1rpm / n_de))
cat(sprintf("DE miRNAs with mean RPM >= 1 in CF: %d of %d (%.1f%%)\n",
            n_de_above_1rpm, n_de, 100 * n_de_above_1rpm / n_de))
cat("\n")

# ==============================================================================
# 6. STATISTICAL TEST: DE RATE BY ABUNDANCE TIER
# ==============================================================================

cat("=======================================================\n")
cat("STATISTICAL TEST: DE RATE VARIES BY ABUNDANCE TIER\n")
cat("=======================================================\n\n")

# Create contingency table
tier_order <- c("Trace (<0.01%)", "Low (0.01-0.1%)", "Minor (0.1-1%)", "Major (≥1%)")
contingency_data <- data.frame(
  Tier = tier_order,
  Total = as.numeric(tier_totals[tier_order]),
  DE = sapply(tier_order, function(t) ifelse(t %in% names(de_by_tier), de_by_tier[t], 0)),
  stringsAsFactors = FALSE
)
contingency_data$NonDE <- contingency_data$Total - contingency_data$DE
contingency_data$Pct_DE <- round(100 * contingency_data$DE / contingency_data$Total, 1)

# Chi-square test
contingency_matrix <- as.matrix(contingency_data[, c("DE", "NonDE")])
rownames(contingency_matrix) <- contingency_data$Tier
chisq_result <- chisq.test(contingency_matrix)

cat("Contingency Table:\n")
print(contingency_data[, c("Tier", "Total", "DE", "NonDE", "Pct_DE")])
cat("\n")

cat("Chi-square test for independence:\n")
cat(sprintf("  X-squared = %.2f, df = %d, p-value = %.2e\n\n", 
            chisq_result$statistic, chisq_result$parameter, chisq_result$p.value))

# Pairwise Fisher's exact tests
cat("Pairwise Fisher's exact tests (each tier vs Major):\n")
major_de <- contingency_data$DE[contingency_data$Tier == "Major (≥1%)"]
major_non_de <- contingency_data$NonDE[contingency_data$Tier == "Major (≥1%)"]

for (tier in c("Trace (<0.01%)", "Low (0.01-0.1%)", "Minor (0.1-1%)")) {
  tier_de <- contingency_data$DE[contingency_data$Tier == tier]
  tier_non_de <- contingency_data$NonDE[contingency_data$Tier == tier]
  
  fisher_matrix <- matrix(c(tier_de, tier_non_de, major_de, major_non_de), 
                          nrow = 2, byrow = TRUE)
  fisher_result <- fisher.test(fisher_matrix)
  
  cat(sprintf("  %s vs Major: OR = %.2f, p = %.4f\n", 
              tier, fisher_result$estimate, fisher_result$p.value))
}
cat("\n")

# ==============================================================================
# 7. CUMULATIVE ABUNDANCE OF DE miRNAs
# ==============================================================================

cat("=======================================================\n")
cat("CUMULATIVE ABUNDANCE OF DE miRNAs\n")
cat("=======================================================\n\n")

total_cargo_de <- sum(de_mirnas$mean_pct_overall, na.rm = TRUE)
cat(sprintf("Total cargo from DE miRNAs: %.2f%%\n", total_cargo_de))

de_up <- de_mirnas[de_mirnas$log2FoldChange > 0, ]
de_down <- de_mirnas[de_mirnas$log2FoldChange < 0, ]
cat(sprintf("  - Upregulated in CF: %.2f%%\n", sum(de_up$mean_pct_overall, na.rm = TRUE)))
cat(sprintf("  - Downregulated in CF: %.2f%%\n", sum(de_down$mean_pct_overall, na.rm = TRUE)))
cat("\n")

non_de <- res_df[!res_df$is_DE | is.na(res_df$is_DE), ]
total_cargo_non_de <- sum(non_de$mean_pct_overall, na.rm = TRUE)
cat(sprintf("Total cargo from non-DE miRNAs: %.2f%%\n", total_cargo_non_de))
cat("\n")

# ==============================================================================
# 8. VOLCANO PLOT WITH ABUNDANCE COLORING
# ==============================================================================

plot_data <- res_df %>%
  filter(!is.na(log2FoldChange) & !is.na(pvalue)) %>%
  mutate(
    neg_log10_pval = -log10(pvalue),
    is_significant = !is.na(padj) & padj < 0.05,
    significance = case_when(
      is_significant & log2FoldChange > 1 ~ "Up (FC>2, FDR<0.05)",
      is_significant & log2FoldChange < -1 ~ "Down (FC>2, FDR<0.05)",
      is_significant ~ "DE (FDR<0.05)",
      TRUE ~ "Not significant"
    ),
    abundance_tier = factor(abundance_tier, 
                            levels = c("Major (≥1%)", "Minor (0.1-1%)", 
                                       "Low (0.01-0.1%)", "Trace (<0.01%)"))
  )

# Create volcano plot colored by abundance tier
p_volcano <- ggplot(plot_data, aes(x = log2FoldChange, y = neg_log10_pval)) +
  geom_point(aes(color = abundance_tier), alpha = 0.7, size = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray40") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray40") +
  scale_color_manual(
    values = c("Major (≥1%)" = "#E41A1C",
               "Minor (0.1-1%)" = "#377EB8",
               "Low (0.01-0.1%)" = "#4DAF4A",
               "Trace (<0.01%)" = "gray70"),
    name = "Abundance Tier"
  ) +
  labs(
    title = "Differential Expression: CF vs HC BAL Exposure",
    subtitle = "Points colored by fractional abundance tier",
    x = expression(log[2]~"Fold Change (CF / HC)"),
    y = expression(-log[10]~"p-value")
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

save_figure(p_volcano, "volcano_CF_vs_HC_by_abundance", width = 8, height = 6)

# ==============================================================================
# 9. FACETED VOLCANO PLOT BY ABUNDANCE TIER
# ==============================================================================

p_facet <- ggplot(plot_data, aes(x = log2FoldChange, y = neg_log10_pval)) +
  geom_point(aes(color = is_significant), alpha = 0.6, size = 1.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray40") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray40") +
  scale_color_manual(
    values = c("TRUE" = "#E41A1C", "FALSE" = "gray60"),
    labels = c("TRUE" = "FDR < 0.05", "FALSE" = "Not significant"),
    name = ""
  ) +
  facet_wrap(~ abundance_tier, scales = "free_y") +
  labs(
    title = "Differential Expression by Abundance Tier",
    subtitle = sprintf("Only %d of 13 Major miRNAs (≥1%% of cargo) show significant differential expression", 
                       sum(plot_data$abundance_tier == "Major (≥1%)" & plot_data$is_significant, na.rm = TRUE)),
    x = expression(log[2]~"Fold Change (CF / HC)"),
    y = expression(-log[10]~"p-value")
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "gray90", color = NA),
    strip.text = element_text(face = "bold")
  )

# Add labels for DE Major miRNAs
de_major <- plot_data %>%
  filter(abundance_tier == "Major (≥1%)" & is_significant)

if (nrow(de_major) > 0) {
  p_facet <- p_facet +
    geom_text_repel(
      data = de_major,
      aes(label = miRNA),
      size = 3,
      box.padding = 0.5,
      max.overlaps = 20
    )
}

save_figure(p_facet, "volcano_faceted_by_abundance", width = 10, height = 8)

# ==============================================================================
# 10. STATUS OF MAJOR miRNAs
# ==============================================================================

cat("=======================================================\n")
cat("CRITICAL FINDING: STATUS OF MAJOR miRNAs\n")
cat("=======================================================\n\n")

major_mirnas <- res_df %>%
  filter(abundance_tier == "Major (≥1%)") %>%
  arrange(desc(mean_pct_overall))

n_major_significant <- sum(major_mirnas$padj < 0.05, na.rm = TRUE)
cat(sprintf("Of 13 Major miRNAs, %d show significant differential expression (FDR < 0.05):\n\n",
            n_major_significant))
print(major_mirnas %>% 
        dplyr::select(miRNA, log2FoldChange, pvalue, padj, mean_pct_overall) %>%
        mutate(significant = ifelse(padj < 0.05, "YES", "no")),
      row.names = FALSE)

cat("\n")
cat("Cumulative abundance of Major miRNAs:", 
    sprintf("%.1f%%", sum(major_mirnas$mean_pct_overall)), "\n")
cat("Major miRNAs with FDR < 0.05:", n_major_significant, "of 13\n\n")

# ==============================================================================
# 11. EXPORT RESULTS
# ==============================================================================

# Sort by adjusted p-value
res_df_sorted <- res_df %>%
  arrange(padj) %>%
  dplyr::select(miRNA, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj,
         mean_pct_overall, abundance_tier, is_DE)

# Export full results
write.csv(res_df_sorted, file.path(output_dir, "DESeq2_results_CF_vs_HC.csv"), row.names = FALSE)
cat("Full results exported: DESeq2_results_CF_vs_HC.csv\n")

# Export DE rate by tier table
tier_table <- data.frame(
  Abundance_Tier = tier_order,
  Total_miRNAs = as.numeric(tier_totals[tier_order]),
  DE_FDR_0.05 = sapply(tier_order, function(t) ifelse(t %in% names(de_by_tier), de_by_tier[t], 0)),
  stringsAsFactors = FALSE
)
tier_table$Pct_DE <- round(100 * tier_table$DE_FDR_0.05 / tier_table$Total_miRNAs, 1)
write.csv(tier_table, file.path(output_dir, "DE_rate_by_abundance_tier.csv"), row.names = FALSE)
cat("DE rate by tier exported: DE_rate_by_abundance_tier.csv\n")

# Export DE miRNAs only
write.csv(res_df_sorted[res_df_sorted$is_DE, ], 
          file.path(output_dir, "DESeq2_DE_miRNAs_CF_vs_HC.csv"), row.names = FALSE)
cat("DE miRNAs exported: DESeq2_DE_miRNAs_CF_vs_HC.csv\n\n")

# ==============================================================================
# 12. COMPARISON: DESeq2 vs COMPOSITIONAL ANALYSIS
# ==============================================================================

cat("=======================================================\n")
cat("COMPARISON: DESeq2 vs COMPOSITIONAL ANALYSIS\n")
cat("=======================================================\n\n")

# Look for compositional analysis results from minor miRNA analysis
compositional_file <- here("results", "02_minor_miRNA", "differential_expression_CF_vs_HC.csv")

if (file.exists(compositional_file)) {
  compositional_results <- read.csv(compositional_file, stringsAsFactors = FALSE)
  
  cat("Columns in compositional results:\n")
  print(colnames(compositional_results))
  cat("\n")
  
  # Identify column names
  mirna_col <- grep("miRNA|Gene|name", colnames(compositional_results), ignore.case = TRUE, value = TRUE)[1]
  log2fc_col <- grep("log2|FC|fold", colnames(compositional_results), ignore.case = TRUE, value = TRUE)[1]
  pval_col <- grep("^p\\.?val|pvalue|P_value", colnames(compositional_results), ignore.case = TRUE, value = TRUE)[1]
  
  cat(sprintf("Using columns: miRNA='%s', log2FC='%s', pvalue='%s'\n\n", 
              mirna_col, log2fc_col, pval_col))
  
  # Define DE in compositional analysis (|log2FC| > 1 and p < 0.05)
  compositional_results$DE_compositional <- 
    abs(compositional_results[[log2fc_col]]) > 1 & 
    compositional_results[[pval_col]] < 0.05
  
  cat(sprintf("Compositional analysis: %d miRNAs DE (|log2FC|>1, p<0.05)\n\n",
              sum(compositional_results$DE_compositional, na.rm = TRUE)))
  
  # Merge with DESeq2 results
  comparison_df <- merge(
    res_df[, c("miRNA", "log2FoldChange", "pvalue", "padj", "is_DE", "mean_pct_overall", "abundance_tier")],
    compositional_results[, c(mirna_col, log2fc_col, pval_col, "DE_compositional")],
    by.x = "miRNA", by.y = mirna_col,
    all.x = TRUE
  )
  
  colnames(comparison_df)[colnames(comparison_df) == log2fc_col] <- "log2FC_compositional"
  colnames(comparison_df)[colnames(comparison_df) == pval_col] <- "pvalue_compositional"
  
  comparison_df$DE_DESeq2 <- comparison_df$is_DE
  comparison_df$DE_compositional[is.na(comparison_df$DE_compositional)] <- FALSE
  
  # Concordance summary
  n_de_both <- sum(comparison_df$DE_DESeq2 & comparison_df$DE_compositional, na.rm = TRUE)
  n_de_deseq_only <- sum(comparison_df$DE_DESeq2 & !comparison_df$DE_compositional, na.rm = TRUE)
  n_de_comp_only <- sum(!comparison_df$DE_DESeq2 & comparison_df$DE_compositional, na.rm = TRUE)
  n_de_neither <- sum(!comparison_df$DE_DESeq2 & !comparison_df$DE_compositional, na.rm = TRUE)
  
  cat("Concordance between DESeq2 (FDR < 0.05) and compositional (|log2FC|>1, p<0.05):\n\n")
  cat(sprintf("  DE by both methods:           %d\n", n_de_both))
  cat(sprintf("  DE by DESeq2 only:            %d\n", n_de_deseq_only))
  cat(sprintf("  DE by compositional only:     %d\n", n_de_comp_only))
  cat(sprintf("  Not DE by either:             %d\n", n_de_neither))
  cat("\n")
  
  # Export comparison
  write.csv(comparison_df, file.path(output_dir, "DESeq2_vs_compositional_comparison.csv"), row.names = FALSE)
  cat("Comparison exported: DESeq2_vs_compositional_comparison.csv\n\n")
  
} else {
  cat("Note: Compositional analysis file not found at:", compositional_file, "\n")
  cat("Skipping comparison. Run 02_minor_miRNA analysis first if comparison is desired.\n\n")
}

# ==============================================================================
# 13. MANUSCRIPT SUMMARY
# ==============================================================================

cat("=======================================================\n")
cat("MANUSCRIPT SUMMARY\n")
cat("=======================================================\n\n")

n_major_de <- sum(de_mirnas$abundance_tier == "Major (≥1%)", na.rm = TRUE)
n_minor_de <- sum(de_mirnas$abundance_tier == "Minor (0.1-1%)", na.rm = TRUE)
n_low_de <- sum(de_mirnas$abundance_tier == "Low (0.01-0.1%)", na.rm = TRUE)
n_trace_de <- sum(de_mirnas$abundance_tier == "Trace (<0.01%)", na.rm = TRUE)

cat(sprintf('DESeq2 identified %d miRNAs as significantly changed (FDR < 0.05):
  - %d upregulated in CF
  - %d downregulated in CF

By abundance tier:
  - Major (≥1%%): %d DE

  - Minor (0.1-1%%): %d DE
  - Low (0.01-0.1%%): %d DE
  - Trace (<0.01%%): %d DE

Cumulative cargo from DE miRNAs: %.2f%%
Cumulative cargo from non-DE miRNAs: %.2f%%
\n',
n_de, n_up, n_down, 
n_major_de, n_minor_de, n_low_de, n_trace_de,
total_cargo_de, total_cargo_non_de))

cat("\n=== All outputs saved to: ===\n")
cat(output_dir, "\n")

cat("\n=== SESSION INFO ===\n")
sessionInfo()
