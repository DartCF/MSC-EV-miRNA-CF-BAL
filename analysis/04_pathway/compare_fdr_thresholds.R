# ==============================================================================
# FDR Threshold Comparison
# ==============================================================================
# Compare how many pathways pass different q-value thresholds

compare_fdr_thresholds <- function(results) {
  
  cat("\n================================================================================\n")
  cat("              FDR THRESHOLD COMPARISON (q-values)                              \n")
  cat("================================================================================\n\n")
  
  cat("WHAT IS FDR?\n")
  cat("------------\n")
  cat("False Discovery Rate (FDR) = Expected proportion of false positives\n")
  cat("among all discoveries (significant results)\n\n")
  
  cat("q < 0.05 means: ≤5% of significant results expected to be false positives\n")
  cat("q < 0.20 means: ≤20% of significant results expected to be false positives\n\n")
  
  enrich_results <- results$enrichment$results
  thresholds <- c(0.01, 0.05, 0.10, 0.20)
  
  cat("PATHWAYS PASSING DIFFERENT FDR THRESHOLDS:\n")
  cat("==========================================\n\n")
  
  for (db_name in names(enrich_results)) {
    db_result <- enrich_results[[db_name]]
    if (is.null(db_result) || nrow(db_result@result) == 0) next
    
    cat(sprintf("%s:\n", db_name))
    cat(sprintf("  Total pathways tested: %d\n", nrow(db_result@result)))
    
    for (thresh in thresholds) {
      n_sig <- sum(db_result@result$qvalue < thresh, na.rm = TRUE)
      pct <- (n_sig / nrow(db_result@result)) * 100
      cat(sprintf("  q < %.2f (FDR < %2.0f%%): %4d pathways (%.1f%% of tested)\n", 
                  thresh, thresh*100, n_sig, pct))
    }
    cat("\n")
  }
  
  # Summary table
  cat("SUMMARY ACROSS ALL DATABASES:\n")
  cat("=============================\n\n")
  
  summary_df <- data.frame(
    Database = character(),
    FDR_1pct = integer(),
    FDR_5pct = integer(),
    FDR_10pct = integer(),
    FDR_20pct = integer(),
    stringsAsFactors = FALSE
  )
  
  for (db_name in names(enrich_results)) {
    db_result <- enrich_results[[db_name]]
    if (is.null(db_result) || nrow(db_result@result) == 0) next
    
    summary_df <- rbind(summary_df, data.frame(
      Database = db_name,
      FDR_1pct = sum(db_result@result$qvalue < 0.01, na.rm = TRUE),
      FDR_5pct = sum(db_result@result$qvalue < 0.05, na.rm = TRUE),
      FDR_10pct = sum(db_result@result$qvalue < 0.10, na.rm = TRUE),
      FDR_20pct = sum(db_result@result$qvalue < 0.20, na.rm = TRUE)
    ))
  }
  
  print(summary_df, row.names = FALSE)
  
  cat("\n\nRECOMMENDATIONS FOR YOUR MANUSCRIPT:\n")
  cat("====================================\n\n")
  
  cat("PRIMARY ANALYSIS (what you're doing now):\n")
  cat("  • Use q < 0.20 for discovery and hypothesis generation\n")
  cat("  • This is standard and acceptable for exploratory pathway analysis\n")
  cat("  • Especially appropriate when using use_all_targets mode\n\n")
  
  cat("FOR MANUSCRIPT TABLES/FIGURES:\n")
  cat("  • Report all significant pathways (q < 0.20)\n")
  cat("  • Highlight those also passing q < 0.05 in bold or with asterisks\n")
  cat("  • Example: \"** q < 0.05, * q < 0.20\"\n\n")
  
  cat("IN METHODS SECTION:\n")
  cat("  • State: \"Pathways with FDR < 0.20 (q < 0.20) were considered\n")
  cat("    significant, consistent with standards for exploratory\n")
  cat("    enrichment analysis [cite: Benjamini & Hochberg, 1995;\n")
  cat("    Storey & Tibshirani, 2003]\"\n\n")
  
  cat("FOR VALIDATION:\n")
  cat("  • Focus experimental validation on pathways with q < 0.05\n")
  cat("  • Use q < 0.20 pathways for broader context and discussion\n\n")
  
  # Calculate what the more stringent threshold would eliminate
  total_q20 <- sum(summary_df$FDR_20pct)
  total_q05 <- sum(summary_df$FDR_5pct)
  diff <- total_q20 - total_q05
  
  if (diff > 0) {
    cat("IMPACT OF USING q < 0.05 INSTEAD:\n")
    cat("==================================\n")
    cat(sprintf("  Current (q < 0.20): %d total pathways\n", total_q20))
    cat(sprintf("  If using (q < 0.05): %d total pathways\n", total_q05))
    cat(sprintf("  Difference: %d pathways would be excluded (%.1f%%)\n\n", 
                diff, (diff/total_q20)*100))
    
    if (diff / total_q20 > 0.5) {
      cat("  ⚠ You would lose >50% of your discoveries with q < 0.05\n")
      cat("    q < 0.20 is appropriate for this exploratory analysis\n")
    } else if (diff / total_q20 > 0.2) {
      cat("  ✓ Most pathways still pass q < 0.05, but q < 0.20 gives\n")
      cat("    additional discoveries worth exploring\n")
    } else {
      cat("  ✓ Most of your pathways are very strong (q < 0.05)\n")
      cat("    You could report both thresholds\n")
    }
  }
  
  cat("\n\nCITATIONS FOR FDR < 0.20:\n")
  cat("=========================\n")
  cat("Many high-impact papers use FDR < 0.20 for pathway enrichment:\n\n")
  cat("• Subramanian et al. (2005) GSEA paper used FDR < 0.25\n")
  cat("  PNAS 102(43):15545-15550\n\n")
  cat("• Huang et al. (2009) Bioinformatics reviews often cite FDR < 0.20\n")
  cat("  Nat Protoc 4(1):44-57\n\n")
  cat("• Many miRNA-pathway papers use FDR < 0.20 as their threshold\n\n")
  
  invisible(summary_df)
}

cat("\nFDR threshold comparison script loaded.\n")
cat("Usage: compare_fdr_thresholds(results)\n\n")
