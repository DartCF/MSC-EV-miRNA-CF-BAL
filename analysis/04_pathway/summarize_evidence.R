# ==============================================================================
# Evidence Analysis Functions
# ==============================================================================

summarize_target_evidence <- function(results, save_to_file = FALSE) {
  hc_df <- results$high_confidence_targets$high_confidence_df
  
  if (is.null(hc_df) || nrow(hc_df) == 0) {
    cat("No targets found.\n")
    return(NULL)
  }
  
  filtering_mode <- results$high_confidence_targets$filtering_mode
  
  cat("\n================================================================================\n")
  cat("                         EVIDENCE SUMMARY                                      \n")
  cat("================================================================================\n\n")
  
  cat("OVERALL STATISTICS\n")
  cat("------------------\n")
  cat("Filtering mode:", filtering_mode, "\n")
  cat("Total target pairs:", nrow(hc_df), "\n")
  cat("Unique genes:", length(unique(hc_df$target_symbol)), "\n\n")
  
  # Evidence breakdown
  has_n_databases <- "n_databases" %in% colnames(hc_df)
  
  if (has_n_databases) {
    cat("DATABASE SUPPORT\n")
    cat("----------------\n")
    db_table <- table(hc_df$n_databases)
    for (i in names(db_table)) {
      cat(sprintf("  %s database(s): %d targets\n", i, db_table[i]))
    }
    cat("\n")
  }
  
  # Q-value diagnostics
  if (!is.null(results$qvalue_diagnostics)) {
    cat("Q-VALUE DIAGNOSTICS\n")
    cat("-------------------\n")
    for (db in names(results$qvalue_diagnostics)) {
      diag <- results$qvalue_diagnostics[[db]]
      cat(sprintf("%s:\n", db))
      cat(sprintf("  Method: %s\n", diag$method_used))
      if (!is.na(diag$pi0)) {
        cat(sprintf("  π₀ = %.3f (%.1f%% estimated true nulls)\n", diag$pi0, diag$pi0 * 100))
      }
    }
    cat("\n")
  }
  
  cat("================================================================================\n\n")
  
  # Save if requested
  if (save_to_file) {
    output_dir <- results$parameters$output_dir
    output_prefix <- results$parameters$output_prefix
    
    f <- file.path(output_dir, paste0(output_prefix, "_evidence_summary.csv"))
    write.csv(hc_df, f, row.names = FALSE)
    cat("✓ Saved to:", f, "\n\n")
  }
  
  return(list(high_confidence_data = hc_df, filtering_mode = filtering_mode))
}

cat("Evidence analysis functions loaded.\n")
cat("Usage: evidence <- summarize_target_evidence(results, save_to_file = TRUE)\n\n")
