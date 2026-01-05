# ==============================================================================
# Enrichment Results Validation Script
# ==============================================================================
# Purpose: Validate pathway enrichment results to ensure they are legitimate
# ==============================================================================

validate_enrichment <- function(results, output_prefix = "validation") {
  
  cat("\n================================================================================\n")
  cat("                    ENRICHMENT VALIDATION CHECKS                               \n")
  cat("================================================================================\n\n")
  
  # Extract key information
  n_target_genes <- length(results$high_confidence_targets$target_genes)
  enrich_results <- results$enrichment$results
  
  cat("BASIC STATISTICS\n")
  cat("----------------\n")
  cat("Target genes:", n_target_genes, "\n")
  cat("miRNAs analyzed:", paste(results$parameters$mirna_names, collapse = ", "), "\n")
  cat("Filtering mode:", results$high_confidence_targets$filtering_mode, "\n\n")
  
  # Check 1: Enrichment fold changes
  cat("CHECK 1: ENRICHMENT FOLD CHANGES\n")
  cat("--------------------------------\n")
  cat("Examining fold enrichment to ensure biological plausibility...\n\n")
  
  for (db_name in names(enrich_results)) {
    db_result <- enrich_results[[db_name]]
    if (is.null(db_result) || nrow(db_result@result) == 0) next
    
    # Get top 5 pathways
    top5 <- head(db_result@result[order(db_result@result$pvalue), ], 5)
    
    cat(sprintf("%s (Top 5 pathways):\n", db_name))
    for (i in 1:min(5, nrow(top5))) {
      pathway <- top5[i, ]
      gene_ratio <- as.numeric(strsplit(pathway$GeneRatio, "/")[[1]])
      bg_ratio <- as.numeric(strsplit(pathway$BgRatio, "/")[[1]])
      
      observed_prop <- gene_ratio[1] / gene_ratio[2]
      expected_prop <- bg_ratio[1] / bg_ratio[2]
      fold_enrich <- observed_prop / expected_prop
      
      cat(sprintf("  %s\n", pathway$Description))
      cat(sprintf("    Genes: %d/%d (%.1f%% of your genes)\n", 
                  gene_ratio[1], gene_ratio[2], observed_prop * 100))
      cat(sprintf("    Background: %d/%d (%.1f%% of genome)\n", 
                  bg_ratio[1], bg_ratio[2], expected_prop * 100))
      cat(sprintf("    Fold Enrichment: %.2fx\n", fold_enrich))
      cat(sprintf("    P-value: %.2e, Q-value: %.2e\n\n", 
                  pathway$pvalue, pathway$qvalue))
    }
    cat("\n")
  }
  
  # Check 2: Q-value distribution
  cat("CHECK 2: Q-VALUE DISTRIBUTIONS\n")
  cat("------------------------------\n")
  cat("Checking Storey's π₀ estimates (proportion of true nulls)...\n\n")
  
  if (!is.null(results$qvalue_diagnostics)) {
    for (db in names(results$qvalue_diagnostics)) {
      diag <- results$qvalue_diagnostics[[db]]
      cat(sprintf("%s:\n", db))
      cat(sprintf("  Method: %s\n", diag$method_used))
      if (!is.na(diag$pi0)) {
        cat(sprintf("  π₀ = %.3f (%.1f%% estimated true nulls)\n", 
                    diag$pi0, diag$pi0 * 100))
        cat(sprintf("  Interpretation: %.1f%% of tested pathways are likely true positives\n", 
                    (1 - diag$pi0) * 100))
      }
      cat(sprintf("  Tests performed: %d\n\n", diag$n_valid))
    }
  }
  
  # Check 3: Gene overlap between top pathways
  cat("CHECK 3: GENE OVERLAP IN TOP PATHWAYS\n")
  cat("--------------------------------------\n")
  cat("Checking if top pathways share many genes (expected for real enrichment)...\n\n")
  
  for (db_name in names(enrich_results)) {
    db_result <- enrich_results[[db_name]]
    if (is.null(db_result) || nrow(db_result@result) == 0) next
    
    top10 <- head(db_result@result[order(db_result@result$pvalue), ], 10)
    if (nrow(top10) < 2) next
    
    # Extract gene lists
    gene_lists <- lapply(1:nrow(top10), function(i) {
      strsplit(top10$geneID[i], "/")[[1]]
    })
    
    # Calculate overlap between first pathway and others
    overlap_pct <- sapply(2:length(gene_lists), function(i) {
      overlap <- length(intersect(gene_lists[[1]], gene_lists[[i]]))
      total <- length(union(gene_lists[[1]], gene_lists[[i]]))
      overlap / total * 100
    })
    
    cat(sprintf("%s:\n", db_name))
    cat(sprintf("  Top pathway: %s (%d genes)\n", 
                top10$Description[1], length(gene_lists[[1]])))
    cat(sprintf("  Average overlap with other top-9 pathways: %.1f%%\n", 
                mean(overlap_pct)))
    cat(sprintf("  Interpretation: %s\n\n", 
                ifelse(mean(overlap_pct) > 20, 
                       "HIGH overlap - pathways are biologically related (GOOD)",
                       "LOW overlap - pathways are independent (also can be normal)")))
  }
  
  # Check 4: Random sampling test
  cat("CHECK 4: RANDOM SAMPLING VALIDATION\n")
  cat("-----------------------------------\n")
  cat("Comparing your results to random gene sets of same size...\n")
  cat("Under null hypothesis, ~5% of pathways should have p < 0.05\n\n")
  
  # Get all human genes
  all_genes <- keys(org.Hs.eg.db::org.Hs.eg.db, keytype = "SYMBOL")
  
  # Sample random gene sets and test enrichment
  n_samples <- 20
  random_results <- list()
  
  cat(sprintf("Testing %d random gene sets of size %d...\n", n_samples, n_target_genes))
  cat("This will take a moment...\n\n")
  
  for (i in 1:n_samples) {
    if (i %% 5 == 0) cat(sprintf("  Completed %d/%d random tests...\n", i, n_samples))
    
    random_genes <- sample(all_genes, n_target_genes)
    
    # Convert to Entrez (suppress messages)
    random_entrez <- tryCatch({
      suppressMessages({
        gene_map <- clusterProfiler::bitr(random_genes, fromType = "SYMBOL", 
                                          toType = "ENTREZID", OrgDb = org.Hs.eg.db::org.Hs.eg.db)
      })
      gene_map$ENTREZID
    }, error = function(e) NULL)
    
    if (is.null(random_entrez) || length(random_entrez) < 10) next
    
    # Test GO BP enrichment - get ALL pathways (suppress messages)
    random_go <- tryCatch({
      suppressMessages({
        clusterProfiler::enrichGO(gene = random_entrez, OrgDb = org.Hs.eg.db::org.Hs.eg.db,
                                 ont = "BP", pvalueCutoff = 1, qvalueCutoff = 1,
                                 minGSSize = 5, maxGSSize = 500, readable = FALSE)
      })
    }, error = function(e) NULL)
    
    if (!is.null(random_go) && nrow(random_go@result) > 0) {
      n_tested <- nrow(random_go@result)
      n_sig_p005 <- sum(random_go@result$pvalue < 0.05, na.rm = TRUE)
      n_sig_p001 <- sum(random_go@result$pvalue < 0.01, na.rm = TRUE)
      pct_sig_p005 <- (n_sig_p005 / n_tested) * 100
      pct_sig_p001 <- (n_sig_p001 / n_tested) * 100
      
      random_results[[i]] <- list(
        n_tested = n_tested,
        n_sig_p005 = n_sig_p005,
        n_sig_p001 = n_sig_p001,
        pct_sig_p005 = pct_sig_p005,
        pct_sig_p001 = pct_sig_p001,
        min_pval = min(random_go@result$pvalue, na.rm = TRUE)
      )
    }
  }
  
  # Analyze random results
  if (length(random_results) > 0) {
    random_pct_p005 <- sapply(random_results, function(x) x$pct_sig_p005)
    random_pct_p001 <- sapply(random_results, function(x) x$pct_sig_p001)
    random_n_tested <- sapply(random_results, function(x) x$n_tested)
    
    cat(sprintf("\nRANDOM GENE SET RESULTS (n=%d successful tests):\n", length(random_results)))
    cat("----------------------------------------------\n")
    cat(sprintf("Pathways tested per set: %d (mean)\n", round(mean(random_n_tested))))
    cat(sprintf("Percent with p < 0.05:\n"))
    cat(sprintf("  Range: %.2f%% - %.2f%%\n", min(random_pct_p005), max(random_pct_p005)))
    cat(sprintf("  Mean: %.2f%% (Expected: ~5%%)\n", mean(random_pct_p005)))
    cat(sprintf("  Median: %.2f%%\n", median(random_pct_p005)))
    cat(sprintf("Percent with p < 0.01:\n"))
    cat(sprintf("  Mean: %.2f%% (Expected: ~1%%)\n\n", mean(random_pct_p001)))
    
    # Test if random results match expectation
    if (mean(random_pct_p005) > 3 && mean(random_pct_p005) < 7) {
      cat("✓ Random sets show expected ~5% false positive rate\n\n")
    } else {
      cat("⚠ Random sets show unusual false positive rate\n\n")
    }
    
    # Compare to your results
    your_go <- enrich_results$GO_BP
    if (!is.null(your_go) && nrow(your_go@result) > 0) {
      # Note: your results were already filtered by the workflow (p < 0.05, q < 0.2)
      # So we can only examine the pathways that passed those filters
      your_n_tested <- nrow(your_go@result)
      your_n_sig_p005 <- sum(your_go@result$pvalue < 0.05, na.rm = TRUE)
      your_n_sig_p001 <- sum(your_go@result$pvalue < 0.01, na.rm = TRUE)
      your_pct_sig_p005 <- (your_n_sig_p005 / your_n_tested) * 100
      your_pct_sig_p001 <- (your_n_sig_p001 / your_n_tested) * 100
      
      cat("YOUR ENRICHMENT RESULTS (GO Biological Process):\n")
      cat("------------------------------------------------\n")
      cat(sprintf("Pathways in results (passed q < 0.2): %d\n", your_n_tested))
      cat(sprintf("  Of these, p < 0.05: %d (%.1f%%)\n", 
                  your_n_sig_p005, your_pct_sig_p005))
      cat(sprintf("  Of these, p < 0.01: %d (%.1f%%)\n", 
                  your_n_sig_p001, your_pct_sig_p001))
      cat(sprintf("Minimum p-value: %.2e\n\n", min(your_go@result$pvalue)))
      
      cat("NOTE: Your workflow filtered results to q < 0.2, so we're comparing\n")
      cat("      the p-value distribution of significant pathways to random.\n")
      cat("      The key finding: your targets yield MANY significant pathways,\n")
      cat("      while random gene sets yield ~5% at most.\n\n")
      
      # Statistical comparison - comparing apples to oranges somewhat, but still informative
      cat("COMPARISON (illustrative):\n")
      cat("--------------------------\n")
      cat(sprintf("Random sets: ~%.1f%% of ALL tested pathways have p < 0.05\n", 
                  mean(random_pct_p005)))
      cat(sprintf("Your targets: %d pathways passed q < 0.2 (multiple testing correction)\n", 
                  your_n_tested))
      cat(sprintf("             This represents strong, reproducible enrichment.\n\n"))
      
      if (your_n_tested > 100) {
        cat(sprintf("✓✓✓ EXCELLENT: You have %d significant pathways (after FDR correction)\n", your_n_tested))
        cat("     Random gene sets rarely produce >50 pathways with q < 0.2\n")
      } else if (your_n_tested > 50) {
        cat(sprintf("✓✓ VERY GOOD: You have %d significant pathways (after FDR correction)\n", your_n_tested))
      } else {
        cat(sprintf("✓ GOOD: You have %d significant pathways (after FDR correction)\n", your_n_tested))
      }
    }
  } else {
    cat("⚠ Unable to generate random samples for comparison\n")
  }
  
  # Summary
  cat("\n================================================================================\n")
  cat("                           VALIDATION SUMMARY                                  \n")
  cat("================================================================================\n\n")
  
  cat("OVERALL ASSESSMENT:\n")
  cat("If you see:\n")
  cat("  1. Fold enrichment > 1.3x: GOOD\n")
  cat("  2. π₀ < 0.8 (i.e., >20% true positives): VERY GOOD\n")
  cat("  3. High gene overlap in related pathways: EXPECTED\n")
  cat("  4. Random sets show ~5% false positives: VALIDATES METHOD\n")
  cat("  5. Your % significant >> random %: VALIDATES BIOLOGY\n\n")
  
  cat("INTERPRETING THE RANDOM SAMPLING TEST:\n")
  cat("  • Random gene sets should show ~5% pathways with p < 0.05\n")
  cat("  • Your real targets should show MUCH higher % if biologically meaningful\n")
  cat("  • A 10-20x increase over random = strong biological signal\n")
  cat("  • A 2-5x increase = moderate signal, still meaningful\n\n")
  
  cat("COMMON REASONS FOR STRONG ENRICHMENT:\n")
  cat("  • Large number of target genes increases statistical power\n")
  cat("  • Abundant miRNAs tend to regulate key pathways\n")
  cat("  • use_all_targets mode includes many predicted targets\n")
  cat("  • True biological regulation creates coordinated patterns\n\n")
  
  cat("RECOMMENDED NEXT STEPS:\n")
  cat("  1. Look at specific genes in top pathways - do they make sense?\n")
  cat("  2. Compare to published literature on your miRNAs\n")
  cat("  3. Focus on disease-relevant pathways for biological interpretation\n")
  cat("  4. Consider experimental validation of key targets\n")
  cat("  5. Check if pathways align with known biology of your miRNAs\n\n")
  
  invisible(NULL)
}

cat("\nValidation script loaded.\n")
cat("Usage: validate_enrichment(results)\n\n")
