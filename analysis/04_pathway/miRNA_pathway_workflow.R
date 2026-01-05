# ==============================================================================
# miRNA Pathway Analysis Workflow
# ==============================================================================
# Complete functional workflow with:
# - Storey's q-value for FDR control
# - Three filtering modes (including use_all_targets for less-studied miRNAs)
# - disease_keywords parameter (works with any disease)
# ==============================================================================

# Load Required Packages
required_packages <- c("multiMiR", "clusterProfiler", "org.Hs.eg.db", 
                      "ReactomePA", "enrichplot", "qvalue", 
                      "dplyr", "ggplot2", "tidyr")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste0("Package '", pkg, "' required but not installed.\n",
                "Install: BiocManager::install('", pkg, "')"))
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

cat("miRNA Pathway Analysis Workflow loaded successfully.\n\n")

# ==============================================================================
# Storey's q-value Function
# ==============================================================================

calculate_storey_qvalue <- function(pvalues, lambda = seq(0, 0.95, 0.05), 
                                   pi0.method = "smoother") {
  valid_idx <- !is.na(pvalues)
  valid_pvals <- pvalues[valid_idx]
  n_valid <- length(valid_pvals)
  
  if (n_valid < 2) {
    warning("Fewer than 2 p-values. Returning as q-values.")
    return(list(qvalues = pvalues, pi0 = 1, method_used = "none", 
               n_tests = length(pvalues), n_valid = n_valid))
  }
  
  tryCatch({
    qobj <- qvalue::qvalue(p = valid_pvals, lambda = lambda, pi0.method = pi0.method)
    qvalues_full <- rep(NA_real_, length(pvalues))
    qvalues_full[valid_idx] <- qobj$qvalues
    
    cat(sprintf("✓ Storey's q-value: π₀ = %.3f (%d tests, %.1f%% estimated true nulls)\n", 
                qobj$pi0, n_valid, qobj$pi0 * 100))
    
    return(list(qvalues = qvalues_full, pi0 = qobj$pi0, method_used = "storey",
               n_tests = length(pvalues), n_valid = n_valid, qobj = qobj))
  }, error = function(e) {
    warning("Storey's q-value failed. Using Benjamini-Hochberg.")
    return(list(qvalues = p.adjust(pvalues, "BH"), pi0 = NA, method_used = "BH",
               n_tests = length(pvalues), n_valid = n_valid))
  })
}

# ==============================================================================
# Fetch Targets from multiMiR
# ==============================================================================

fetch_mirna_targets <- function(mirna_names, skip_validation = FALSE) {
  cat("\n=== FETCHING TARGETS ===\n")
  cat("miRNAs:", paste(mirna_names, collapse = ", "), "\n\n")
  
  mirna_names_formatted <- ifelse(grepl("^hsa-", mirna_names), mirna_names, 
                                 paste0("hsa-", mirna_names))
  
  cat("Fetching predicted targets...\n")
  predicted_result <- tryCatch({
    multiMiR::get_multimir(org = "hsa", mirna = mirna_names_formatted, 
                          table = "predicted", summary = FALSE)
  }, error = function(e) {
    cat("Error:", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(predicted_result) || nrow(predicted_result@data) == 0) {
    stop("No predicted targets found")
  }
  
  predicted_targets <- predicted_result@data
  cat(sprintf("  Found %d predicted target pairs\n", nrow(predicted_targets)))
  
  validated_targets <- NULL
  if (!skip_validation) {
    cat("Fetching validated targets...\n")
    validated_result <- tryCatch({
      multiMiR::get_multimir(org = "hsa", mirna = mirna_names_formatted,
                            table = "validated", summary = FALSE)
    }, error = function(e) NULL)
    
    if (!is.null(validated_result) && nrow(validated_result@data) > 0) {
      validated_targets <- validated_result@data
      cat(sprintf("  Found %d validated target pairs\n", nrow(validated_targets)))
    } else {
      cat("  No validated targets found\n")
    }
  }
  
  cat("\n")
  return(list(predicted = predicted_targets, validated = validated_targets))
}

# ==============================================================================
# Filter Targets - Three Modes
# ==============================================================================

filter_high_confidence_targets <- function(predicted_targets, validated_targets = NULL,
                                           skip_validation = FALSE, use_all_targets = FALSE) {
  cat("\n=== TARGET FILTERING ===\n")
  
  if (use_all_targets) {
    cat("\n*** MODE 3: ALL TARGETS (Maximum Sensitivity) ***\n")
    cat("Including ALL predicted targets without filtering.\n")
    cat("WARNING: Higher false positive rate - validate experimentally!\n\n")
    
    all_predicted <- predicted_targets %>%
      dplyr::select(mature_mirna_id, target_symbol, target_entrez) %>%
      dplyr::distinct()
    
    pred_with_counts <- predicted_targets %>%
      dplyr::group_by(mature_mirna_id, target_symbol, target_entrez) %>%
      dplyr::summarise(n_databases = dplyr::n_distinct(database),
                       databases = paste(unique(database), collapse = ";"),
                       .groups = "drop")
    
    all_predicted <- all_predicted %>%
      dplyr::left_join(pred_with_counts, 
                       by = c("mature_mirna_id", "target_symbol", "target_entrez")) %>%
      dplyr::mutate(evidence_type = "predicted_all", confidence = "exploratory")
    
    # FIX 1: Add validated targets if available
    if (!skip_validation && !is.null(validated_targets) && nrow(validated_targets) > 0) {
      validated_clean <- validated_targets %>%
        # Filter out negative evidence (no detected interaction)
        dplyr::filter(tolower(trimws(support_type)) != "negative" | is.na(support_type)) %>%
        # Filter out records with missing identifiers
        dplyr::filter(!is.na(target_symbol) & target_symbol != "" &
                        !is.na(target_entrez) & target_entrez != "") %>%
        dplyr::select(mature_mirna_id, target_symbol, target_entrez) %>%
        dplyr::distinct() %>%
        dplyr::mutate(evidence_type = "validated", confidence = "experimental",
                      n_databases = NA_integer_, databases = "validated")
      
      # FIXED: Actually merge validated targets with predicted targets
      all_predicted <- dplyr::bind_rows(all_predicted, validated_clean) %>%
        dplyr::distinct(mature_mirna_id, target_symbol, target_entrez, .keep_all = TRUE)
    }
    
    cat(sprintf("Total targets: %d\n", nrow(all_predicted)))
    cat(sprintf("Unique genes: %d\n\n", length(unique(all_predicted$target_symbol))))
    
    return(list(high_confidence_df = all_predicted, 
                target_genes = unique(all_predicted$target_symbol),
                filtering_mode = "all_targets"))
  }
  
  # MODE 1 & 2: HIGH-CONFIDENCE
  cat("\n*** MODE 1/2: HIGH-CONFIDENCE FILTERING ***\n")
  
  # Check if pubmed_id column exists
  has_pmid <- "pubmed_id" %in% colnames(predicted_targets)
  
  if (has_pmid) {
    predicted_counts <- predicted_targets %>%
      dplyr::group_by(mature_mirna_id, target_symbol, target_entrez) %>%
      dplyr::summarise(
        n_databases = dplyr::n_distinct(database),
        n_pmids = dplyr::n_distinct(pubmed_id[!is.na(pubmed_id) & pubmed_id != ""]),
        databases = paste(unique(database), collapse = ";"),
        .groups = "drop") %>%
      dplyr::filter(n_databases >= 2 | n_pmids >= 2) %>%
      dplyr::mutate(evidence_type = "predicted", confidence = "multiple_evidence")
  } else {
    # No pubmed_id column - use only database count
    predicted_counts <- predicted_targets %>%
      dplyr::group_by(mature_mirna_id, target_symbol, target_entrez) %>%
      dplyr::summarise(
        n_databases = dplyr::n_distinct(database),
        databases = paste(unique(database), collapse = ";"),
        .groups = "drop") %>%
      dplyr::filter(n_databases >= 2) %>%
      dplyr::mutate(evidence_type = "predicted", confidence = "multiple_evidence")
  }
  
  cat(sprintf("Predicted targets with multiple evidence: %d\n", nrow(predicted_counts)))
  
  # FIX 2: Corrected validation handling for MODE 1 & 2
  if (!skip_validation && !is.null(validated_targets) && nrow(validated_targets) > 0) {
    validated_counts <- validated_targets %>%
      # Only filter out "negative" support types (no need to filter functional mti!)
      dplyr::filter(tolower(trimws(support_type)) != "negative" | is.na(support_type)) %>%
      # Filter out records with missing identifiers
      dplyr::filter(!is.na(target_symbol) & target_symbol != "" &
                      !is.na(target_entrez) & target_entrez != "") %>%
      dplyr::group_by(mature_mirna_id, target_symbol, target_entrez) %>%
      dplyr::summarise(n_low_throughput = dplyr::n(), 
                       databases = "validated", 
                       .groups = "drop") %>%
      dplyr::mutate(evidence_type = "validated", 
                    confidence = "experimental",
                    n_databases = NA_integer_)
    
    cat(sprintf("Validated targets (low-throughput): %d\n", nrow(validated_counts)))
    high_conf <- dplyr::bind_rows(predicted_counts, validated_counts) %>%
      dplyr::distinct(mature_mirna_id, target_symbol, target_entrez, .keep_all = TRUE)
  } else {
    high_conf <- predicted_counts
  }
  
  cat(sprintf("\nTotal high-confidence targets: %d\n", nrow(high_conf)))
  cat(sprintf("Unique genes: %d\n\n", length(unique(high_conf$target_symbol))))
  
  return(list(high_confidence_df = high_conf, 
              target_genes = unique(high_conf$target_symbol),
              filtering_mode = if(skip_validation) "predictions_only" else "high_confidence"))
}
# ==============================================================================
# Pathway Enrichment with Storey's q-value
# ==============================================================================

run_pathway_enrichment <- function(gene_list, min_gene_count = 5, max_gene_count = 500,
                                  p_cutoff = 0.05, q_cutoff = 0.2) {
  cat("\n=== PATHWAY ENRICHMENT ===\n")
  cat("Using Storey's q-value for FDR control\n\n")
  
  gene_entrez <- clusterProfiler::bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", 
                                      OrgDb = org.Hs.eg.db::org.Hs.eg.db)
  
  if (nrow(gene_entrez) == 0) stop("No genes mapped to Entrez IDs")
  
  cat(sprintf("Mapped %d/%d genes to Entrez IDs\n\n", nrow(gene_entrez), length(gene_list)))
  
  entrez_ids <- gene_entrez$ENTREZID
  results <- list()
  qvalue_diagnostics <- list()
  
  # GO BP
  cat("Running GO Biological Process...\n")
  tryCatch({
    go_bp <- clusterProfiler::enrichGO(gene = entrez_ids, OrgDb = org.Hs.eg.db::org.Hs.eg.db,
                                      ont = "BP", pvalueCutoff = 1, qvalueCutoff = 1,
                                      minGSSize = min_gene_count, maxGSSize = max_gene_count,
                                      readable = TRUE)
    
    if (!is.null(go_bp) && nrow(go_bp@result) > 0) {
      qval_result <- calculate_storey_qvalue(go_bp@result$pvalue)
      go_bp@result$qvalue <- qval_result$qvalues
      go_bp@result$p.adjust <- qval_result$qvalues
      qvalue_diagnostics$GO_BP <- qval_result
      
      go_bp@result <- go_bp@result %>% dplyr::filter(pvalue < p_cutoff & qvalue < q_cutoff)
      results$GO_BP <- go_bp
      cat(sprintf("  Found %d significant pathways\n", nrow(go_bp@result)))
    }
  }, error = function(e) cat("  Error:", e$message, "\n"))
  
  # GO MF
  cat("Running GO Molecular Function...\n")
  tryCatch({
    go_mf <- clusterProfiler::enrichGO(gene = entrez_ids, OrgDb = org.Hs.eg.db::org.Hs.eg.db,
                                      ont = "MF", pvalueCutoff = 1, qvalueCutoff = 1,
                                      minGSSize = min_gene_count, maxGSSize = max_gene_count,
                                      readable = TRUE)
    
    if (!is.null(go_mf) && nrow(go_mf@result) > 0) {
      qval_result <- calculate_storey_qvalue(go_mf@result$pvalue)
      go_mf@result$qvalue <- qval_result$qvalues
      go_mf@result$p.adjust <- qval_result$qvalues
      qvalue_diagnostics$GO_MF <- qval_result
      
      go_mf@result <- go_mf@result %>% dplyr::filter(pvalue < p_cutoff & qvalue < q_cutoff)
      results$GO_MF <- go_mf
      cat(sprintf("  Found %d significant pathways\n", nrow(go_mf@result)))
    }
  }, error = function(e) cat("  Error:", e$message, "\n"))
  
  # GO CC
  cat("Running GO Cellular Component...\n")
  tryCatch({
    go_cc <- clusterProfiler::enrichGO(gene = entrez_ids, OrgDb = org.Hs.eg.db::org.Hs.eg.db,
                                      ont = "CC", pvalueCutoff = 1, qvalueCutoff = 1,
                                      minGSSize = min_gene_count, maxGSSize = max_gene_count,
                                      readable = TRUE)
    
    if (!is.null(go_cc) && nrow(go_cc@result) > 0) {
      qval_result <- calculate_storey_qvalue(go_cc@result$pvalue)
      go_cc@result$qvalue <- qval_result$qvalues
      go_cc@result$p.adjust <- qval_result$qvalues
      qvalue_diagnostics$GO_CC <- qval_result
      
      go_cc@result <- go_cc@result %>% dplyr::filter(pvalue < p_cutoff & qvalue < q_cutoff)
      results$GO_CC <- go_cc
      cat(sprintf("  Found %d significant pathways\n", nrow(go_cc@result)))
    }
  }, error = function(e) cat("  Error:", e$message, "\n"))
  
  # KEGG
  cat("Running KEGG pathways...\n")
  tryCatch({
    kegg <- clusterProfiler::enrichKEGG(gene = entrez_ids, organism = "hsa",
                                       pvalueCutoff = 1, qvalueCutoff = 1,
                                       minGSSize = min_gene_count, maxGSSize = max_gene_count)
    
    if (!is.null(kegg) && nrow(kegg@result) > 0) {
      cat(sprintf("  Retrieved %d KEGG pathways\n", nrow(kegg@result)))
      
      qval_result <- calculate_storey_qvalue(kegg@result$pvalue)
      kegg@result$qvalue <- qval_result$qvalues
      kegg@result$p.adjust <- qval_result$qvalues
      qvalue_diagnostics$KEGG <- qval_result
      
      n_before <- nrow(kegg@result)
      kegg@result <- kegg@result %>% dplyr::filter(pvalue < p_cutoff & qvalue < q_cutoff)
      kegg_readable <- clusterProfiler::setReadable(kegg, OrgDb = org.Hs.eg.db::org.Hs.eg.db, keyType = "ENTREZID")
      results$KEGG <- kegg_readable
      
      cat(sprintf("  Significant (p < %.2f, q < %.2f): %d\n", p_cutoff, q_cutoff, nrow(kegg_readable@result)))
      if (nrow(kegg_readable@result) == 0 && n_before > 0) {
        cat(sprintf("  NOTE: %d pathways found but none significant. Try q_cutoff = %.1f\n", n_before, q_cutoff * 2))
      }
    } else {
      cat("  No KEGG pathways retrieved\n")
      cat("  NOTE: KEGG requires internet. Check connectivity if this persists.\n")
      results$KEGG <- NULL
    }
  }, error = function(e) {
    cat("  ERROR:", e$message, "\n")
    cat("  KEGG enrichment failed. Check internet connection.\n")
    results$KEGG <- NULL
  })
  
  # Reactome
  cat("Running Reactome pathways...\n")
  tryCatch({
    reactome <- ReactomePA::enrichPathway(gene = entrez_ids, organism = "human",
                                         pvalueCutoff = 1, qvalueCutoff = 1,
                                         minGSSize = min_gene_count, maxGSSize = max_gene_count,
                                         readable = TRUE)
    
    if (!is.null(reactome) && nrow(reactome@result) > 0) {
      qval_result <- calculate_storey_qvalue(reactome@result$pvalue)
      reactome@result$qvalue <- qval_result$qvalues
      reactome@result$p.adjust <- qval_result$qvalues
      qvalue_diagnostics$Reactome <- qval_result
      
      reactome@result <- reactome@result %>% dplyr::filter(pvalue < p_cutoff & qvalue < q_cutoff)
      results$Reactome <- reactome
      cat(sprintf("  Found %d significant pathways\n", nrow(reactome@result)))
    }
  }, error = function(e) cat("  Error:", e$message, "\n"))
  
  cat("\nEnrichment complete.\n")
  return(list(results = results, gene_mapping = gene_entrez, qvalue_diagnostics = qvalue_diagnostics))
}

# ==============================================================================
# Disease-Specific Pathway Filtering
# ==============================================================================

extract_disease_pathways <- function(enrichment_results, disease_keywords, q_cutoff = 0.05) {
  cat("\n=== EXTRACTING DISEASE-RELEVANT PATHWAYS ===\n")
  cat("Keywords:", paste(disease_keywords, collapse = ", "), "\n\n")
  
  results_list <- enrichment_results$results
  filtered_by_db <- list()
  disease_qval_diagnostics <- list()
  
  for (db_name in names(results_list)) {
    db_result <- results_list[[db_name]]
    if (is.null(db_result) || nrow(db_result@result) == 0) next
    
    cat(sprintf("\n%s: %d total pathways\n", db_name, nrow(db_result@result)))
    
    pattern <- paste(disease_keywords, collapse = "|")
    disease_pathways <- db_result@result %>%
      dplyr::filter(grepl(pattern, Description, ignore.case = TRUE))
    
    if (nrow(disease_pathways) == 0) {
      cat("  No disease-relevant pathways\n")
      next
    }
    
    cat(sprintf("  Disease-relevant: %d\n", nrow(disease_pathways)))
    
    # Recalculate q-values
    qval_result <- calculate_storey_qvalue(disease_pathways$pvalue)
    disease_pathways$qvalue_disease <- qval_result$qvalues
    disease_pathways$qvalue_original <- disease_pathways$qvalue
    disease_pathways$significant_disease <- disease_pathways$qvalue_disease < q_cutoff
    disease_pathways$significant_original <- disease_pathways$qvalue_original < q_cutoff
    disease_pathways$Database <- db_name
    
    disease_qval_diagnostics[[db_name]] <- qval_result
    
    newly_sig <- sum(disease_pathways$significant_disease & !disease_pathways$significant_original)
    if (newly_sig > 0) cat(sprintf("  Newly significant: %d\n", newly_sig))
    
    cat(sprintf("  Significant (disease FDR < %.2f): %d\n", q_cutoff, 
               sum(disease_pathways$significant_disease)))
    
    filtered_by_db[[db_name]] <- disease_pathways
  }
  
  if (length(filtered_by_db) > 0) {
    combined <- dplyr::bind_rows(filtered_by_db) %>% dplyr::arrange(qvalue_disease)
    cat("\n=== SUMMARY ===\n")
    cat(sprintf("Total disease-relevant: %d\n", nrow(combined)))
    cat(sprintf("Significant (FDR < %.2f): %d\n\n", q_cutoff, sum(combined$significant_disease)))
  } else {
    combined <- NULL
    cat("\nNo disease-relevant pathways found.\n\n")
  }
  
  return(list(disease_pathways = combined, pathways_by_database = filtered_by_db,
             disease_qvalue_diagnostics = disease_qval_diagnostics,
             keywords_used = disease_keywords))
}

# ==============================================================================
# Save Results
# ==============================================================================

save_results <- function(results, output_dir, output_prefix) {
  cat("\n=== SAVING RESULTS ===\n")
  files_created <- c()
  
  # High-confidence targets
  if (!is.null(results$high_confidence_targets$high_confidence_df)) {
    f <- file.path(output_dir, paste0(output_prefix, "_high_confidence_targets.csv"))
    write.csv(results$high_confidence_targets$high_confidence_df, f, row.names = FALSE)
    files_created <- c(files_created, f)
    cat("✓", basename(f), "\n")
  }
  
  # Target genes
  if (!is.null(results$high_confidence_targets$target_genes)) {
    f <- file.path(output_dir, paste0(output_prefix, "_target_genes.csv"))
    write.csv(data.frame(gene = results$high_confidence_targets$target_genes), f, row.names = FALSE)
    files_created <- c(files_created, f)
    cat("✓", basename(f), "\n")
  }
  
  # Enrichment results
  for (db_name in names(results$enrichment$results)) {
    db_result <- results$enrichment$results[[db_name]]
    if (!is.null(db_result) && nrow(db_result@result) > 0) {
      f <- file.path(output_dir, paste0(output_prefix, "_", db_name, "_enrichment.csv"))
      write.csv(db_result@result, f, row.names = FALSE)
      files_created <- c(files_created, f)
      cat("✓", basename(f), "\n")
    }
  }
  
  # Disease pathways
  if (!is.null(results$disease_pathways$disease_pathways)) {
    f <- file.path(output_dir, paste0(output_prefix, "_disease_relevant_pathways.csv"))
    write.csv(results$disease_pathways$disease_pathways, f, row.names = FALSE)
    files_created <- c(files_created, f)
    cat("✓", basename(f), "\n")
  }
  
  # Q-value diagnostics
  if (!is.null(results$qvalue_diagnostics)) {
    qval_df <- do.call(rbind, lapply(names(results$qvalue_diagnostics), function(db) {
      diag <- results$qvalue_diagnostics[[db]]
      data.frame(Database = db, Method = diag$method_used,
                Pi0 = ifelse(is.null(diag$pi0) || is.na(diag$pi0), NA, diag$pi0),
                N_tests = diag$n_valid)
    }))
    f <- file.path(output_dir, paste0(output_prefix, "_qvalue_diagnostics.csv"))
    write.csv(qval_df, f, row.names = FALSE)
    files_created <- c(files_created, f)
    cat("✓", basename(f), "\n")
  }
  
  cat("\n")
  return(files_created)
}

# ==============================================================================
# Create Visualizations
# ==============================================================================

create_pathway_plots <- function(results, output_dir, output_prefix) {
  cat("\n=== CREATING VISUALIZATIONS ===\n")
  
  plots_created <- c()
  
  # Function to safely create and save dotplot
  create_dotplot <- function(enrich_obj, db_name, title) {
    if (is.null(enrich_obj) || nrow(enrich_obj@result) == 0) {
      return(NULL)
    }
    
    tryCatch({
      # Create dotplot showing top pathways
      n_show <- min(20, nrow(enrich_obj@result))
      
      p <- enrichplot::dotplot(enrich_obj, showCategory = n_show, 
                              title = title) +
        ggplot2::theme(axis.text.y = ggplot2::element_text(size = 10))
      
      # Save plot
      # MANUSCRIPT: Fig_5B - GO BP dotplot (search: GO_BP_dotplot)
      filename <- file.path(output_dir, paste0(output_prefix, "_", db_name, "_dotplot.png"))
      ggplot2::ggsave(filename, p, width = 10, height = 8, dpi = 300)
      cat("✓", basename(filename), "\n")
      # Also save PDF
      filename_pdf <- file.path(output_dir, paste0(output_prefix, "_", db_name, "_dotplot.pdf"))
      ggplot2::ggsave(filename_pdf, p, width = 10, height = 8)
      return(filename)
    }, error = function(e) {
      cat("  Error creating", db_name, "plot:", e$message, "\n")
      return(NULL)
    })
  }
  
  # Create dotplots for each database
  for (db_name in names(results$enrichment$results)) {
    db_result <- results$enrichment$results[[db_name]]
    if (!is.null(db_result) && nrow(db_result@result) > 0) {
      title <- switch(db_name,
                     "GO_BP" = "GO Biological Process",
                     "GO_MF" = "GO Molecular Function", 
                     "GO_CC" = "GO Cellular Component",
                     "KEGG" = "KEGG Pathways",
                     "Reactome" = "Reactome Pathways",
                     db_name)
      
      plot_file <- create_dotplot(db_result, db_name, title)
      if (!is.null(plot_file)) {
        plots_created <- c(plots_created, plot_file)
      }
    }
  }
  
  # Create combined pathway plot if multiple databases have results
  tryCatch({
    databases_with_results <- names(results$enrichment$results)[
      sapply(results$enrichment$results, function(x) !is.null(x) && nrow(x@result) > 0)
    ]
    
    if (length(databases_with_results) >= 2) {
      # Combine top pathways from each database
      combined_df <- data.frame()
      
      for (db_name in databases_with_results) {
        db_result <- results$enrichment$results[[db_name]]@result
        db_top <- db_result %>%
          dplyr::arrange(qvalue) %>%
          head(10) %>%
          dplyr::mutate(Database = db_name,
                       GeneRatio_numeric = sapply(strsplit(GeneRatio, "/"), 
                                                 function(x) as.numeric(x[1])/as.numeric(x[2])))
        combined_df <- rbind(combined_df, db_top)
      }
      
      if (nrow(combined_df) > 0) {
        # Create combined plot
        p_combined <- ggplot2::ggplot(combined_df, 
                                      ggplot2::aes(x = GeneRatio_numeric, 
                                                  y = reorder(Description, GeneRatio_numeric),
                                                  color = qvalue, size = Count)) +
          ggplot2::geom_point() +
          ggplot2::scale_color_gradient(low = "red", high = "blue") +
          ggplot2::facet_wrap(~Database, scales = "free_y", ncol = 1) +
          ggplot2::labs(title = "Top Enriched Pathways Across Databases",
                       x = "Gene Ratio", y = "Pathway", color = "Q-value", size = "Gene Count") +
          ggplot2::theme_bw() +
          ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8),
                        strip.text = ggplot2::element_text(face = "bold"))
        
        combined_file <- file.path(output_dir, paste0(output_prefix, "_combined_pathways.png"))
        ggplot2::ggsave(combined_file, p_combined, width = 12, height = 14, dpi = 300)
        plots_created <- c(plots_created, combined_file)
        cat("✓", basename(combined_file), "\n")
      }
    }
  }, error = function(e) {
    cat("  Error creating combined plot:", e$message, "\n")
  })
  
  # Create disease-specific plot if available
  if (!is.null(results$disease_pathways$disease_pathways)) {
    tryCatch({
      disease_df <- results$disease_pathways$disease_pathways %>%
        dplyr::filter(significant_disease) %>%
        dplyr::arrange(qvalue_disease) %>%
        head(30)
      
      if (nrow(disease_df) > 0) {
        disease_df <- disease_df %>%
          dplyr::mutate(GeneRatio_numeric = sapply(strsplit(GeneRatio, "/"), 
                                                   function(x) as.numeric(x[1])/as.numeric(x[2])))
        
        p_disease <- ggplot2::ggplot(disease_df,
                                     ggplot2::aes(x = GeneRatio_numeric,
                                                 y = reorder(Description, GeneRatio_numeric),
                                                 color = qvalue_disease, size = Count)) +
          ggplot2::geom_point() +
          ggplot2::scale_color_gradient(low = "red", high = "blue") +
          ggplot2::facet_wrap(~Database, scales = "free_y", ncol = 1) +
          ggplot2::labs(title = "Top Disease-Relevant Pathways",
                       x = "Gene Ratio", y = "Pathway", 
                       color = "Q-value (disease)", size = "Gene Count") +
          ggplot2::theme_bw() +
          ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8),
                        strip.text = ggplot2::element_text(face = "bold"))
        
        disease_file <- file.path(output_dir, paste0(output_prefix, "_disease_pathways.png"))
        ggplot2::ggsave(disease_file, p_disease, width = 12, height = 16, dpi = 300)
        plots_created <- c(plots_created, disease_file)
        cat("✓", basename(disease_file), "\n")
      }
    }, error = function(e) {
      cat("  Error creating disease plot:", e$message, "\n")
    })
  }
  
  cat("\n")
  return(plots_created)
}

# ==============================================================================
# MAIN FUNCTION
# ==============================================================================

analyze_mirna_pathways <- function(mirna_names,
                                  output_prefix = "miRNA_analysis",
                                  output_dir = ".",
                                  min_gene_count = 5,
                                  max_gene_count = 500,
                                  p_cutoff = 0.05,
                                  q_cutoff = 0.2,
                                  disease_keywords = NULL,
                                  cf_keywords = NULL,  # backward compatibility
                                  skip_validation = FALSE,
                                  use_all_targets = FALSE) {
  
  cat("\n================================================================================\n")
  cat("                    miRNA Pathway Analysis Workflow                            \n")
  cat("================================================================================\n\n")
  
  # Backward compatibility
  if (!is.null(cf_keywords) && is.null(disease_keywords)) {
    warning("'cf_keywords' is deprecated. Use 'disease_keywords' instead.")
    disease_keywords <- cf_keywords
  }
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat("Created output directory:", output_dir, "\n\n")
  }
  
  # Parameters
  params <- list(mirna_names = mirna_names, output_prefix = output_prefix,
                output_dir = output_dir, disease_keywords = disease_keywords,
                skip_validation = skip_validation, use_all_targets = use_all_targets,
                timestamp = Sys.time())
  
  cat("Analysis Parameters:\n")
  cat("  miRNAs:", paste(mirna_names, collapse = ", "), "\n")
  cat("  Filtering mode:", if(use_all_targets) "ALL TARGETS (Mode 3)" 
      else if(skip_validation) "PREDICTIONS ONLY (Mode 2)" 
      else "HIGH-CONFIDENCE (Mode 1)", "\n")
  cat("  FDR method: Storey's q-value\n")
  if (!is.null(disease_keywords)) {
    cat("  Disease keywords:", paste(disease_keywords, collapse = ", "), "\n")
  }
  cat("\n")
  
  # Fetch targets
  targets <- fetch_mirna_targets(mirna_names, skip_validation)
  
  # Filter targets
  filtered_targets <- filter_high_confidence_targets(
    predicted_targets = targets$predicted,
    validated_targets = targets$validated,
    skip_validation = skip_validation,
    use_all_targets = use_all_targets
  )
  
  # Pathway enrichment
  enrichment_results <- run_pathway_enrichment(
    gene_list = filtered_targets$target_genes,
    min_gene_count = min_gene_count,
    max_gene_count = max_gene_count,
    p_cutoff = p_cutoff,
    q_cutoff = q_cutoff
  )
  
  # Disease-specific filtering
  disease_pathways <- NULL
  if (!is.null(disease_keywords)) {
    disease_pathways <- extract_disease_pathways(enrichment_results, disease_keywords, q_cutoff)
  }
  
  # Compile results
  results <- list(
    parameters = params,
    targets_data = targets,
    high_confidence_targets = filtered_targets,
    enrichment = enrichment_results,
    disease_pathways = disease_pathways,
    qvalue_diagnostics = enrichment_results$qvalue_diagnostics
  )
  
  # Save files
  files_created <- save_results(results, output_dir, output_prefix)
  results$files_created <- files_created
  
  # Create visualizations
  plots_created <- create_pathway_plots(results, output_dir, output_prefix)
  results$plots_created <- plots_created
  
  # Backward compatibility
  results$cf_pathways <- disease_pathways
  
  cat("\n================================================================================\n")
  cat("                           ANALYSIS COMPLETE                                    \n")
  cat("================================================================================\n\n")
  cat("Output directory:", output_dir, "\n")
  cat("CSV files created:", length(files_created), "\n")
  cat("Plots created:", length(plots_created), "\n\n")
  
  return(results)
}

cat("Workflow ready. Type ?analyze_mirna_pathways for help.\n\n")
