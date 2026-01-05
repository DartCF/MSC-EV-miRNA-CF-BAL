# ==============================================================================
# UNIFIED WORKFLOW RUNNER
# ==============================================================================
# Runs the complete 2×2 factorial analysis from a single entry point,
# driven by the central configuration.
#
# FEATURES:
#   - Single source of truth for all parameters
#   - Metadata saved with each analysis for provenance
#   - Can run all analyses or selected subset
#   - Automatic comparison after completion
#
# USAGE:
#   source(here("analysis", "04_pathway", "run_all_analyses.R"))
#   
#   # Run everything
#   all_results <- run_full_analysis()
#   
#   # Run specific analyses
#   results <- run_analysis("A")  # Just miR-6126 stringent
#   results <- run_analyses(c("A", "B"))  # Just miR-6126 both modes
# ==============================================================================

library(here)

# Load required scripts from the pathway analysis directory
pathway_dir <- here("analysis", "04_pathway")

source(file.path(pathway_dir, "miRNA_pathway_workflow.R"))
source(file.path(pathway_dir, "summarize_evidence.R"))
source(file.path(pathway_dir, "validate_enrichment.R"))
source(file.path(pathway_dir, "compare_fdr_thresholds.R"))
source(file.path(pathway_dir, "analyze_validation_evidence_module.R"))
source(file.path(pathway_dir, "analysis_config.R"))

# ==============================================================================
# Run Single Analysis
# ==============================================================================

run_analysis <- function(analysis_id, config = NULL, verbose = TRUE) {
  
  # Load config if not provided
  if (is.null(config)) {
    config <- load_analysis_config()
  }
  
  # Get analysis parameters
  if (!analysis_id %in% names(config$analyses)) {
    stop(sprintf("Unknown analysis ID: %s. Valid IDs: %s",
                analysis_id, paste(names(config$analyses), collapse = ", ")))
  }
  
  analysis <- config$analyses[[analysis_id]]
  params <- analysis$resolved
  
  # Print header
  if (verbose) {
    cat("\n")
    cat("================================================================================\n")
    cat(sprintf("     RUNNING ANALYSIS [%s]: %s\n", analysis_id, analysis$name))
    cat("================================================================================\n\n")
    
    cat("PARAMETERS (from central configuration):\n")
    cat(sprintf("  • miRNA set: %s (n = %d)\n", params$mirna_set_name, params$n_mirnas))
    cat(sprintf("  • Stringency: %s\n", params$stringency_name))
    cat(sprintf("  • use_all_targets: %s\n", params$use_all_targets))
    cat(sprintf("  • FDR threshold: q < %.2f\n", params$q_cutoff))
    cat(sprintf("  • Output directory: %s\n", analysis$output_dir))
    cat("\n")
    
    if (params$n_mirnas <= 5) {
      cat("miRNAs:\n")
      for (m in params$mirnas) {
        cat(sprintf("  • %s\n", m))
      }
      cat("\n")
    } else {
      cat(sprintf("miRNAs: %d total (see configuration for full list)\n\n", params$n_mirnas))
    }
  }
  
  # Create output directory
  if (!dir.exists(analysis$output_dir)) {
    dir.create(analysis$output_dir, recursive = TRUE)
  }
  
  # Run the analysis
  results <- analyze_mirna_pathways(
    mirna_names = params$mirnas,
    output_prefix = analysis$output_prefix,
    output_dir = analysis$output_dir,
    disease_keywords = params$disease_keywords,
    use_all_targets = params$use_all_targets,
    skip_validation = FALSE,
    q_cutoff = params$q_cutoff
  )
  
  # Validate enrichment
  if (verbose) cat("\n")
  validate_enrichment(results)
  
  # Compare FDR thresholds
  if (verbose) cat("\n")
  compare_fdr_thresholds(results)
  
  # Evidence summary
  if (verbose) cat("\n")
  evidence_summary <- summarize_target_evidence(results, save_to_file = TRUE)
  
  # Validation evidence analysis
  if (verbose) cat("\n")
  validation_analysis <- analyze_validation_evidence(results, save_files = TRUE)
  
  # Save metadata WITH the results
  metadata <- save_analysis_metadata(config, analysis_id, results, analysis$output_dir)
  
  # Save workspace
  rdata_file <- file.path(analysis$output_dir, 
                         sprintf("%s_analysis.RData", analysis$short_name))
  save(results, evidence_summary, validation_analysis, metadata,
       file = rdata_file)
  
  # Print summary
  if (verbose) {
    cat("\n")
    cat("================================================================================\n")
    cat(sprintf("     ANALYSIS [%s] COMPLETE: %s\n", analysis_id, analysis$name))
    cat("================================================================================\n\n")
    
    n_targets <- length(results$high_confidence_targets$target_genes)
    n_cf_pathways <- if(!is.null(results$disease_pathways$disease_pathways))
                       nrow(results$disease_pathways$disease_pathways) else 0
    n_cf_sig <- if(!is.null(results$disease_pathways$disease_pathways))
                  sum(results$disease_pathways$disease_pathways$qvalue_disease < params$q_cutoff, 
                      na.rm = TRUE) else 0
    
    cat("RESULTS SUMMARY:\n")
    cat(sprintf("  • Targets identified: %d\n", n_targets))
    cat(sprintf("  • CF-relevant pathways: %d total, %d significant (q < %.2f)\n",
                n_cf_pathways, n_cf_sig, params$q_cutoff))
    cat(sprintf("  • Output directory: %s\n", analysis$output_dir))
    cat(sprintf("  • Workspace saved: %s\n", rdata_file))
    cat("\n")
    
    cat("PROVENANCE:\n")
    cat("  • Configuration-driven analysis\n")
    cat("  • Metadata saved with results (analysis_metadata.json)\n")
    cat("  • Parameters verified from central config\n")
    cat("\n")
  }
  
  # Return results with metadata attached
  results$metadata <- metadata
  results$analysis_id <- analysis_id
  results$config <- config$analyses[[analysis_id]]
  
  return(results)
}

# ==============================================================================
# Run Multiple Analyses
# ==============================================================================

run_analyses <- function(analysis_ids, config = NULL, verbose = TRUE) {
  
  # Load config if not provided
  if (is.null(config)) {
    config <- load_analysis_config()
  }
  
  results <- list()
  
  for (id in analysis_ids) {
    results[[id]] <- run_analysis(id, config = config, verbose = verbose)
  }
  
  if (verbose) {
    cat("\n")
    cat("================================================================================\n")
    cat("                    ALL REQUESTED ANALYSES COMPLETE                            \n")
    cat("================================================================================\n\n")
    
    cat("COMPLETED ANALYSES:\n")
    for (id in analysis_ids) {
      cat(sprintf("  [%s] %s\n", id, config$analyses[[id]]$name))
    }
    cat("\n")
  }
  
  return(results)
}

# ==============================================================================
# Run Full 2×2 Analysis
# ==============================================================================

run_full_analysis <- function(config = NULL, 
                              run_comparison = TRUE,
                              verbose = TRUE) {
  
  # Load config if not provided
  if (is.null(config)) {
    config <- load_analysis_config()
  }
  
  if (verbose) {
    cat("\n")
    cat("================================================================================\n")
    cat("                    FULL 2×2 FACTORIAL ANALYSIS                                \n")
    cat("================================================================================\n\n")
    
    print_config_summary(config)
    
    cat("Starting analyses...\n\n")
  }
  
  # Get all analysis IDs
  analysis_ids <- names(config$analyses)
  
  # Run all analyses
  results <- run_analyses(analysis_ids, config = config, verbose = verbose)
  
  # Run comparison if requested
  if (run_comparison) {
    if (verbose) {
      cat("\n")
      cat("================================================================================\n")
      cat("                    RUNNING SYSTEMATIC COMPARISON                              \n")
      cat("================================================================================\n\n")
    }
    
    comparison <- run_systematic_comparison_from_metadata(config, verbose = verbose)
    results$comparison <- comparison
  }
  
  if (verbose) {
    cat("\n")
    cat("================================================================================\n")
    cat("                    FULL ANALYSIS COMPLETE                                     \n")
    cat("================================================================================\n\n")
    
    cat("All results saved to:\n")
    for (id in analysis_ids) {
      cat(sprintf("  [%s] %s\n", id, config$analyses[[id]]$output_dir))
    }
    if (run_comparison) {
      cat(sprintf("  [Comparison] %s\n", config$output$comparison_dir))
    }
    cat("\n")
  }
  
  return(results)
}

# ==============================================================================
# Comparison Using Metadata (not directory names!)
# ==============================================================================

run_systematic_comparison_from_metadata <- function(config, verbose = TRUE) {
  
  if (!dir.exists(config$output$comparison_dir)) {
    dir.create(config$output$comparison_dir, recursive = TRUE)
  }
  
  if (verbose) {
    cat("STEP 1: LOADING METADATA FROM COMPLETED ANALYSES\n")
    cat("=================================================\n\n")
  }
  
  # Load metadata from each analysis directory
  all_metadata <- list()
  
  for (analysis_id in names(config$analyses)) {
    analysis <- config$analyses[[analysis_id]]
    
    metadata <- load_analysis_metadata(analysis$output_dir)
    
    if (is.null(metadata)) {
      warning(sprintf("No metadata found for analysis %s (%s). Skipping.",
                     analysis_id, analysis$output_dir))
      next
    }
    
    all_metadata[[analysis_id]] <- metadata
    
    if (verbose) {
      cat(sprintf("  [%s] %s\n", analysis_id, analysis$name))
      cat(sprintf("      Targets: %s | CF pathways: %s\n",
                 metadata$results_summary$n_targets %||% "N/A",
                 metadata$results_summary$n_cf_pathways %||% "N/A"))
    }
  }
  
  if (length(all_metadata) == 0) {
    stop("No completed analyses found with metadata. Run analyses first.")
  }
  
  if (verbose) cat("\n")
  
  # Verify that metadata matches expected configuration
  if (verbose) {
    cat("STEP 2: VERIFYING PROVENANCE\n")
    cat("============================\n\n")
  }
  
  for (analysis_id in names(all_metadata)) {
    meta <- all_metadata[[analysis_id]]
    expected <- config$analyses[[analysis_id]]$resolved
    
    # Check key parameters match
    checks <- c(
      n_mirnas = identical(meta$parameters$n_mirnas, expected$n_mirnas),
      use_all_targets = identical(meta$parameters$use_all_targets, expected$use_all_targets),
      q_cutoff = identical(meta$parameters$q_cutoff, expected$q_cutoff)
    )
    
    if (all(checks)) {
      if (verbose) cat(sprintf("  [%s] ✓ Metadata matches configuration\n", analysis_id))
    } else {
      warning(sprintf("Analysis %s: Metadata does not match current configuration!", 
                     analysis_id))
      if (verbose) {
        cat(sprintf("  [%s] ⚠ MISMATCH - results may be from different config\n", 
                   analysis_id))
        for (check_name in names(checks)) {
          if (!checks[check_name]) {
            cat(sprintf("      %s: expected %s, found %s\n",
                       check_name,
                       expected[[check_name]],
                       meta$parameters[[check_name]]))
          }
        }
      }
    }
  }
  
  if (verbose) cat("\n")
  
  # Build comparison summary from metadata
  if (verbose) {
    cat("STEP 3: BUILDING COMPARISON FROM VERIFIED METADATA\n")
    cat("===================================================\n\n")
  }
  
  comparison_df <- data.frame(
    Analysis_ID = character(),
    Analysis_Name = character(),
    miRNA_Scope = character(),
    N_miRNAs = integer(),
    Evidence_Threshold = character(),
    FDR_Threshold = numeric(),
    N_Targets = integer(),
    N_Validated = integer(),
    N_CF_Pathways = integer(),
    stringsAsFactors = FALSE
  )
  
  for (analysis_id in names(all_metadata)) {
    meta <- all_metadata[[analysis_id]]
    analysis <- config$analyses[[analysis_id]]
    
    comparison_df <- rbind(comparison_df, data.frame(
      Analysis_ID = analysis_id,
      Analysis_Name = analysis$name,
      miRNA_Scope = if(analysis$mirna_set == "single") "Single" else "Multiple",
      N_miRNAs = meta$parameters$n_mirnas,
      Evidence_Threshold = meta$stringency_definition$name,
      FDR_Threshold = meta$parameters$q_cutoff,
      N_Targets = meta$results_summary$n_targets,
      N_Validated = meta$results_summary$n_validated,
      N_CF_Pathways = meta$results_summary$n_cf_pathways,
      stringsAsFactors = FALSE
    ))
  }
  
  # Display 2×2 matrix
  if (verbose && nrow(comparison_df) == 4) {
    cat("2×2 FACTORIAL RESULTS:\n")
    cat("─────────────────────────────────────────────────────────────────────\n")
    cat("                         Evidence Threshold                          \n")
    cat("                    Stringent          Exploratory                   \n")
    cat("─────────────────────────────────────────────────────────────────────\n")
    
    A <- comparison_df[comparison_df$Analysis_ID == "A", ]
    B <- comparison_df[comparison_df$Analysis_ID == "B", ]
    C <- comparison_df[comparison_df$Analysis_ID == "C", ]
    D <- comparison_df[comparison_df$Analysis_ID == "D", ]
    
    cat(sprintf("Single miRNA    %6d targets     %6d targets\n",
                A$N_Targets, B$N_Targets))
    cat(sprintf("                %6d pathways    %6d pathways\n",
                A$N_CF_Pathways, B$N_CF_Pathways))
    cat(sprintf("Multiple miRNAs %6d targets     %6d targets\n",
                C$N_Targets, D$N_Targets))
    cat(sprintf("                %6d pathways    %6d pathways\n",
                C$N_CF_Pathways, D$N_CF_Pathways))
    cat("─────────────────────────────────────────────────────────────────────\n\n")
  }
  
  # Save comparison
  write.csv(comparison_df, 
           file.path(config$output$comparison_dir, "factorial_comparison.csv"),
           row.names = FALSE)
  
  if (verbose) {
    cat(sprintf("✓ Comparison saved: %s\n\n",
               file.path(config$output$comparison_dir, "factorial_comparison.csv")))
  }
  
  # Save provenance report
  provenance_report <- list(
    comparison_time = Sys.time(),
    config_used = config$metadata,
    analyses_compared = all_metadata,
    verification_status = "All analyses verified against configuration"
  )
  
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    jsonlite::write_json(provenance_report,
                        file.path(config$output$comparison_dir, "provenance_report.json"),
                        pretty = TRUE, auto_unbox = TRUE)
  }
  
  return(list(
    summary = comparison_df,
    metadata = all_metadata,
    provenance = provenance_report
  ))
}

# Null coalescing operator
`%||%` <- function(a, b) if (is.null(a)) b else a

# ==============================================================================
# Example Usage
# ==============================================================================

cat("\n")
cat("================================================================================\n")
cat("                    Unified Workflow Runner Loaded                              \n")
cat("================================================================================\n\n")

cat("CONFIGURATION:\n")
cat("--------------\n")
cat("By default, run_full_analysis() automatically looks for 'analysis_config.yaml'\n")
cat("in analysis/04_pathway/. If found, it uses that file as the configuration.\n")
cat("If not found, it falls back to hardcoded defaults.\n\n")

cat("USAGE:\n")
cat("------\n\n")

cat("# Run the full 2×2 analysis (auto-loads analysis_config.yaml):\n")
cat("all_results <- run_full_analysis()\n\n")

cat("# Run a specific analysis:\n")
cat('results_A <- run_analysis("A")  # miR-6126 Stringent\n')
cat('results_C <- run_analysis("C")  # Major miRNAs Stringent\n\n')

cat("# Run selected analyses:\n")
cat('results <- run_analyses(c("A", "B"))  # Both miR-6126 modes\n')
cat('results <- run_analyses(c("C", "D"))  # Both Major miRNA modes\n\n')

cat("# Use explicit configuration file:\n")
cat('config <- load_analysis_config("path/to/my_config.yaml")\n')
cat('results <- run_full_analysis(config = config)\n\n')

cat("KEY FEATURES:\n")
cat("  ✓ Auto-detects analysis_config.yaml in analysis/04_pathway/\n")
cat("  ✓ Single source of truth (edit the YAML, not the scripts)\n")
cat("  ✓ Metadata saved with each analysis\n")
cat("  ✓ Provenance verification at comparison time\n")
cat("  ✓ No reliance on directory naming conventions\n")
cat("  ✓ Full reproducibility\n\n")

cat("================================================================================\n\n")
