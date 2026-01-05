# ==============================================================================
# Validation Evidence Analysis Functions (CORRECTED VERSION)
# ==============================================================================
# Purpose: Analyze and summarize experimental validation evidence
# Integrates with miRNA_pathway_workflow.R output structure
#
# FIXES APPLIED:
# 1. Excludes "negative" support types (experiments that found NO interaction)
# 2. Filters out targets with missing identifiers
# 3. REMOVED unnecessary "functional mti" filter (functional MTI is positive evidence!)
# ==============================================================================

library(multiMiR)
library(dplyr)

analyze_validation_evidence <- function(results, save_files = TRUE) {
  
  cat("\n================================================================================\n")
  cat("           VALIDATION EVIDENCE ANALYSIS                                        \n")
  cat("================================================================================\n\n")
  
  # Extract parameters
  mirna_names <- results$parameters$mirna_names
  output_dir <- results$parameters$output_dir
  output_prefix <- results$parameters$output_prefix
  
  # Check if we have validated targets
  if (is.null(results$targets_data$validated) || 
      nrow(results$targets_data$validated) == 0) {
    cat("⚠ No validated targets found in results.\n")
    cat("  This analysis was likely run with skip_validation = TRUE\n")
    cat("  or no experimental validations exist for this miRNA.\n\n")
    return(NULL)
  }
  
  validated_data <- results$targets_data$validated
  cat(sprintf("Analyzing validation evidence for: %s\n", 
              paste(mirna_names, collapse = ", ")))
  cat(sprintf("Total validation records: %d\n\n", nrow(validated_data)))
  
  # ==============================================================================
  # 1. DATABASE SOURCES
  # ==============================================================================
  
  cat("DATABASE SOURCES OF VALIDATION\n")
  cat("================================\n")
  
  db_summary <- validated_data %>%
    group_by(database) %>%
    summarise(
      n_records = n(),
      n_unique_targets = n_distinct(target_symbol),
      .groups = "drop"
    ) %>%
    arrange(desc(n_records))
  
  print(db_summary)
  cat("\n")
  
  # ==============================================================================
  # 2. EXPERIMENTAL EVIDENCE TYPES
  # ==============================================================================
  
  cat("EXPERIMENTAL EVIDENCE TYPES\n")
  cat("============================\n\n")
  
  available_cols <- names(validated_data)
  
  exp_summary <- NULL
  support_summary <- NULL
  
  if ("experiment" %in% available_cols) {
    exp_summary <- validated_data %>%
      filter(!is.na(experiment) & experiment != "") %>%
      group_by(experiment) %>%
      summarise(
        n_records = n(),
        n_targets = n_distinct(target_symbol),
        databases = paste(unique(database), collapse = ", "),
        .groups = "drop"
      ) %>%
      arrange(desc(n_records))
    
    cat("EXPERIMENTAL METHODS (top 15):\n")
    print(head(exp_summary, 15))
    cat("\n")
  }
  
  if ("support_type" %in% available_cols) {
    support_summary <- validated_data %>%
      filter(!is.na(support_type) & support_type != "") %>%
      group_by(support_type) %>%
      summarise(
        n_records = n(),
        n_targets = n_distinct(target_symbol),
        databases = paste(unique(database), collapse = ", "),
        .groups = "drop"
      ) %>%
      arrange(desc(n_records))
    
    cat("SUPPORT TYPES:\n")
    print(support_summary)
    cat("\n")
  }
  
  # ==============================================================================
  # 3. PER-TARGET EVIDENCE SUMMARY (CORRECTED)
  # ==============================================================================
  
  cat("PER-TARGET VALIDATION EVIDENCE\n")
  cat("================================\n\n")
  
  # CORRECTED: Only exclude "negative" support types
  # Functional MTI, Non-Functional MTI, and positive are all KEPT as valid evidence
  target_evidence <- validated_data %>%
    filter(tolower(trimws(support_type)) != "negative" | is.na(support_type)) %>%
    filter(!is.na(target_symbol) & target_symbol != "" & 
           !is.na(target_entrez) & target_entrez != "") %>%
    group_by(target_symbol, target_entrez) %>%
    summarise(
      n_validation_records = n(),
      databases = paste(unique(database), collapse = "; "),
      n_databases = n_distinct(database),
      experiments = if("experiment" %in% names(validated_data)) {
        paste(unique(experiment[!is.na(experiment) & experiment != ""]), collapse = "; ")
      } else {
        NA_character_
      },
      # Only exclude negative support types from summary
      support_types = if("support_type" %in% names(validated_data)) {
        valid_support <- support_type[!is.na(support_type) & 
                                      support_type != "" & 
                                      tolower(trimws(support_type)) != "negative"]
        if(length(valid_support) > 0) {
          paste(unique(valid_support), collapse = "; ")
        } else {
          NA_character_
        }
      } else {
        NA_character_
      },
      .groups = "drop"
    ) %>%
    arrange(desc(n_validation_records))
  
  cat("TOP 20 TARGETS BY NUMBER OF VALIDATION RECORDS:\n")
  cat("(Excluding negative evidence only)\n\n")
  print(head(target_evidence, 20), n = 20)
  cat("\n")
  
  # Report filtering statistics
  n_total_records <- nrow(validated_data)
  n_negative <- sum(tolower(trimws(validated_data$support_type)) == "negative", na.rm = TRUE)
  n_missing_ids <- sum(is.na(validated_data$target_symbol) | validated_data$target_symbol == "" |
                       is.na(validated_data$target_entrez) | validated_data$target_entrez == "")
  
  cat("FILTERING STATISTICS:\n")
  cat(sprintf("  Total validation records: %d\n", n_total_records))
  cat(sprintf("  Records with negative evidence: %d (excluded)\n", n_negative))
  cat(sprintf("  Records with missing identifiers: %d (excluded)\n", n_missing_ids))
  cat(sprintf("  Records included in summary: %d\n\n", nrow(target_evidence)))
  
  # ==============================================================================
  # 4. CLASSIFICATION OF EVIDENCE TYPES
  # ==============================================================================
  
  cat("CLASSIFICATION OF EVIDENCE\n")
  cat("===========================\n\n")
  
  low_throughput_methods <- c(
    "Luciferase reporter assay", "Reporter assay",
    "Western blot", "Immunoblot",
    "qRT-PCR", "qPCR", "RT-PCR",
    "Northern blot",
    "Immunoprecipitation",
    "ELISA"
  )
  
  high_throughput_methods <- c(
    "Microarray",
    "RNA-Seq", "Sequencing",
    "PAR-CLIP", "HITS-CLIP", "CLASH", "CLIP",
    "Proteomics", "pSILAC"
  )
  
  enhanced_summary <- NULL
  
  if ("experiment" %in% available_cols) {
    # CORRECTED: Only filter out negative evidence and missing IDs
    evidence_classified <- validated_data %>%
      filter(tolower(trimws(support_type)) != "negative" | is.na(support_type)) %>%
      filter(!is.na(target_symbol) & target_symbol != "" & 
             !is.na(target_entrez) & target_entrez != "") %>%
      filter(!is.na(experiment) & experiment != "") %>%
      mutate(
        evidence_class = case_when(
          grepl(paste(low_throughput_methods, collapse = "|"), 
                experiment, ignore.case = TRUE) ~ "Low-throughput",
          grepl(paste(high_throughput_methods, collapse = "|"), 
                experiment, ignore.case = TRUE) ~ "High-throughput",
          TRUE ~ "Other/Unclear"
        )
      )
    
    class_summary <- evidence_classified %>%
      group_by(evidence_class) %>%
      summarise(
        n_records = n(),
        n_targets = n_distinct(target_symbol),
        .groups = "drop"
      )
    
    cat("EVIDENCE CLASSIFICATION:\n")
    print(class_summary)
    cat("\n")
    
    # Create enhanced summary with evidence classification
    enhanced_summary <- evidence_classified %>%
      mutate(
        is_low_throughput = grepl(paste(low_throughput_methods, collapse = "|"), 
                                   experiment, ignore.case = TRUE),
        is_high_throughput = grepl(paste(high_throughput_methods, collapse = "|"), 
                                    experiment, ignore.case = TRUE)
      ) %>%
      group_by(target_symbol, target_entrez) %>%
      summarise(
        n_total_evidence = n(),
        n_low_throughput = sum(is_low_throughput, na.rm = TRUE),
        n_high_throughput = sum(is_high_throughput, na.rm = TRUE),
        n_databases = n_distinct(database),
        databases = paste(unique(database), collapse = "; "),
        experiment_types = paste(unique(experiment[!is.na(experiment) & experiment != ""]), 
                                 collapse = "; "),
        .groups = "drop"
      ) %>%
      arrange(desc(n_total_evidence))
    
    cat("TOP 10 TARGETS WITH DETAILED EVIDENCE:\n")
    cat("(Excluding negative evidence only)\n\n")
    print(head(enhanced_summary, 10), n = 10)
    cat("\n")
  }
  
  # ==============================================================================
  # 5. SAVE FILES
  # ==============================================================================
  
  files_created <- character()
  
  if (save_files) {
    cat("SAVING FILES\n")
    cat("=============\n")
    
    # Detailed evidence file
    f1 <- file.path(output_dir, paste0(output_prefix, "_detailed_validation_evidence.csv"))
    write.csv(target_evidence, f1, row.names = FALSE)
    files_created <- c(files_created, f1)
    cat("✓", basename(f1), "\n")
    
    # Enhanced summary with classification
    if (!is.null(enhanced_summary)) {
      f2 <- file.path(output_dir, paste0(output_prefix, "_enhanced_validation_summary.csv"))
      write.csv(enhanced_summary, f2, row.names = FALSE)
      files_created <- c(files_created, f2)
      cat("✓", basename(f2), "\n")
    }
    
    # Database summary
    f3 <- file.path(output_dir, paste0(output_prefix, "_validation_database_summary.csv"))
    write.csv(db_summary, f3, row.names = FALSE)
    files_created <- c(files_created, f3)
    cat("✓", basename(f3), "\n")
    
    # Experiment summary (if available)
    if (!is.null(exp_summary)) {
      f4 <- file.path(output_dir, paste0(output_prefix, "_validation_experiment_types.csv"))
      write.csv(exp_summary, f4, row.names = FALSE)
      files_created <- c(files_created, f4)
      cat("✓", basename(f4), "\n")
    }
    
    cat("\n")
  }
  
  # ==============================================================================
  # 6. CITATION INFORMATION
  # ==============================================================================
  
  cat("CITATION INFORMATION FOR MANUSCRIPT\n")
  cat("=====================================\n\n")
  
  cat("PRIMARY CITATION (REQUIRED):\n")
  cat("Ru Y, Kechris KJ, Tabakoff B, et al. (2014) The multiMiR R package and database:\n")
  cat("integration of microRNA-target interactions along with their disease and drug associations.\n")
  cat("Nucleic Acids Res, 42(17):e133. DOI: 10.1093/nar/gku631\n\n")
  
  cat("DATABASE-SPECIFIC CITATIONS (cite those present in your data):\n\n")
  
  databases_present <- unique(db_summary$database)
  
  if ("mirtarbase" %in% tolower(databases_present)) {
    cat("miRTarBase:\n")
    cat("Huang HY, Lin YC, Li J, et al. (2020) miRTarBase 2020: updates to the\n")
    cat("experimentally validated microRNA-target interaction database.\n")
    cat("Nucleic Acids Res, 48(D1):D148-D154. DOI: 10.1093/nar/gkz896\n\n")
  }
  
  if ("mirecords" %in% tolower(databases_present)) {
    cat("miRecords:\n")
    cat("Xiao F, Zuo Z, Cai G, et al. (2009) miRecords: an integrated resource\n")
    cat("for microRNA-target interactions.\n")
    cat("Nucleic Acids Res, 37:D105-D110. DOI: 10.1093/nar/gkn851\n\n")
  }
  
  if ("tarbase" %in% tolower(databases_present)) {
    cat("TarBase:\n")
    cat("Karagkouni D, Paraskevopoulou MD, Chatzopoulos S, et al. (2018)\n")
    cat("DIANA-TarBase v8: a decade-long collection of experimentally supported\n")
    cat("miRNA-gene interactions. Nucleic Acids Res, 46(D1):D239-D245.\n")
    cat("DOI: 10.1093/nar/gkx1141\n\n")
  }
  
  cat("SUGGESTED METHODS TEXT:\n")
  cat('"Experimentally validated miRNA-target interactions were retrieved using multiMiR\n')
  cat("(Ru et al., 2014), which integrates validated interactions from miRecords, miRTarBase,\n")
  cat('and TarBase. Only interactions with positive experimental evidence were retained;\n')
  cat('records marked as "negative" (indicating no detectable interaction) were excluded\n')
  cat('from the analysis. The remaining validated targets were confirmed through experimental\n')
  cat('methods including luciferase reporter assays, Western blotting, qRT-PCR, and high-throughput\n')
  cat('techniques such as CLIP-seq."\n\n')
  
  # ==============================================================================
  # RETURN RESULTS
  # ==============================================================================
  
  cat("================================================================================\n")
  cat("                    VALIDATION ANALYSIS COMPLETE                               \n")
  cat("================================================================================\n\n")
  
  return(list(
    target_evidence = target_evidence,
    enhanced_summary = enhanced_summary,
    database_summary = db_summary,
    experiment_summary = exp_summary,
    support_summary = support_summary,
    files_created = files_created,
    databases_present = databases_present,
    filtering_stats = list(
      n_total = n_total_records,
      n_negative = n_negative,
      n_missing_ids = n_missing_ids,
      n_included = nrow(target_evidence)
    )
  ))
}

cat("Validation evidence analysis functions loaded (CORRECTED VERSION).\n")
cat("Usage: validation_results <- analyze_validation_evidence(results, save_files = TRUE)\n\n")
cat("\nCORRECTIONS APPLIED:\n")
cat("  1. Negative evidence excluded from support type summaries\n")
cat("  2. Records with missing target symbols or entrez IDs excluded\n")
cat("  3. REMOVED unnecessary 'functional mti' filter - all positive evidence now included\n")
cat("  4. Updated methods text reflects accurate filtering approach\n\n")
