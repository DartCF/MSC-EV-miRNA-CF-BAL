# ==============================================================================
# MASTER CONFIGURATION: miRNA Pathway Analysis
# ==============================================================================
# This file defines ALL parameters for the 2×2 factorial analysis design.
# It is the SINGLE SOURCE OF TRUTH for:
#   - Disease keywords
#   - Stringency definitions
#   - miRNA lists
#   - Output paths
#   - Analysis metadata
#
# USAGE:
#   config <- load_analysis_config()
#   # or
#   config <- load_analysis_config("path/to/custom_config.yaml")
# ==============================================================================

library(here)

# ==============================================================================
# Default Configuration (as R list)
# ==============================================================================

create_default_config <- function() {
  
  list(
    # =========================================================================
    # METADATA
    # =========================================================================
    metadata = list(
      project_name = "MSC-EV miRNA Pathway Analysis",
      project_id = "CF_EV_miRNA_2024",
      version = "1.0.0",
      created = Sys.time(),
      author = "",
      description = "2×2 factorial analysis varying miRNA scope and target stringency"
    ),
    
    # =========================================================================
    # DISEASE CONTEXT
    # =========================================================================
    disease = list(
      name = "Cystic Fibrosis",
      abbreviation = "CF",
      keywords = c(
        # Inflammation
        "inflamm", "inflammatory", "acute-phase",
        
        # Immune system
        "immune", "immunity", "immunological", "lymphocyte", "leukocyte",
        "T cell", "B cell", "macrophage", "neutrophil", "eosinophil",
        
        # Cytokines & signaling
        "cytokine", "interleukin", "IL-", "TNF", "interferon", "chemokine",
        "NF-kappa", "NFκB", "MAPK", "JAK-STAT",
        
        # Infection & defense
        "infection", "bacterial", "pathogen", "antimicrobial", "antibacterial",
        "antiviral", "toll-like", "TLR", "innate immunity", "adaptive immunity",
        
        # CF-specific
        "CFTR", "cystic fibrosis", "mucus", "mucociliary", "airway",
        
        # Tissue damage & repair
        "fibrosis", "wound", "tissue repair", "extracellular matrix",
        "collagen", "TGF-beta", "epithelial"
      )
    ),
    
    # =========================================================================
    # STRINGENCY DEFINITIONS
    # =========================================================================
    stringency = list(
      # Stringent mode: high-confidence targets
      stringent = list(
        name = "Stringent",
        description = "High-confidence targets only",
        use_all_targets = FALSE,  # Require multiple evidence
        min_databases = 2,        # Minimum prediction databases
        min_publications = 2,     # OR minimum publications
        q_cutoff = 0.05,          # Traditional FDR
        use_case = "Manuscript main figures, validation priorities"
      ),
      
      # Exploratory mode: maximum sensitivity
      exploratory = list(
        name = "Exploratory",
        description = "All predicted targets for maximum sensitivity",
        use_all_targets = TRUE,   # Include all predictions
        min_databases = 1,        # Any single database sufficient
        min_publications = 1,     # Any single publication sufficient
        q_cutoff = 0.20,          # Lenient FDR for discovery
        use_case = "Hypothesis generation, comprehensive coverage"
      )
    ),
    
    # =========================================================================
    # miRNA SETS
    # =========================================================================
    mirna_sets = list(
      # Single miRNA analysis
      single = list(
        name = "miR-6126",
        description = "Most abundant miRNA (33% of reads)",
        mirnas = c("miR-6126"),
        source = "manual",
        source_file = NULL,
        abundance_threshold = NULL
      ),
      
      # Major miRNAs analysis (≥1% abundance)
      major = list(
        name = "Major miRNAs",
        description = "All miRNAs with ≥1% average abundance",
        mirnas = NULL,  # Will be loaded from file
        source = "compositional_analysis",
        source_file = "results/01_compositional/major_miRNAs_1pct_summary.csv",
        abundance_threshold = 1.0  # Percent
      )
    ),
    
    # =========================================================================
    # 2×2 FACTORIAL DESIGN
    # =========================================================================
    analyses = list(
      A = list(
        id = "A",
        name = "miR-6126 Stringent",
        short_name = "6126_str",
        mirna_set = "single",
        stringency_mode = "stringent",
        output_dir = "results/04_pathway/mir6126/stringent",
        output_prefix = "CF_mir6126_stringent"
      ),
      B = list(
        id = "B",
        name = "miR-6126 Exploratory",
        short_name = "6126_exp",
        mirna_set = "single",
        stringency_mode = "exploratory",
        output_dir = "results/04_pathway/mir6126/exploratory",
        output_prefix = "CF_mir6126_exploratory"
      ),
      C = list(
        id = "C",
        name = "Major miRNAs Stringent",
        short_name = "major_str",
        mirna_set = "major",
        stringency_mode = "stringent",
        output_dir = "results/04_pathway/major/stringent",
        output_prefix = "CF_major_stringent"
      ),
      D = list(
        id = "D",
        name = "Major miRNAs Exploratory",
        short_name = "major_exp",
        mirna_set = "major",
        stringency_mode = "exploratory",
        output_dir = "results/04_pathway/major/exploratory",
        output_prefix = "CF_major_exploratory"
      )
    ),
    
    # =========================================================================
    # ENRICHMENT SETTINGS
    # =========================================================================
    enrichment = list(
      min_gene_count = 5,
      max_gene_count = 500,
      databases = c("GO_BP", "GO_MF", "GO_CC", "KEGG", "Reactome"),
      fdr_method = "storey"  # or "BH" for Benjamini-Hochberg
    ),
    
    # =========================================================================
    # OUTPUT SETTINGS
    # =========================================================================
    output = list(
      base_dir = "results/04_pathway",
      comparison_dir = "results/04_pathway/comparison",
      save_rdata = TRUE,
      save_metadata = TRUE,  # Save config with each analysis
      metadata_format = "json"  # or "yaml" or "rds"
    )
  )
}

# ==============================================================================
# Configuration Loading and Validation
# ==============================================================================

load_analysis_config <- function(config_file = NULL) {
  
  # Default config file locations to search (in order of priority)
  # Using here() for portable paths
  default_config_files <- c(
    here("analysis", "04_pathway", "analysis_config.yaml"),
    here("analysis", "04_pathway", "analysis_config.yml"),
    here("analysis_config.yaml"),
    here("config.yaml")
  )
  
  if (is.null(config_file)) {
    # Look for a config file in standard locations
    found_config <- NULL
    for (candidate in default_config_files) {
      if (file.exists(candidate)) {
        found_config <- candidate
        break
      }
    }
    
    if (!is.null(found_config)) {
      config_file <- found_config
      cat(sprintf("Found configuration file: %s\n", config_file))
    } else {
      # No config file found - use hardcoded defaults
      cat("No configuration file found. Using hardcoded defaults.\n")
      cat(sprintf("  (Looked for: %s)\n", paste(basename(default_config_files), collapse = ", ")))
      cat("  To use a config file, create 'analysis_config.yaml' in analysis/04_pathway/.\n\n")
      config <- create_default_config()
      config <- resolve_mirna_lists(config)
      validate_config(config)
      config <- resolve_analysis_parameters(config)
      return(config)
    }
  }
  
  if (endsWith(config_file, ".yaml") || endsWith(config_file, ".yml")) {
    # Load from YAML
    if (!requireNamespace("yaml", quietly = TRUE)) {
      stop("Package 'yaml' required for YAML config files")
    }
    config <- yaml::read_yaml(config_file)
    cat(sprintf("Loaded configuration from: %s\n", config_file))
  } else if (endsWith(config_file, ".rds")) {
    # Load from RDS
    config <- readRDS(config_file)
    cat(sprintf("Loaded configuration from: %s\n", config_file))
  } else {
    stop("Unsupported config file format. Use .yaml, .yml, or .rds")
  }
  
  # Load miRNA lists from source files if specified
  config <- resolve_mirna_lists(config)
  
  # Validate configuration
  validate_config(config)
  
  # Add resolved parameters to each analysis
  config <- resolve_analysis_parameters(config)
  
  return(config)
}

# Resolve miRNA lists from source files
resolve_mirna_lists <- function(config) {
  
  for (set_name in names(config$mirna_sets)) {
    set <- config$mirna_sets[[set_name]]
    
    if (set$source == "compositional_analysis" && !is.null(set$source_file)) {
      # Use here() to resolve the path relative to project root
      source_path <- here(set$source_file)
      
      if (file.exists(source_path)) {
        mirna_table <- read.csv(source_path, stringsAsFactors = FALSE)
        config$mirna_sets[[set_name]]$mirnas <- mirna_table$miRNA
        config$mirna_sets[[set_name]]$mirna_table <- mirna_table
        cat(sprintf("Loaded %d miRNAs from %s\n", 
                   length(mirna_table$miRNA), source_path))
      } else {
        warning(sprintf("miRNA source file not found: %s\n  (resolved to: %s)", 
                       set$source_file, source_path))
      }
    }
  }
  
  return(config)
}

# Validate configuration
validate_config <- function(config) {
  
  errors <- character()
  
  # Check required sections exist
  required_sections <- c("metadata", "disease", "stringency", "mirna_sets", "analyses")
  for (section in required_sections) {
    if (is.null(config[[section]])) {
      errors <- c(errors, sprintf("Missing required section: %s", section))
    }
  }
  
  # Check each analysis references valid miRNA set and stringency mode
  for (analysis_id in names(config$analyses)) {
    analysis <- config$analyses[[analysis_id]]
    
    if (!analysis$mirna_set %in% names(config$mirna_sets)) {
      errors <- c(errors, sprintf("Analysis %s references unknown miRNA set: %s",
                                 analysis_id, analysis$mirna_set))
    }
    
    if (!analysis$stringency_mode %in% names(config$stringency)) {
      errors <- c(errors, sprintf("Analysis %s references unknown stringency mode: %s",
                                 analysis_id, analysis$stringency_mode))
    }
  }
  
  # Check miRNA lists are populated
  for (set_name in names(config$mirna_sets)) {
    if (is.null(config$mirna_sets[[set_name]]$mirnas) || 
        length(config$mirna_sets[[set_name]]$mirnas) == 0) {
      errors <- c(errors, sprintf("miRNA set '%s' has no miRNAs defined", set_name))
    }
  }
  
  if (length(errors) > 0) {
    stop(paste("Configuration validation failed:\n", 
               paste(" -", errors, collapse = "\n")))
  }
  
  cat("✓ Configuration validated successfully\n")
  invisible(TRUE)
}

# Resolve full parameters for each analysis
resolve_analysis_parameters <- function(config) {
  
  for (analysis_id in names(config$analyses)) {
    analysis <- config$analyses[[analysis_id]]
    
    # Get miRNA set parameters
    mirna_set <- config$mirna_sets[[analysis$mirna_set]]
    
    # Get stringency parameters
    stringency <- config$stringency[[analysis$stringency_mode]]
    
    # Resolve output_dir to absolute path using here()
    resolved_output_dir <- here(analysis$output_dir)
    
    # Add resolved parameters
    config$analyses[[analysis_id]]$resolved <- list(
      mirnas = mirna_set$mirnas,
      n_mirnas = length(mirna_set$mirnas),
      mirna_set_name = mirna_set$name,
      use_all_targets = stringency$use_all_targets,
      q_cutoff = stringency$q_cutoff,
      stringency_name = stringency$name,
      disease_keywords = config$disease$keywords
    )
    
    # Update output_dir to resolved path
    config$analyses[[analysis_id]]$output_dir <- resolved_output_dir
  }
  
  # Also resolve comparison_dir
  config$output$comparison_dir <- here(config$output$comparison_dir)
  config$output$base_dir <- here(config$output$base_dir)
  
  return(config)
}

# ==============================================================================
# Save Configuration/Metadata with Results
# ==============================================================================

save_analysis_metadata <- function(config, analysis_id, results, output_dir) {
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Create metadata object
  metadata <- list(
    # Analysis identification
    analysis_id = analysis_id,
    analysis_name = config$analyses[[analysis_id]]$name,
    
    # Timestamp
    run_time = Sys.time(),
    r_version = R.version.string,
    
    # Parameters used (from config)
    parameters = config$analyses[[analysis_id]]$resolved,
    
    # Stringency definition used
    stringency_definition = config$stringency[[config$analyses[[analysis_id]]$stringency_mode]],
    
    # Disease context
    disease = config$disease,
    
    # Results summary
    results_summary = list(
      n_targets = length(results$high_confidence_targets$target_genes),
      n_validated = sum(results$high_confidence_targets$high_confidence_df$evidence_type == "validated", 
                       na.rm = TRUE),
      n_cf_pathways = if(!is.null(results$disease_pathways$disease_pathways)) 
                        nrow(results$disease_pathways$disease_pathways) else 0
    ),
    
    # Full config for complete reproducibility
    full_config = config
  )
  
  # Save based on format preference
  format <- config$output$metadata_format
  
  if (format == "json") {
    if (requireNamespace("jsonlite", quietly = TRUE)) {
      jsonlite::write_json(metadata, 
                          file.path(output_dir, "analysis_metadata.json"),
                          pretty = TRUE, auto_unbox = TRUE)
    } else {
      format <- "rds"  # Fallback
    }
  }
  
  if (format == "yaml") {
    if (requireNamespace("yaml", quietly = TRUE)) {
      yaml::write_yaml(metadata, 
                      file.path(output_dir, "analysis_metadata.yaml"))
    } else {
      format <- "rds"  # Fallback
    }
  }
  
  if (format == "rds") {
    saveRDS(metadata, file.path(output_dir, "analysis_metadata.rds"))
  }
  
  # Also save a simple CSV summary for quick inspection
  summary_df <- data.frame(
    Field = c("Analysis ID", "Analysis Name", "Run Time", 
              "miRNA Set", "N miRNAs", "Stringency Mode",
              "use_all_targets", "q_cutoff",
              "N Targets", "N Validated", "N CF Pathways"),
    Value = c(analysis_id, 
              config$analyses[[analysis_id]]$name,
              as.character(Sys.time()),
              config$analyses[[analysis_id]]$mirna_set,
              metadata$parameters$n_mirnas,
              config$analyses[[analysis_id]]$stringency_mode,
              metadata$parameters$use_all_targets,
              metadata$parameters$q_cutoff,
              metadata$results_summary$n_targets,
              metadata$results_summary$n_validated,
              metadata$results_summary$n_cf_pathways),
    stringsAsFactors = FALSE
  )
  
  write.csv(summary_df, file.path(output_dir, "analysis_summary.csv"), row.names = FALSE)
  
  cat(sprintf("✓ Metadata saved to: %s\n", output_dir))
  
  invisible(metadata)
}

# ==============================================================================
# Load Metadata from Results (for comparison)
# ==============================================================================

load_analysis_metadata <- function(results_dir) {
  
  # Try JSON first
  json_file <- file.path(results_dir, "analysis_metadata.json")
  if (file.exists(json_file) && requireNamespace("jsonlite", quietly = TRUE)) {
    return(jsonlite::read_json(json_file))
  }

  # Try YAML
  yaml_file <- file.path(results_dir, "analysis_metadata.yaml")
  if (file.exists(yaml_file) && requireNamespace("yaml", quietly = TRUE)) {
    return(yaml::read_yaml(yaml_file))
  }
  
  # Try RDS
  rds_file <- file.path(results_dir, "analysis_metadata.rds")
  if (file.exists(rds_file)) {
    return(readRDS(rds_file))
  }
  
  # Try CSV summary as fallback
  csv_file <- file.path(results_dir, "analysis_summary.csv")
  if (file.exists(csv_file)) {
    df <- read.csv(csv_file, stringsAsFactors = FALSE)
    # Convert to list
    metadata <- as.list(setNames(df$Value, df$Field))
    metadata$source <- "csv_summary"
    return(metadata)
  }
  
  warning(sprintf("No metadata found in: %s", results_dir))
  return(NULL)
}

# ==============================================================================
# Print Configuration Summary
# ==============================================================================

print_config_summary <- function(config) {
  
  cat("\n")
  cat("================================================================================\n")
  cat("                    ANALYSIS CONFIGURATION SUMMARY                             \n")
  cat("================================================================================\n\n")
  
  cat("PROJECT:\n")
  cat(sprintf("  Name: %s\n", config$metadata$project_name))
  cat(sprintf("  Version: %s\n", config$metadata$version))
  cat("\n")
  
  cat("DISEASE CONTEXT:\n")
  cat(sprintf("  Disease: %s (%s)\n", config$disease$name, config$disease$abbreviation))
  cat(sprintf("  Keywords: %d terms defined\n", length(config$disease$keywords)))
  cat("\n")
  
  cat("STRINGENCY DEFINITIONS:\n")
  for (mode in names(config$stringency)) {
    s <- config$stringency[[mode]]
    cat(sprintf("  %s:\n", toupper(s$name)))
    cat(sprintf("    use_all_targets: %s\n", s$use_all_targets))
    cat(sprintf("    q_cutoff: %.2f\n", s$q_cutoff))
    cat(sprintf("    Use case: %s\n", s$use_case))
  }
  cat("\n")
  
  cat("miRNA SETS:\n")
  for (set_name in names(config$mirna_sets)) {
    set <- config$mirna_sets[[set_name]]
    cat(sprintf("  %s: %d miRNAs\n", set$name, length(set$mirnas)))
  }
  cat("\n")
  
  cat("2×2 FACTORIAL DESIGN:\n")
  cat("─────────────────────────────────────────────────────────────────────\n")
  cat("                          Evidence Threshold                         \n")
  cat("                    Stringent          Exploratory                   \n")
  cat("─────────────────────────────────────────────────────────────────────\n")
  
  A <- config$analyses$A$resolved
  B <- config$analyses$B$resolved
  C <- config$analyses$C$resolved
  D <- config$analyses$D$resolved
  
  cat(sprintf("%-13s    [A] %d miRNA(s)    [B] %d miRNA(s)\n",
              config$mirna_sets$single$name,
              A$n_mirnas, B$n_mirnas))
  cat(sprintf("%-13s    q < %.2f           q < %.2f\n", "", A$q_cutoff, B$q_cutoff))
  cat(sprintf("%-13s    [C] %d miRNA(s)    [D] %d miRNA(s)\n",
              config$mirna_sets$major$name,
              C$n_mirnas, D$n_mirnas))
  cat(sprintf("%-13s    q < %.2f           q < %.2f\n", "", C$q_cutoff, D$q_cutoff))
  cat("─────────────────────────────────────────────────────────────────────\n")
  cat("\n")
  
  cat("OUTPUT:\n")
  for (analysis_id in names(config$analyses)) {
    a <- config$analyses[[analysis_id]]
    cat(sprintf("  [%s] %s → %s\n", analysis_id, a$name, a$output_dir))
  }
  cat("\n")
  
  cat("================================================================================\n\n")
}

# ==============================================================================
# Export Configuration to YAML
# ==============================================================================

export_config_yaml <- function(config, output_file = "analysis_config.yaml") {
  
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' required. Install with: install.packages('yaml')")
  }
  
  # Remove resolved parameters (they're derived, not source)
  export_config <- config
  for (analysis_id in names(export_config$analyses)) {
    export_config$analyses[[analysis_id]]$resolved <- NULL
  }
  
  # Remove loaded miRNA tables (keep just the file reference)
  for (set_name in names(export_config$mirna_sets)) {
    export_config$mirna_sets[[set_name]]$mirna_table <- NULL
    # Keep mirnas for single set, remove for file-based
    if (export_config$mirna_sets[[set_name]]$source == "compositional_analysis") {
      export_config$mirna_sets[[set_name]]$mirnas <- NULL
    }
  }
  
  yaml::write_yaml(export_config, output_file)
  cat(sprintf("✓ Configuration exported to: %s\n", output_file))
}

# ==============================================================================
# Example Usage
# ==============================================================================

cat("\n")
cat("================================================================================\n")
cat("              Analysis Configuration System Loaded                              \n")
cat("================================================================================\n\n")

cat("CONFIGURATION LOADING BEHAVIOR:\n")
cat("-------------------------------\n")
cat("When you call load_analysis_config() with no arguments, it searches for:\n")
cat("  1. analysis/04_pathway/analysis_config.yaml\n")
cat("  2. analysis/04_pathway/analysis_config.yml\n")
cat("  3. analysis_config.yaml (project root)\n")
cat("  4. config.yaml (project root)\n")
cat("If none found, falls back to hardcoded defaults.\n\n")

cat("USAGE:\n")
cat("------\n")
cat("# Load configuration (auto-detects analysis_config.yaml)\n")
cat("config <- load_analysis_config()\n\n")

cat("# Print summary\n")
cat("print_config_summary(config)\n\n")

cat("# Access parameters for a specific analysis\n")
cat("config$analyses$A$resolved$mirnas      # miRNAs for analysis A\n")
cat("config$analyses$A$resolved$q_cutoff    # FDR threshold for analysis A\n")
cat("config$disease$keywords                # Disease keywords\n\n")

cat("# Export to YAML (if you modified the hardcoded defaults)\n")
cat("export_config_yaml(config, 'my_analysis_config.yaml')\n\n")

cat("# Load from explicit path\n")
cat("config <- load_analysis_config('path/to/my_config.yaml')\n\n")

cat("================================================================================\n\n")
