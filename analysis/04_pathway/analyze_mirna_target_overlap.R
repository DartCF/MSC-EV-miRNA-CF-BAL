# ==============================================================================
# miRNA Target Overlap Analysis
# ==============================================================================
# Analyzes overlap between miRNA target sets using Jaccard index
# Identifies unique targets for each miRNA
# Creates publication-ready visualizations
#
# USAGE:
#   source(here("analysis", "04_pathway", "analyze_mirna_target_overlap.R"))
#   
#   # Basic usage (auto-detects input file):
#   results <- analyze_mirna_overlap()
#   
#   # With custom input file:
#   results <- analyze_mirna_overlap("path/to/targets.csv")
#   
#   # Focus on a specific miRNA:
#   results <- analyze_mirna_overlap(focus_mirna = "hsa-miR-21-5p")
#   
#   # Access results:
#   results$jaccard_matrix        # Pairwise Jaccard indices
#   results$unique_summary        # Unique targets per miRNA
#   results$unique_targets        # List of unique target genes
#   results$focus_analysis        # Detailed analysis of focus miRNA
#
# INPUT:
#   results/04_pathway/major/stringent/CF_major_stringent_high_confidence_targets.csv
#
# OUTPUT:
#   results/04_pathway/major/stringent/mirna_overlap_summary.png (Fig 4)
#   results/04_pathway/major/stringent/mirna_jaccard_heatmap.png
#   results/04_pathway/major/stringent/mirna_unique_vs_shared.png
#   results/04_pathway/major/stringent/mirna_pct_unique.png
#   results/04_pathway/major/stringent/*.csv (various data tables)
#
# ==============================================================================

library(here)
library(dplyr)
library(tidyr)
library(ggplot2)

# ==============================================================================
# CORE FUNCTIONS
# ==============================================================================

#' Calculate Jaccard Index between two sets
#' 
#' @param set1 Character vector of items in first set
#' @param set2 Character vector of items in second set
#' @return Numeric Jaccard index (0-1)
calc_jaccard <- function(set1, set2) {
  set1 <- unique(set1)
  set2 <- unique(set2)
  
  if (length(set1) == 0 && length(set2) == 0) return(NA_real_)
  
  intersection <- length(intersect(set1, set2))
  union <- length(union(set1, set2))
  
  if (union == 0) return(0)
  return(intersection / union)
}

#' Calculate overlap coefficient (intersection / smaller set)
#' 
#' @param set1 Character vector of items in first set
#' @param set2 Character vector of items in second set
#' @return Numeric overlap coefficient (0-1)
calc_overlap_coef <- function(set1, set2) {
  set1 <- unique(set1)
  set2 <- unique(set2)
  
  if (length(set1) == 0 || length(set2) == 0) return(0)
  
  intersection <- length(intersect(set1, set2))
  smaller <- min(length(set1), length(set2))
  
  return(intersection / smaller)
}

#' Extract target sets for each miRNA from a data frame
#' 
#' @param df Data frame with columns: mature_mirna_id, target_symbol
#' @param mirna_col Name of miRNA column (default: "mature_mirna_id")
#' @param target_col Name of target column (default: "target_symbol")
#' @return Named list of character vectors (target sets per miRNA)
get_mirna_target_sets <- function(df, 
                                   mirna_col = "mature_mirna_id",
                                   target_col = "target_symbol") {
  
  # Remove empty/NA targets
 df_clean <- df %>%
    filter(!is.na(.data[[target_col]]) & .data[[target_col]] != "")
  
  # Get unique miRNAs
  mirnas <- unique(df_clean[[mirna_col]])
  
  # Build list of target sets
  target_sets <- lapply(mirnas, function(m) {
    targets <- df_clean %>%
      filter(.data[[mirna_col]] == m) %>%
      pull(.data[[target_col]]) %>%
      unique()
    return(targets)
  })
  names(target_sets) <- mirnas
  
  return(target_sets)
}

#' Calculate pairwise Jaccard matrix for all miRNA target sets
#' 
#' @param target_sets Named list of target sets (from get_mirna_target_sets)
#' @return List with jaccard_matrix, intersection_matrix, and union_matrix
calc_pairwise_jaccard <- function(target_sets) {
  
  mirnas <- names(target_sets)
  n <- length(mirnas)
  
  # Initialize matrices
  jaccard_mat <- matrix(NA_real_, nrow = n, ncol = n,
                        dimnames = list(mirnas, mirnas))
  intersection_mat <- matrix(0L, nrow = n, ncol = n,
                             dimnames = list(mirnas, mirnas))
  union_mat <- matrix(0L, nrow = n, ncol = n,
                      dimnames = list(mirnas, mirnas))
  
  # Calculate pairwise metrics
 for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      set_i <- target_sets[[mirnas[i]]]
      set_j <- target_sets[[mirnas[j]]]
      
      jaccard_mat[i, j] <- calc_jaccard(set_i, set_j)
      intersection_mat[i, j] <- length(intersect(set_i, set_j))
      union_mat[i, j] <- length(union(set_i, set_j))
    }
  }
  
  return(list(
    jaccard = jaccard_mat,
    intersection = intersection_mat,
    union = union_mat
  ))
}

#' Find unique targets for each miRNA (not targeted by any other)
#' 
#' @param target_sets Named list of target sets
#' @return List with summary data frame and detailed target lists
find_unique_targets <- function(target_sets) {
  
  mirnas <- names(target_sets)
  all_targets <- unique(unlist(target_sets))
  
  unique_targets <- list()
  summary_data <- list()
  
  for (mirna in mirnas) {
    # Get targets of this miRNA
    this_targets <- target_sets[[mirna]]
    
    # Get targets of all OTHER miRNAs
    other_mirnas <- setdiff(mirnas, mirna)
    other_targets <- unique(unlist(target_sets[other_mirnas]))
    
    # Find unique targets (in this, not in others)
    unique_to_this <- setdiff(this_targets, other_targets)
    unique_targets[[mirna]] <- sort(unique_to_this)
    
    # Summary stats
    summary_data[[mirna]] <- data.frame(
      miRNA = mirna,
      miRNA_short = gsub("hsa-", "", mirna),
      Total_Targets = length(this_targets),
      Unique_Targets = length(unique_to_this),
      Shared_Targets = length(this_targets) - length(unique_to_this),
      Pct_Unique = 100 * length(unique_to_this) / length(this_targets),
      stringsAsFactors = FALSE
    )
  }
  
  summary_df <- bind_rows(summary_data) %>%
    arrange(desc(Unique_Targets))
  
  return(list(
    summary = summary_df,
    unique_targets = unique_targets
  ))
}

#' Analyze overlap for a specific focus miRNA
#' 
#' @param target_sets Named list of target sets
#' @param focus_mirna Name of miRNA to focus on
#' @return List with overlap details
analyze_focus_mirna <- function(target_sets, focus_mirna) {
  
  if (!focus_mirna %in% names(target_sets)) {
    warning(sprintf("Focus miRNA '%s' not found in data", focus_mirna))
    return(NULL)
  }
  
  focus_targets <- target_sets[[focus_mirna]]
  other_mirnas <- setdiff(names(target_sets), focus_mirna)
  
  # Calculate overlap with each other miRNA
  overlaps <- lapply(other_mirnas, function(other) {
    other_targets <- target_sets[[other]]
    shared <- intersect(focus_targets, other_targets)
    
    data.frame(
      miRNA = other,
      miRNA_short = gsub("hsa-", "", other),
      Shared_Targets = length(shared),
      Pct_of_Focus = 100 * length(shared) / length(focus_targets),
      Jaccard = calc_jaccard(focus_targets, other_targets),
      Shared_Genes = paste(sort(shared), collapse = ", "),
      stringsAsFactors = FALSE
    )
  })
  
  overlap_df <- bind_rows(overlaps) %>%
    arrange(desc(Shared_Targets))
  
  # Get unique targets for focus miRNA
  all_other_targets <- unique(unlist(target_sets[other_mirnas]))
  unique_to_focus <- setdiff(focus_targets, all_other_targets)
  shared_with_any <- intersect(focus_targets, all_other_targets)
  
  return(list(
    focus_mirna = focus_mirna,
    n_total = length(focus_targets),
    n_unique = length(unique_to_focus),
    n_shared = length(shared_with_any),
    pct_unique = 100 * length(unique_to_focus) / length(focus_targets),
    unique_genes = sort(unique_to_focus),
    shared_genes = sort(shared_with_any),
    overlap_with_others = overlap_df
  ))
}

# ==============================================================================
# VISUALIZATION FUNCTIONS
# ==============================================================================

#' Create Jaccard similarity heatmap
#' 
#' @param jaccard_matrix Jaccard matrix from calc_pairwise_jaccard
#' @param title Plot title
#' @param max_scale Maximum value for color scale (default: auto)
#' @return ggplot object
plot_jaccard_heatmap <- function(jaccard_matrix, 
                                  title = "miRNA Target Similarity (Jaccard Index)",
                                  max_scale = NULL) {
  
  # Convert to long format
  mirnas <- rownames(jaccard_matrix)
  short_names <- gsub("hsa-", "", mirnas)
  
  jaccard_long <- expand.grid(
    miRNA1 = short_names,
    miRNA2 = short_names,
    stringsAsFactors = FALSE
  )
  jaccard_long$Jaccard <- as.vector(jaccard_matrix)
  
  # Set factor levels for ordering
  jaccard_long$miRNA1 <- factor(jaccard_long$miRNA1, levels = short_names)
  jaccard_long$miRNA2 <- factor(jaccard_long$miRNA2, levels = rev(short_names))
  
  # Determine scale maximum
  if (is.null(max_scale)) {
    off_diag <- jaccard_matrix[row(jaccard_matrix) != col(jaccard_matrix)]
    max_scale <- max(off_diag, na.rm = TRUE) * 1.2
    max_scale <- max(max_scale, 0.1)  # Minimum of 0.1
  }
  
  p <- ggplot(jaccard_long, aes(x = miRNA1, y = miRNA2, fill = Jaccard)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f", Jaccard)), size = 3) +
    scale_fill_gradient2(low = "white", mid = "#FFEDA0", high = "#E31A1C",
                        midpoint = max_scale / 2,
                        limits = c(0, max_scale),
                        name = "Jaccard\nIndex") +
    labs(title = title, x = "", y = "") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      panel.grid = element_blank()
    ) +
    coord_fixed()
  
  return(p)
}

#' Create unique targets bar plot
#' 
#' @param unique_summary Summary data frame from find_unique_targets
#' @param highlight_mirna miRNA to highlight (optional)
#' @return ggplot object
plot_unique_targets <- function(unique_summary, highlight_mirna = NULL) {
  
  # Prepare data
  plot_data <- unique_summary %>%
    arrange(Total_Targets) %>%
    mutate(miRNA_short = factor(miRNA_short, levels = miRNA_short))
  
  # Reshape for stacked bar
  plot_long <- plot_data %>%
    pivot_longer(cols = c(Unique_Targets, Shared_Targets),
                 names_to = "Type",
                 values_to = "Count") %>%
    mutate(Type = ifelse(Type == "Unique_Targets", "Unique", "Shared"))
  
  p <- ggplot(plot_long, aes(x = miRNA_short, y = Count, fill = Type)) +
    geom_bar(stat = "identity", alpha = 0.85) +
    coord_flip() +
    scale_fill_manual(values = c("Unique" = "#2E86AB", "Shared" = "#A23B72"),
                      name = "Target Type") +
    labs(title = "Unique vs Shared Targets per miRNA",
         x = "miRNA",
         y = "Number of Target Genes") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.text.y = element_text(size = 10),
      panel.grid.major.y = element_blank(),
      legend.position = "bottom"
    )
  
  return(p)
}

#' Create percent unique bar plot
#' 
#' @param unique_summary Summary data frame from find_unique_targets
#' @param highlight_mirna miRNA to highlight (optional)
#' @return ggplot object
plot_pct_unique <- function(unique_summary, highlight_mirna = NULL) {
  
  plot_data <- unique_summary %>%
    arrange(Pct_Unique) %>%
    mutate(miRNA_short = factor(miRNA_short, levels = miRNA_short))
  
  # Set highlight color
  if (!is.null(highlight_mirna)) {
    highlight_short <- gsub("hsa-", "", highlight_mirna)
    plot_data$Highlight <- ifelse(plot_data$miRNA_short == highlight_short, 
                                  "Focus", "Other")
  } else {
    plot_data$Highlight <- "Other"
  }
  
  p <- ggplot(plot_data, aes(x = miRNA_short, y = Pct_Unique, fill = Highlight)) +
    geom_bar(stat = "identity", alpha = 0.85) +
    geom_hline(yintercept = 50, linetype = "dashed", color = "red", alpha = 0.7) +
    geom_text(aes(label = sprintf("%.1f%%", Pct_Unique)), 
              hjust = -0.1, size = 3.5) +
    coord_flip() +
    scale_fill_manual(values = c("Focus" = "#EFC000", "Other" = "#0073C2"),
                      guide = "none") +
    scale_y_continuous(limits = c(0, 105), expand = c(0, 0)) +
    labs(title = "Percentage of Unique Targets",
         subtitle = if (!is.null(highlight_mirna)) 
           sprintf("%s highlighted", gsub("hsa-", "", highlight_mirna)) else NULL,
         x = "miRNA",
         y = "% of Targets that are Unique") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10, color = "gray40"),
      axis.text.y = element_text(size = 10),
      panel.grid.major.y = element_blank()
    )
  
  return(p)
}

#' Create focus miRNA overlap bar plot
#' 
#' @param focus_analysis Results from analyze_focus_mirna
#' @return ggplot object
plot_focus_overlap <- function(focus_analysis) {
  
  if (is.null(focus_analysis)) return(NULL)
  
  plot_data <- focus_analysis$overlap_with_others %>%
    filter(Shared_Targets > 0) %>%
    arrange(Shared_Targets) %>%
    mutate(miRNA_short = factor(miRNA_short, levels = miRNA_short))
  
  if (nrow(plot_data) == 0) {
    return(NULL)
  }
  
  focus_short <- gsub("hsa-", "", focus_analysis$focus_mirna)
  
  p <- ggplot(plot_data, aes(x = miRNA_short, y = Shared_Targets)) +
    geom_bar(stat = "identity", fill = "#E74C3C", alpha = 0.85) +
    geom_text(aes(label = Shared_Targets), hjust = -0.2, size = 4) +
    coord_flip() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(title = sprintf("%s Target Overlap with Other miRNAs", focus_short),
         subtitle = sprintf("%d/%d targets shared (%.1f%% unique)",
                           focus_analysis$n_shared,
                           focus_analysis$n_total,
                           focus_analysis$pct_unique),
         x = "miRNA",
         y = "Number of Shared Targets") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10, color = "gray40"),
      axis.text.y = element_text(size = 10),
      panel.grid.major.y = element_blank()
    )
  
  return(p)
}

#' Create combined summary figure (requires patchwork)
#' 
#' @param results Results list from analyze_mirna_overlap
#' @return Combined ggplot object
plot_summary_figure <- function(results) {
  
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    warning("Package 'patchwork' required for combined figure")
    return(NULL)
  }
  
  library(patchwork)
  
  p1 <- plot_jaccard_heatmap(results$jaccard_matrix)
  p2 <- plot_unique_targets(results$unique_summary)
  p3 <- plot_pct_unique(results$unique_summary, results$focus_mirna)
  p4 <- plot_focus_overlap(results$focus_analysis)
  
  if (is.null(p4)) {
    combined <- (p1 | p2) / p3 +
      plot_annotation(title = "miRNA Target Overlap Analysis Summary")
  } else {
    combined <- (p1 | p2) / (p3 | p4) +
      plot_annotation(title = "miRNA Target Overlap Analysis Summary")
  }
  
  return(combined)
}

# ==============================================================================
# MAIN ANALYSIS FUNCTION
# ==============================================================================

#' Run complete miRNA target overlap analysis
#' 
#' @param input_file Path to CSV file with miRNA-target data
#' @param output_dir Directory for output files (default: same as input)
#' @param mirna_col Name of miRNA column (default: "mature_mirna_id")
#' @param target_col Name of target column (default: "target_symbol")
#' @param focus_mirna miRNA to analyze in detail (default: "hsa-miR-6126")
#' @param save_outputs Save CSV and plot files (default: TRUE)
#' @param verbose Print progress messages (default: TRUE)
#' @return List with all analysis results
analyze_mirna_overlap <- function(input_file = NULL,
                                   output_dir = NULL,
                                   mirna_col = "mature_mirna_id",
                                   target_col = "target_symbol",
                                   focus_mirna = "hsa-miR-6126",
                                   save_outputs = TRUE,
                                   verbose = TRUE) {
  
  # --------------------------------------------------------------------------
  # Setup
  # --------------------------------------------------------------------------
  
  if (verbose) {
    cat("\n")
    cat("================================================================================\n")
    cat("              miRNA TARGET OVERLAP ANALYSIS                                     \n")
    cat("================================================================================\n\n")
  }
  
  # Auto-detect input file if not provided
  if (is.null(input_file)) {
    default_files <- c(
      here("results", "04_pathway", "major", "stringent", 
           "CF_major_stringent_high_confidence_targets.csv"),
      "CF_major_stringent_high_confidence_targets.csv"
    )
    for (f in default_files) {
      if (file.exists(f)) {
        input_file <- f
        break
      }
    }
    if (is.null(input_file)) {
      stop("No input file specified and no default file found.\n",
           "Expected: results/04_pathway/major/stringent/CF_major_stringent_high_confidence_targets.csv")
    }
  }
  
  if (!file.exists(input_file)) {
    stop(sprintf("Input file not found: %s", input_file))
  }
  
  # Set output directory
  if (is.null(output_dir)) {
    output_dir <- dirname(input_file)
  }
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  if (verbose) {
    cat(sprintf("Input file: %s\n", input_file))
    cat(sprintf("Output directory: %s\n\n", output_dir))
  }
  
  # --------------------------------------------------------------------------
  # Load and clean data
  # --------------------------------------------------------------------------
  
  df <- read.csv(input_file, stringsAsFactors = FALSE)
  
  # Check for required columns
  if (!mirna_col %in% colnames(df)) {
    stop(sprintf("Column '%s' not found in data", mirna_col))
  }
  if (!target_col %in% colnames(df)) {
    stop(sprintf("Column '%s' not found in data", target_col))
  }
  
  # Report data quality
  n_total <- nrow(df)
  n_empty <- sum(df[[target_col]] == "" | is.na(df[[target_col]]))
  
  if (verbose) {
    cat("DATA QUALITY:\n")
    cat("-------------\n")
    cat(sprintf("Total rows: %d\n", n_total))
    cat(sprintf("Rows with empty target: %d (excluded)\n", n_empty))
    cat(sprintf("Valid rows: %d\n\n", n_total - n_empty))
  }
  
  # Get target sets
  target_sets <- get_mirna_target_sets(df, mirna_col, target_col)
  n_mirnas <- length(target_sets)
  n_unique_targets <- length(unique(unlist(target_sets)))
  
  if (verbose) {
    cat(sprintf("Number of miRNAs: %d\n", n_mirnas))
    cat(sprintf("Total unique targets: %d\n\n", n_unique_targets))
  }
  
  # --------------------------------------------------------------------------
  # Calculate Jaccard matrix
  # --------------------------------------------------------------------------
  
  if (verbose) cat("Calculating pairwise Jaccard indices...\n")
  
  pairwise_results <- calc_pairwise_jaccard(target_sets)
  jaccard_matrix <- pairwise_results$jaccard
  
  if (verbose) {
    cat("✓ Jaccard matrix calculated\n\n")
    
    # Print summary
    cat("JACCARD SIMILARITY SUMMARY:\n")
    cat("---------------------------\n")
    off_diag <- jaccard_matrix[row(jaccard_matrix) != col(jaccard_matrix)]
    cat(sprintf("Range (off-diagonal): %.4f - %.4f\n", 
                min(off_diag, na.rm = TRUE), max(off_diag, na.rm = TRUE)))
    cat(sprintf("Mean (off-diagonal): %.4f\n", mean(off_diag, na.rm = TRUE)))
    cat(sprintf("Median (off-diagonal): %.4f\n\n", median(off_diag, na.rm = TRUE)))
  }
  
  # --------------------------------------------------------------------------
  # Find unique targets
  # --------------------------------------------------------------------------
  
  if (verbose) cat("Finding unique targets per miRNA...\n")
  
  unique_results <- find_unique_targets(target_sets)
  
  if (verbose) {
    cat("✓ Unique targets identified\n\n")
    cat("UNIQUE TARGETS SUMMARY:\n")
    cat("-----------------------\n")
    print(unique_results$summary[, c("miRNA_short", "Total_Targets", 
                                      "Unique_Targets", "Pct_Unique")], 
          row.names = FALSE, digits = 1)
    cat("\n")
  }
  
  # --------------------------------------------------------------------------
  # Focus miRNA analysis
  # --------------------------------------------------------------------------
  
  focus_analysis <- NULL
  if (!is.null(focus_mirna)) {
    if (verbose) cat(sprintf("Analyzing focus miRNA: %s\n", focus_mirna))
    
    focus_analysis <- analyze_focus_mirna(target_sets, focus_mirna)
    
    if (!is.null(focus_analysis) && verbose) {
      cat("✓ Focus analysis complete\n\n")
      cat(sprintf("FOCUS miRNA: %s\n", focus_mirna))
      cat("---------------------------\n")
      cat(sprintf("Total targets: %d\n", focus_analysis$n_total))
      cat(sprintf("Unique targets: %d (%.1f%%)\n",
                  focus_analysis$n_unique, focus_analysis$pct_unique))
      cat(sprintf("Shared targets: %d\n\n", focus_analysis$n_shared))
      
      if (focus_analysis$n_unique > 0) {
        cat("Unique targets:\n")
        genes <- focus_analysis$unique_genes
        if (length(genes) <= 20) {
          cat(paste(genes, collapse = ", "), "\n")
        } else {
          cat(paste(genes[1:20], collapse = ", "), 
              sprintf(", ... and %d more\n", length(genes) - 20))
        }
      } else {
        cat("⚠ No unique targets - all targets shared with other miRNAs\n")
      }
      cat("\n")
    }
  }
  
  # --------------------------------------------------------------------------
  # Save outputs
  # --------------------------------------------------------------------------
  
  if (save_outputs) {
    if (verbose) cat("Saving outputs...\n")
    
    # Jaccard matrix
    jaccard_df <- as.data.frame(jaccard_matrix)
    jaccard_df$miRNA <- rownames(jaccard_df)
    jaccard_df <- jaccard_df[, c("miRNA", setdiff(names(jaccard_df), "miRNA"))]
    write.csv(jaccard_df, 
              file.path(output_dir, "mirna_jaccard_matrix.csv"),
              row.names = FALSE)
    if (verbose) cat("  ✓ mirna_jaccard_matrix.csv\n")
    
    # Unique targets summary
    write.csv(unique_results$summary,
              file.path(output_dir, "mirna_unique_targets_summary.csv"),
              row.names = FALSE)
    if (verbose) cat("  ✓ mirna_unique_targets_summary.csv\n")
    
    # Detailed unique target lists
    unique_list <- lapply(names(unique_results$unique_targets), function(m) {
      targets <- unique_results$unique_targets[[m]]
      if (length(targets) > 0) {
        data.frame(
          miRNA = m,
          miRNA_short = gsub("hsa-", "", m),
          Unique_Target = targets,
          stringsAsFactors = FALSE
        )
      } else {
        NULL
      }
    })
    unique_list_df <- bind_rows(unique_list)
    write.csv(unique_list_df,
              file.path(output_dir, "mirna_unique_targets_detailed.csv"),
              row.names = FALSE)
    if (verbose) cat("  ✓ mirna_unique_targets_detailed.csv\n")
    
    # Focus miRNA analysis
    if (!is.null(focus_analysis)) {
      focus_short <- gsub("hsa-", "", focus_mirna)
      write.csv(focus_analysis$overlap_with_others,
                file.path(output_dir, sprintf("%s_overlap_with_others.csv", focus_short)),
                row.names = FALSE)
      if (verbose) cat(sprintf("  ✓ %s_overlap_with_others.csv\n", focus_short))
    }
    
    # Create plots
    if (verbose) cat("\nCreating visualizations...\n")
    
    # Jaccard heatmap
    p1 <- plot_jaccard_heatmap(jaccard_matrix)
    ggsave(file.path(output_dir, "mirna_jaccard_heatmap.pdf"), p1,
           width = 10, height = 9, dpi = 300)
    ggsave(file.path(output_dir, "mirna_jaccard_heatmap.png"), p1,
           width = 10, height = 9, dpi = 300)
    if (verbose) cat("  ✓ mirna_jaccard_heatmap.pdf/png\n")
    
    # Unique targets plot
    p2 <- plot_unique_targets(unique_results$summary)
    ggsave(file.path(output_dir, "mirna_unique_vs_shared.pdf"), p2,
           width = 10, height = 8, dpi = 300)
    ggsave(file.path(output_dir, "mirna_unique_vs_shared.png"), p2,
           width = 10, height = 8, dpi = 300)
    if (verbose) cat("  ✓ mirna_unique_vs_shared.pdf/png\n")
    
    # Percent unique plot
    p3 <- plot_pct_unique(unique_results$summary, focus_mirna)
    ggsave(file.path(output_dir, "mirna_pct_unique.pdf"), p3,
           width = 10, height = 8, dpi = 300)
    ggsave(file.path(output_dir, "mirna_pct_unique.png"), p3,
           width = 10, height = 8, dpi = 300)
    if (verbose) cat("  ✓ mirna_pct_unique.pdf/png\n")
    
    # Focus miRNA overlap plot
    if (!is.null(focus_analysis)) {
      p4 <- plot_focus_overlap(focus_analysis)
      if (!is.null(p4)) {
        focus_short <- gsub("hsa-", "", focus_mirna)
        ggsave(file.path(output_dir, sprintf("%s_overlap_barplot.pdf", focus_short)), p4,
               width = 9, height = 6, dpi = 300)
        ggsave(file.path(output_dir, sprintf("%s_overlap_barplot.png", focus_short)), p4,
               width = 9, height = 6, dpi = 300)
        if (verbose) cat(sprintf("  ✓ %s_overlap_barplot.pdf/png\n", focus_short))
      }
    }
    
    # Multi-panel summary (if patchwork available)
    if (requireNamespace("patchwork", quietly = TRUE)) {
      if (verbose) cat("  Creating combined summary figure...\n")
      tryCatch({
        p_combined <- plot_summary_figure(list(
          jaccard_matrix = jaccard_matrix,
          unique_summary = unique_results$summary,
          focus_mirna = focus_mirna,
          focus_analysis = focus_analysis
        ))
        ggsave(file.path(output_dir, "mirna_overlap_summary.pdf"), p_combined,
               width = 16, height = 14, dpi = 300)
        ggsave(file.path(output_dir, "mirna_overlap_summary.png"), p_combined,
               width = 16, height = 14, dpi = 300)
        if (verbose) cat("  ✓ mirna_overlap_summary.pdf/png\n")
      }, error = function(e) {
        if (verbose) cat("  (Could not create combined figure)\n")
      })
    }
  }
  
  # --------------------------------------------------------------------------
  # Return results
  # --------------------------------------------------------------------------
  
  if (verbose) {
    cat("\n")
    cat("================================================================================\n")
    cat("                         ANALYSIS COMPLETE                                      \n")
    cat("================================================================================\n\n")
  }
  
  results <- list(
    # Input info
    input_file = input_file,
    n_mirnas = n_mirnas,
    n_unique_targets = n_unique_targets,
    
    # Target sets
    target_sets = target_sets,
    
    # Jaccard analysis
    jaccard_matrix = jaccard_matrix,
    intersection_matrix = pairwise_results$intersection,
    union_matrix = pairwise_results$union,
    
    # Unique targets
    unique_summary = unique_results$summary,
    unique_targets = unique_results$unique_targets,
    
    # Focus miRNA
    focus_mirna = focus_mirna,
    focus_analysis = focus_analysis
  )
  
  class(results) <- c("mirna_overlap_analysis", class(results))
  
  invisible(results)
}

#' Print method for mirna_overlap_analysis objects
#' 
#' @param x mirna_overlap_analysis object
#' @param ... Additional arguments (ignored)
print.mirna_overlap_analysis <- function(x, ...) {
  cat("miRNA Target Overlap Analysis\n")
  cat("=============================\n")
  cat(sprintf("miRNAs: %d\n", x$n_mirnas))
  cat(sprintf("Unique targets: %d\n", x$n_unique_targets))
  if (!is.null(x$focus_mirna)) {
    cat(sprintf("Focus miRNA: %s\n", x$focus_mirna))
  }
  cat("\nUse $ to access components: jaccard_matrix, unique_summary, unique_targets, focus_analysis\n")
}

# ==============================================================================
# AUTO-RUN WHEN SOURCED
# ==============================================================================

cat("\n")
cat("================================================================================\n")
cat("           miRNA Target Overlap Analysis Script Loaded                          \n")
cat("================================================================================\n\n")

cat("Running analysis with default settings...\n\n")
results <- analyze_mirna_overlap()

cat("\nUSAGE:\n")
cat("------\n")
cat("# Access results:\n")
cat("results$jaccard_matrix        # Pairwise Jaccard indices\n")
cat("results$unique_summary        # Summary table\n")
cat("results$unique_targets        # List of unique genes per miRNA\n")
cat("results$focus_analysis        # Details for focus miRNA\n\n")

cat("# Re-run with different options:\n")
cat('results <- analyze_mirna_overlap(focus_mirna = "hsa-miR-21-5p")\n\n')

cat("================================================================================\n\n")
