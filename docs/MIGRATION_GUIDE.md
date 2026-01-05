# Migration Guide: Adapting Scripts to New Repository Structure

This guide documents the changes needed when copying scripts from the original analysis to the new `MSC-EV-miRNA-CF-BAL` repository.

## Overview of Changes

| Original | New Structure |
|----------|---------------|
| Relative paths (`../raw_count.xlsx`) | `here::here("data", "raw", "raw_count.xlsx")` |
| Outputs in same directory as script | Outputs to `results/{module}/` |
| Scripts in results folders | Scripts only in `analysis/{module}/` |
| Ad-hoc config in scripts | Centralized `config/analysis_config.yaml` |

---

## Step-by-Step Migration

### 1. Add `here` to script headers

**Before:**
```r
library(readxl)
library(dplyr)
```

**After:**
```r
library(here)
library(readxl)
library(dplyr)
```

### 2. Update input file paths

**Before:**
```r
raw_data <- read_excel("../raw_count.xlsx")
# or
raw_data <- read_excel("../../data/raw_count.xlsx")
```

**After:**
```r
raw_data <- read_excel(here("data", "raw", "raw_count.xlsx"))
```

### 3. Update output directories

**Before:**
```r
ggsave("major_miRNAs_barplot.png", p, width = 10, height = 8, dpi = 300)
write.csv(results, "analysis_results.csv", row.names = FALSE)
```

**After:**
```r
# Define output directory at top of script
output_dir <- here("results", "01_compositional")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Use file.path for outputs
ggsave(file.path(output_dir, "major_miRNAs_barplot.png"), p, width = 10, height = 8, dpi = 300)
write.csv(results, file.path(output_dir, "analysis_results.csv"), row.names = FALSE)
```

### 4. Update config loading

**Before:**
```r
source("analysis_config.R")
# or inline definitions:
q_cutoff <- 0.05
min_databases <- 2
```

**After:**
```r
library(yaml)
config <- read_yaml(here("config", "analysis_config.yaml"))
q_cutoff <- config$stringency$stringent$q_cutoff
min_databases <- config$stringency$stringent$min_databases
```

### 5. Update cross-references between scripts

**Before:**
```r
# Load results from another script
prev_results <- read.csv("../miRNA_compositional_analysis/major_miRNAs_1pct_summary.csv")
```

**After:**
```r
# Load results from previous analysis step
prev_results <- read.csv(here("results", "01_compositional", "major_miRNAs_1pct_summary.csv"))
```

### 6. Update sample metadata loading

**Before:**
```r
raw_sheets <- lapply(raw_names, function(x) read_excel("../raw_count.xlsx", sheet = x))
sample_mapping <- raw_sheets$Sample_Group_Mapping
```

**After:**
```r
# Load sample metadata (curated file with treatment groups)
sample_metadata <- read_csv(here("data", "metadata", "sample_metadata.csv"))

# Filter to analysis samples
samples_to_analyze <- sample_metadata %>%
  filter(In_Current_Study == TRUE) %>%
  pull(Sample)
```

---

## File-by-File Migration Checklist

### 00_data_prep

| Original File | New Location | Status |
|---------------|--------------|--------|
| (download from GEO) | `analysis/00_data_prep/00_download_and_preprocess.R` | Created |

### 01_compositional

| Original File | New Location | Changes Needed |
|---------------|--------------|----------------|
| `CompositionalAnalysis_1pct_with_diversity.Rmd` | `analysis/01_compositional/CompositionalAnalysis.Rmd` | Update paths, output_dir |

Key path changes:
```r
# Input
read_excel(here("data", "raw", "GSE282919_...xlsx"), sheet = "Raw")
sample_metadata <- read_csv(here("data", "metadata", "sample_metadata.csv"))

# Output
output_dir <- here("results", "01_compositional")
```

### 02_minor_miRNA

| Original File | New Location | Changes Needed |
|---------------|--------------|----------------|
| `MinorMiRNA_Prevalence_Analysis.Rmd` | `analysis/02_minor_miRNA/MinorMiRNA_Prevalence_Analysis.Rmd` | Update paths |

Key path changes:
```r
# Input (depends on 01_compositional results)
major_mirnas <- read_csv(here("results", "01_compositional", "major_miRNAs_1pct_summary.csv"))

# Output
output_dir <- here("results", "02_minor_miRNA")
```

### 03_deseq2

| Original File | New Location | Changes Needed |
|---------------|--------------|----------------|
| `DESeq2_CF_vs_HC_analysis.R` | `analysis/03_deseq2/DESeq2_CF_vs_HC_analysis.R` | Update paths |

Key path changes:
```r
# Input
count_data <- read_csv(here("data", "processed", "mirna_counts_filtered.csv"))
sample_metadata <- read_csv(here("data", "metadata", "sample_metadata.csv"))

# Output
output_dir <- here("results", "03_deseq2")
```

### 04_pathway

| Original File | New Location | Changes Needed |
|---------------|--------------|----------------|
| `analysis_config.R` | `config/analysis_config.yaml` | Convert to YAML (done) |
| `miRNA_pathway_workflow.R` | `analysis/04_pathway/miRNA_pathway_workflow.R` | Update paths |
| `analyze_mirna_target_overlap.R` | `analysis/04_pathway/analyze_mirna_target_overlap.R` | Update paths |
| `analyze_pathway_similarity.R` | `analysis/04_pathway/analyze_pathway_similarity.R` | Update paths |
| `categorize_pathways.R` | `analysis/04_pathway/categorize_pathways.R` | Update paths |
| `plot_pathway_counts.R` | `analysis/04_pathway/plot_pathway_counts.R` | Update paths |
| `run_analysis.Rmd` | `analysis/04_pathway/run_analysis.Rmd` | Update paths |
| `summarize_evidence.R` | `analysis/04_pathway/summarize_evidence.R` | Update paths |

Key path changes:
```r
# Config
config <- read_yaml(here("config", "analysis_config.yaml"))

# Input (depends on 01_compositional)
major_mirnas <- read_csv(here("results", "01_compositional", "major_miRNAs_1pct_summary.csv"))

# Output
output_dir <- here("results", "04_pathway")
```

---

## Common Patterns

### Pattern 1: Script header template

```r
# ==============================================================================
# [Script Name]
# ==============================================================================
# [Description]
#
# INPUTS:
#   - data/metadata/sample_metadata.csv
#   - results/01_compositional/[file].csv
#
# OUTPUTS:
#   - results/[module]/[outputs]
# ==============================================================================

# Load packages
library(here)
library(tidyverse)
# ... other packages

# Configuration
config <- yaml::read_yaml(here("config", "analysis_config.yaml"))
output_dir <- here("results", "[module]")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load data
sample_metadata <- read_csv(here("data", "metadata", "sample_metadata.csv"))
# ... rest of script
```

### Pattern 2: Rmd YAML header

```yaml
---
title: "[Analysis Title]"
author: "DartCF Bioinformatics Core"
date: "`r Sys.Date()`"
output:
  html_document:
    toc: true
    toc_float: true
    code_folding: hide
params:
  output_dir: !r here::here("results", "[module]")
---
```

### Pattern 3: Saving multiple output types

```r
# Save both PNG and PDF
save_figure <- function(plot, filename, width = 10, height = 8) {
  ggsave(file.path(output_dir, paste0(filename, ".png")), plot, 
         width = width, height = height, dpi = 300)
  ggsave(file.path(output_dir, paste0(filename, ".pdf")), plot, 
         width = width, height = height)
  cat("Saved:", filename, ".png/.pdf\n")
}

# Usage
save_figure(p, "major_miRNAs_barplot")
```

---

## Testing Migration

After migrating each script:

1. **Clear environment**: `rm(list = ls())`
2. **Set working directory to repo root**: `setwd("/path/to/MSC-EV-miRNA-CF-BAL")`
3. **Run script**: Verify it finds all inputs and writes to correct outputs
4. **Check outputs**: Compare to original analysis results

---

## Troubleshooting

### "File not found" errors

- Verify `here()` is returning the repository root: `here::here()`
- Check that `.here` file or `.Rproj` file exists in repo root
- Ensure previous analysis steps have run (check `results/{module}/` for inputs)

### Package version issues

- Run `renv::restore()` to install locked package versions
- If a package fails, check `renv.lock` for the expected version

### Output not appearing in expected location

- Verify `output_dir` is set correctly at top of script
- Check that `dir.create()` ran without errors
- Look for outputs in working directory (may indicate `here()` not working)
