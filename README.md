# MSC-EV-miRNA-CF-BAL

Analysis pipeline for the manuscript: **[Title TBD - Cytotherapy submission]**

This repository contains all code to reproduce the figures and tables from publicly available data.

## GEO Accession

**GSE282919** - Raw miRNA count data is downloaded automatically by the pipeline.

## Quick Start

```r
# Clone the repository and open in RStudio, then:
library(here)
source(here("run_all_analyses.R"))  # ~10-15 minutes

# Create submission archive with all figures/tables
source(here("archive_for_submission.R"))
```

## Prerequisites

### R Version
R ≥ 4.0 recommended

### Package Installation

This project uses `renv` for reproducible package management. After cloning:

```r
install.packages("renv")
renv::restore()
```

This installs all 168 required packages (including Bioconductor dependencies) at the exact versions used in the analysis.

**Key packages:** DESeq2, clusterProfiler, ReactomePA, multiMiR, enrichplot, qvalue, tidyverse, flextable

## Directory Structure

```
MSC-EV-miRNA-CF-BAL/
├── run_all_analyses.R          # Master script - runs entire pipeline
├── archive_for_submission.R    # Creates ZIP with figures/tables
├── manuscript_manifest.csv     # Maps figures/tables to source scripts
│
├── analysis/
│   ├── 00_data_prep/
│   │   └── 00_download_and_preprocess.R    # Downloads from GEO
│   ├── 01_compositional/
│   │   ├── CompositionalAnalysis.Rmd       # Figs 1, 2, 3A, 3B
│   │   └── simpson_diversity_figure.R      # Fig 3C
│   ├── 02_minor_miRNA/
│   │   └── MinorMiRNA_Prevalence_Analysis.Rmd  # Fig 6, Table 4
│   ├── 03_deseq2/
│   │   └── DESeq2_CF_vs_HC_analysis.R      # Fig 7
│   └── 04_pathway/
│       ├── create_pathway_figures.R        # Orchestrates Figs 4, 5A-C
│       ├── run_factorial_analysis.R        # 2×2 factorial pathway analysis
│       ├── miRNA_pathway_workflow.R        # Core enrichment functions
│       ├── plot_pathway_counts.R           # Fig 5A
│       ├── categorize_pathways.R           # Fig 5C, Table 1
│       ├── analyze_mirna_target_overlap.R  # Fig 4
│       └── analyze_pathway_similarity.R    # Tables 2, 3
│
├── data/
│   ├── raw/                    # GEO data (auto-downloaded, not in repo)
│   ├── processed/              # Cleaned count matrix
│   └── metadata/               # Sample annotations (included in repo)
│       ├── sample_metadata.csv
│       ├── geo_to_curated_sample_mapping.csv
│       └── CF_miRNA_annotations.csv
│
└── results/                    # Generated outputs (not in repo)
    ├── 01_compositional/       # Figs 1, 2, 3A, 3B, 3C
    ├── 02_minor_miRNA/         # Fig 6, Table 4
    ├── 03_deseq2/              # Fig 7
    └── 04_pathway/             # Figs 4, 5A, 5B, 5C, Tables 1, 2, 3
```

## Manuscript Outputs

| Figure | Description | Source Script |
|--------|-------------|---------------|
| Fig 1 | Major miRNA stacked bar by sample | CompositionalAnalysis.Rmd |
| Fig 2 | Major miRNA mean abundance | CompositionalAnalysis.Rmd |
| Fig 3A | Treatment differences | CompositionalAnalysis.Rmd |
| Fig 3B | Density distributions | CompositionalAnalysis.Rmd |
| Fig 3C | Simpson diversity | simpson_diversity_figure.R |
| Fig 4 | miRNA target overlap | analyze_mirna_target_overlap.R |
| Fig 5A | Pathway counts | plot_pathway_counts.R |
| Fig 5B | GO BP dotplot | miRNA_pathway_workflow.R |
| Fig 5C | Pathway categories | categorize_pathways.R |
| Fig 6 | CF-relevant volcano | MinorMiRNA_Prevalence_Analysis.Rmd |
| Fig 7 | DESeq2 volcano | DESeq2_CF_vs_HC_analysis.R |

| Table | Description | Source Script |
|-------|-------------|---------------|
| Table 1 | Pathway categories | categorize_pathways.R |
| Table 2 | Inflammation genes | analyze_pathway_similarity.R |
| Table 3 | Remodeling genes | analyze_pathway_similarity.R |
| Table 4 | CF-relevant miRNAs | MinorMiRNA_Prevalence_Analysis.Rmd |
| Table 5 | miR-21 context | Manually curated |

For grep-able search terms, see `manuscript_manifest.csv`.

## Finding Source Code

To find the code that generates a specific figure:

```bash
# Use the search_term from manuscript_manifest.csv
grep -r "volcano_CF_vs_HC" analysis/
```

## Statistical Methods

### FDR Correction
Pathway enrichment uses Storey's q-value method rather than Benjamini-Hochberg. This avoids counterintuitive behavior where the same raw p-value can yield different adjusted p-values depending on list length.

Reference: Storey JD, Tibshirani R. Statistical significance for genomewide studies. PNAS 2003;100(16):9440-9445.

### Compositional Analysis
miRNA abundances are treated as fractional (compositional) data. Major miRNAs are defined as those with ≥1% mean abundance across samples.

## License

[TBD]

## Citation

[TBD - add citation once published]

## Contact

[TBD]
# MSC-EV-miRNA-CF-BAL
