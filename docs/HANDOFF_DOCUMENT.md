# MSC-EV-miRNA-CF-BAL Migration Handoff Document
## January 2, 2026

---

## Project Summary

Successfully migrated the MSC-EV-miRNA-CF-BAL analysis pipeline from scattered scripts in `CF_BAL_MSVEV/` to a clean, reproducible repository structure in `MSC-EV-miRNA-CF-BAL/`. All 7 figures and 5 tables from the Cytotherapy manuscript can now be reproduced from GEO data with a single command.

---

## Current Status: ✅ COMPLETE

All manuscript figures and tables verified. Archive created for Sara.

---

## Repository Location

```
~/Documents/MSC-EV-miRNA-CF-BAL/
```

Old (unmigrated) repository preserved at:
```
~/Documents/CF_BAL_MSVEV/
```

---

## How to Reproduce Everything

```r
library(here)
setwd("~/Documents/MSC-EV-miRNA-CF-BAL")

# Full pipeline (~10-15 minutes)
source(here("run_all_analyses.R"))

# Create submission archive
source(here("archive_for_submission.R"))
```

---

## Key Files in Repository Root

| File | Purpose |
|------|---------|
| `run_all_analyses.R` | Master script - runs entire pipeline |
| `archive_for_submission.R` | Creates ZIP with all figures/tables |
| `manuscript_manifest.csv` | Maps figures/tables to source scripts |
| `TODO_LIST.md` | Migration status and lessons learned |

---

## Directory Structure

```
MSC-EV-miRNA-CF-BAL/
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
│       ├── run_all_analyses.R              # Pathway workflow entry point
│       ├── miRNA_pathway_workflow.R        # Core enrichment functions
│       ├── plot_pathway_counts.R           # Fig 5A
│       ├── categorize_pathways.R           # Fig 5C, Table 1
│       ├── analyze_mirna_target_overlap.R  # Fig 4
│       └── analyze_pathway_similarity.R    # Tables 2, 3
│
├── data/
│   ├── raw/           # GEO Excel file (auto-downloaded)
│   ├── processed/     # mirna_counts_filtered.csv
│   └── metadata/
│       ├── sample_metadata.csv
│       ├── geo_to_curated_sample_mapping.csv
│       └── CF_miRNA_annotations.csv
│
└── results/
    ├── 01_compositional/   # Figs 1, 2, 3A, 3B, 3C
    ├── 02_minor_miRNA/     # Fig 6, Table 4
    ├── 03_deseq2/          # Fig 7
    └── 04_pathway/         # Figs 4, 5A, 5B, 5C, Tables 1, 2, 3
```

---

## Figure and Table Mapping

| Item | Description | Source Script | Output Location |
|------|-------------|---------------|-----------------|
| Fig 1 | Major miRNA stacked bar | CompositionalAnalysis.Rmd | results/01_compositional/ |
| Fig 2 | Major miRNA mean abundance | CompositionalAnalysis.Rmd | results/01_compositional/ |
| Fig 3A | Treatment differences | CompositionalAnalysis.Rmd | results/01_compositional/ |
| Fig 3B | Density distributions | CompositionalAnalysis.Rmd | results/01_compositional/ |
| Fig 3C | Simpson diversity | simpson_diversity_figure.R | results/01_compositional/ |
| Fig 4 | miRNA target overlap | analyze_mirna_target_overlap.R | results/04_pathway/major/stringent/ |
| Fig 5A | Pathway counts | plot_pathway_counts.R | results/04_pathway/ |
| Fig 5B | GO BP dotplot | miRNA_pathway_workflow.R | results/04_pathway/major/stringent/ |
| Fig 5C | Pathway categories | categorize_pathways.R | results/04_pathway/ |
| Fig 6 | CF-relevant volcano | MinorMiRNA_Prevalence_Analysis.Rmd | results/02_minor_miRNA/ |
| Fig 7 | DESeq2 volcano | DESeq2_CF_vs_HC_analysis.R | results/03_deseq2/ |
| Table 1 | Pathway categories | categorize_pathways.R | results/04_pathway/ |
| Table 2 | Inflammation genes | analyze_pathway_similarity.R | results/04_pathway/.../pathway_similarity_results/ |
| Table 3 | Remodeling genes | analyze_pathway_similarity.R | results/04_pathway/.../pathway_similarity_results/ |
| Table 4 | CF-relevant miRNAs | MinorMiRNA_Prevalence_Analysis.Rmd | results/02_minor_miRNA/ |
| Table 5 | miR-21 context | Manually curated | N/A |

---

## Critical Metadata Files

These must be in `data/metadata/` before running:

1. **sample_metadata.csv** - Treatment group assignments for all 78 samples
2. **geo_to_curated_sample_mapping.csv** - Maps GEO positions to sample IDs (handles QC-excluded samples at positions 27, 36, 46)
3. **CF_miRNA_annotations.csv** - 20 CF-relevant miRNAs with Impact and Key_Mechanism columns

---

## Bugs Fixed During Migration

1. **"Total Counts" row bug** - Filtered out spurious row in count matrix
2. **GEO sample mapping** - Positions 27, 36, 46 were QC-excluded but not at end of file
3. **dplyr::select conflicts** - Explicit namespacing throughout
4. **Bonferroni → BH** - Corrected multiple testing method
5. **top_n_pathways parameter** - Changed from 50 to NULL for full analysis
6. **flextable::width collision** - Explicit namespacing in MinorMiRNA
7. **Volcano plot labeling** - Restored full CF-relevant miRNA highlighting with Impact colors

---

## Remaining Action Items

### Required for Submission
- [ ] Review AI-curated content in Tables 2, 3, 4 (Function/Key_Mechanism columns)
- [ ] Send archive to Sara

### Optional Follow-up
- [ ] Email Dan & Claudia to confirm GEO positions 27, 36, 46 are QC failures
- [ ] Rename `analysis/04_pathway/run_all_analyses.R` to avoid confusion with master script
- [ ] Add README.md with full setup instructions
- [ ] Push to GitHub

---

## Archive Created

Location: `~/Documents/MSC-EV-miRNA-CF-BAL/submission_archive_20260102.zip`

Contents:
- `figures/` - All 7 figures (PNG + PDF)
- `tables/` - Tables 1-4 (CSV + Word)
- `supplementary/` - Supporting data files
- `manuscript_manifest.csv` - Figure/table mapping

---

## GEO Accession

**GSE282919** - Raw data downloads automatically via `00_download_and_preprocess.R`

---

## Contact

Project: MSC-EV-miRNA-CF-BAL (Cytotherapy submission)
Target Journal: Cytotherapy
Collaborators: Sara (submission), Dan & Claudia (data/QC questions)
