# MSC-EV-miRNA-CF-BAL To-Do List
## Updated: January 3, 2026

---

## Current Status: Repository Ready for GitHub

Pipeline refactoring complete. All 11 figures and 5 tables reproduce from GEO data.

---

## Completed Tasks ✓

### Submission Archive (January 2, 2026)
- [x] Archive created and sent to Sara
- [x] Email sent to Dan & Claudia re: GEO sample mapping verification

### Pipeline Refactoring (January 3, 2026)
- [x] Renamed `analysis/04_pathway/run_all_analyses.R` → `run_factorial_analysis.R`
- [x] Created `analysis/04_pathway/create_pathway_figures.R` (orchestrates Figs 4, 5A-C, Tables 1-3)
- [x] Simplified root `run_all_analyses.R` to delegate to pathway orchestrator
- [x] Added grep tag for Fig 5B in `miRNA_pathway_workflow.R`
- [x] Updated `manuscript_manifest.csv` with `search_term` column
- [x] Verified full pipeline runs to completion
- [x] All 11 figures confirmed present
- [x] All tables confirmed present
- [x] 41 PDF versions generated

### Bug Fixes (Previous Sessions)
- [x] "Total Counts" row bug - filtered out spurious row in count matrix
- [x] GEO sample mapping - positions 27, 36, 46 are QC-excluded samples
- [x] dplyr::select namespace conflicts - explicit namespacing throughout
- [x] Bonferroni → BH correction - fixed multiple testing method
- [x] flextable::width collision - explicit namespacing in MinorMiRNA
- [x] Volcano plot labeling - restored full CF-relevant miRNA highlighting
- [x] top_n_pathways parameter - changed from 50 to NULL for full analysis

### Script Migration (Complete)
- [x] 00_download_and_preprocess.R
- [x] CompositionalAnalysis.Rmd
- [x] simpson_diversity_figure.R
- [x] MinorMiRNA_Prevalence_Analysis.Rmd
- [x] DESeq2_CF_vs_HC_analysis.R
- [x] miRNA_pathway_workflow.R
- [x] run_factorial_analysis.R (formerly run_all_analyses.R)
- [x] create_pathway_figures.R (new orchestrator)
- [x] analyze_mirna_target_overlap.R
- [x] plot_pathway_counts.R
- [x] categorize_pathways.R
- [x] analyze_pathway_similarity.R

---

## Pending Tasks

### Repository Finalization
- [ ] Create README.md with setup instructions and workflow overview
- [ ] Create .gitignore (exclude data/, results/, *.RData, etc.)
- [ ] Add LICENSE file
- [ ] Initialize git and push to GitHub

### Content Review Before Publication
- [ ] Review AI-curated content in Table 2 (Function column - inflammation genes)
- [ ] Review AI-curated content in Table 3 (Function column - remodeling genes)
- [ ] Review AI-curated content in Table 4 (Key_Mechanism column - CF-relevant miRNAs)

### Awaiting Response
- [ ] Dan & Claudia confirmation of QC-excluded samples (positions 27, 36, 46)

---

## Upcoming Deadlines

### Lab Meeting Presentation
- **Date:** Friday, January 24, 2026 at 1:00 PM
- **Start prep by:** January 17, 2026
- [ ] Create slide deck presenting MSC-EV-miRNA-CF-BAL work

---

## Author Comments to Address (Bruce)

### 1. Q-value vs BH FDR Explanation
Bruce asked about the difference between Storey's q-values and Benjamini-Hochberg FDR.

**Background:** We switched to Storey's q-values after discovering that BH-adjusted p-values can behave counterintuitively—the same raw p-value can yield a *smaller* adjusted p-value in a longer list than in a shorter one. Storey's method estimates π₀ (the proportion of true nulls) and is less susceptible to this issue.

**Action needed:** Add brief explanation to Methods or respond to reviewer comment.

**Reference:** Storey JD, Tibshirani R. Statistical significance for genomewide studies. PNAS 2003;100(16):9440-9445.

### 2. Systematic DESeq2 Analysis of All miRNAs
Bruce suggested a comprehensive differential expression analysis using DESeq2.

**Tom's assessment:** This is not well-suited for this paper because:
- The paper focuses on the most highly abundant miRNAs in hMSC EVs
- These major miRNAs are relatively stable across treatment conditions
- DESeq2 primarily identifies significant changes in low-abundance miRNAs
- Including this would shift focus away from the paper's central message

**Action needed:** Decide whether to add as supplementary analysis or respond explaining why it's not included.

---

## File Locations Quick Reference

### New/Renamed Files (January 3, 2026)
| File | Location |
|------|----------|
| run_factorial_analysis.R | analysis/04_pathway/ |
| create_pathway_figures.R | analysis/04_pathway/ |
| manuscript_manifest.csv | project root |

### Key Metadata Files
| File | Location |
|------|----------|
| sample_metadata.csv | data/metadata/ |
| geo_to_curated_sample_mapping.csv | data/metadata/ |
| CF_miRNA_annotations.csv | data/metadata/ |

---

## Reproduction Commands

```r
# Full pipeline (~10-15 minutes)
library(here)
setwd("~/Documents/MSC-EV-miRNA-CF-BAL")
source(here("run_all_analyses.R"))

# Create submission archive
source(here("archive_for_submission.R"))

# Find source code for any figure
# Use search_term from manuscript_manifest.csv:
# grep -r "volcano_CF_vs_HC" analysis/
```

---

## Lessons Learned

1. **Pipeline dependencies matter.** Upstream fixes require re-running all downstream analyses. Consider {targets} or Makefile for future projects.

2. **Grep-able filenames.** Dynamic filename construction (via `paste0`) breaks searchability. Add comment tags like `# MANUSCRIPT: Fig_5B` for figures with constructed names.

3. **Separate orchestration from computation.** Having `run_factorial_analysis.R` (does the heavy lifting) separate from `create_pathway_figures.R` (orchestrates manuscript outputs) makes the codebase easier to navigate.
