# Analysis Scripts

Scripts are organized in numbered directories reflecting execution order:

| Directory | Description | Key Outputs |
|-----------|-------------|-------------|
| `00_data_prep` | Download from GEO, filter probes | Filtered count matrix |
| `01_compositional` | Fractional abundance analysis | Figures 1-3 |
| `02_minor_miRNA` | Minor miRNA differential analysis | Figures 6, Table 4 |
| `03_deseq2` | Conventional DE analysis | Figure 7 |
| `04_pathway` | Target prediction & enrichment | Figures 4-5, Tables 1-3 |

## Running the Full Analysis

Option 1: Use Make
```bash
make all
```

Option 2: Run scripts in order
```r
source(here::here("analysis", "00_data_prep", "00_download_and_preprocess.R"))
rmarkdown::render(here::here("analysis", "01_compositional", "CompositionalAnalysis.Rmd"))
# etc.
```

## Dependencies

All R package dependencies are managed with `renv`. After cloning:
```r
renv::restore()
```
