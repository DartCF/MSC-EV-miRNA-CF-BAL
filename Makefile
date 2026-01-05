# ==============================================================================
# Makefile for MSC-EV-miRNA-CF-BAL Analysis Pipeline
# ==============================================================================

.PHONY: all clean data compositional minor deseq2 pathway figures help

# Default target
all: data compositional minor deseq2 pathway figures

# Help
help:
	@echo "MSC-EV-miRNA-CF-BAL Analysis Pipeline"
	@echo ""
	@echo "Targets:"
	@echo "  all           Run full analysis pipeline"
	@echo "  data          Download and preprocess data from GEO"
	@echo "  compositional Run compositional analysis (Figures 1-3)"
	@echo "  minor         Run minor miRNA analysis (Figure 6, Table 4)"
	@echo "  deseq2        Run DESeq2 analysis (Figure 7)"
	@echo "  pathway       Run pathway analysis (Figures 4-5, Tables 1-3)"
	@echo "  figures       Copy final figures to figures/"
	@echo "  clean         Remove all generated files"
	@echo ""
	@echo "Prerequisites:"
	@echo "  1. R with renv installed"
	@echo "  2. Run 'make restore' to install R packages"

# Restore R environment
restore:
	Rscript -e "renv::restore()"

# Data preprocessing
data: results/00_data_prep/.done

results/00_data_prep/.done: analysis/00_data_prep/00_download_and_preprocess.R
	Rscript -e "source('analysis/00_data_prep/00_download_and_preprocess.R')"
	touch $@

# Compositional analysis
compositional: results/01_compositional/.done

results/01_compositional/.done: results/00_data_prep/.done analysis/01_compositional/*.Rmd
	Rscript -e "rmarkdown::render('analysis/01_compositional/CompositionalAnalysis.Rmd', output_dir = 'results/01_compositional')"
	touch $@

# Minor miRNA analysis
minor: results/02_minor_miRNA/.done

results/02_minor_miRNA/.done: results/01_compositional/.done analysis/02_minor_miRNA/*.Rmd
	Rscript -e "rmarkdown::render('analysis/02_minor_miRNA/MinorMiRNA_Prevalence_Analysis.Rmd', output_dir = 'results/02_minor_miRNA')"
	touch $@

# DESeq2 analysis
deseq2: results/03_deseq2/.done

results/03_deseq2/.done: results/00_data_prep/.done analysis/03_deseq2/*.R
	Rscript analysis/03_deseq2/DESeq2_CF_vs_HC_analysis.R
	touch $@

# Pathway analysis
pathway: results/04_pathway/.done

results/04_pathway/.done: results/01_compositional/.done analysis/04_pathway/*.R
	cd analysis/04_pathway && Rscript run_analysis.R
	touch $@

# Copy final figures
figures: compositional minor deseq2 pathway
	@echo "Copy final figures from results/ to figures/"
	@echo "(Manual step - review outputs first)"

# Clean generated files
clean:
	rm -f results/*/.done
	rm -f results/**/*.csv
	rm -f results/**/*.png
	rm -f results/**/*.pdf
	rm -f data/processed/*.csv
