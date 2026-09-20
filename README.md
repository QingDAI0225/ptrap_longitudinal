[![DOI](https://zenodo.org/badge/1378424523.svg)](https://doi.org/10.5281/zenodo.22859443)

# P-trap longitudinal microbiome analysis

This repository contains the finalized R analysis workflow associated with the article available at [https://doi.org/10.64898/2026.05.13.724980](https://doi.org/10.64898/2026.05.13.724980), covering longitudinal analysis of full-length 16S sequencing data from p-trap samples, source-support modeling, final-day recovery analysis, and metagenomic validation.

## Workflow

Run the analysis in this order:

1. `01-PacBio_data_analysis.Rmd` — preprocessing, contaminant filtering, diversity analyses, PCA, and longitudinal PCA/GAM analyses.
2. `02-Source_support_compute.Rmd` — source-support and taxon-attribution computations.
3. `03-data_visualization.Rmd` — figures and summaries generated from the source-support results.
4. `final_day_recovery_analysis.Rmd` — projection and recovery analysis for the final sampling day.
5. `metagenomic_validation.Rmd` — comparison of longitudinal 16S and metagenomic patterns.

`longitudinal_pca_functions.R` contains the shared longitudinal PCA and GAM functions used by the first analysis.

## Software

The workflow was reproduced with R 4.3.2. Major dependencies include:

- CRAN: `dplyr`, `tidyr`, `ggplot2`, `vroom`, `stringr`, `tibble`, `purrr`, `readr`, `ggrepel`, `patchwork`, `vegan`, `mgcv`, and `broom`
- Bioconductor: `decontam` and `phyloseq`

## Input data and paths

Input data are not included in this repository. Paths can be configured with environment variables:

- `PTRAP_DATA_DIR` — directory containing the 16S tables and metadata
- `PTRAP_METAGENOME_FILE` — metagenomic relative-abundance table
- `PTRAP_CACHE_DIR` — directory for intermediate RDS files
- `PTRAP_OUTPUT_DIR` — directory for generated tables and figures

For example:

```bash
export PTRAP_DATA_DIR=/path/to/HiFi_16S
export PTRAP_METAGENOME_FILE=/path/to/Mondays-S-Relative-Abundance-Counts.txt
export PTRAP_CACHE_DIR="$PWD/cache"
export PTRAP_OUTPUT_DIR="$PWD/results"
mkdir -p "$PTRAP_CACHE_DIR" "$PTRAP_OUTPUT_DIR"
```

The 16S data directory is expected to contain:

- `vsearch_merged_freq_tax_691.tsv`
- `vsearch_merged_freq_tax_692.tsv`
- `metadata.csv`
- `DNA_concentration.csv`
- `volume_metadata.csv`

## Rendering

Render the R Markdown files sequentially from the repository root. For example:

```bash
Rscript -e 'rmarkdown::render("01-PacBio_data_analysis.Rmd")'
Rscript -e 'rmarkdown::render("02-Source_support_compute.Rmd")'
Rscript -e 'rmarkdown::render("03-data_visualization.Rmd")'
Rscript -e 'rmarkdown::render("final_day_recovery_analysis.Rmd")'
Rscript -e 'rmarkdown::render("metagenomic_validation.Rmd")'
```
