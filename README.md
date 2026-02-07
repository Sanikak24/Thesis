# Thesis analysis repository

This repository contains R scripts for my thesis analyses.

## Contents
- **BCC dataset**: `bcc/Scripts/`
  - `00_initial_pipeline.R`
  - `01_alluvial_plot.R`
- **NSCLC dataset**: `Thesis/Scripts/`
  - `01_fig1.R` … `07_sfig2.R`

## Suggested run order
### BCC
```bash
Rscript bcc/scripts/00_initial_pipeline.R
Rscript bcc/scripts/01_alluvial_plot.R
```

### NSCLC
```bash
Rscript nsclc/scripts/01_fig1.R
Rscript nsclc/scripts/02_fig2.R
Rscript nsclc/scripts/03_fig3.R
Rscript nsclc/scripts/04_fig4.R
Rscript nsclc/scripts/05_fig5.R
Rscript nsclc/scripts/06_sfig1.R
Rscript nsclc/scripts/07_sfig2.R
```

## Outputs
Write generated outputs to:
- `results/bcc/` for BCC results
- `results/nsclc/` for NSCLC results

## Data
Raw data are **not** included. Put inputs under `data/` (see `data/README.md`).
