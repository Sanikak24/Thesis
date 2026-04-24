# Thesis analysis repository

This repository contains R scripts for my thesis analyses.

# Datasets
1) Reference atlas: https://singlecell.mdanderson.org/TCM/
2) (NSCLC) Liu et al: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE179994
3) (BCC) Yost et al: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE123813 
4) Clinical meta data taken from the Supplementary tables of the respective papers:
1. Liu, B., Hu, X., Feng, K. et al. Temporal single-cell tracing reveals clonal revival and expansion of precursor exhausted T cells during anti-PD-1 therapy in lung cancer. Nat Cancer 3, 108–121 (2022). https://doi.org/10.1038/s43018-021-00292-8
2. Yost, K.E., Satpathy, A.T., Wells, D.K. et al. Clonal replacement of tumor-specific T cells following PD-1 blockade. Nat Med 25, 1251–1259 (2019). https://doi.org/10.1038/s41591-019-0522-3

## Contents
- **BCC dataset**: `bcc/Scripts/`
  - `Mapping.R`
  - `Alluvial_plot.R`
- **NSCLC dataset**: `Thesis/Scripts/`
  - `01_fig1.R` … `07_sfig2.R`

## Suggested run order
### BCC
```bash
Rscript bcc/scripts/Mapping.R
Rscript bcc/scripts/Alluvial_plot.R
```

### NSCLC
```bash
Rscript nsclc/scripts/fig1.R
Rscript nsclc/scripts/fig2.R
Rscript nsclc/scripts/fig3.R
Rscript nsclc/scripts/fig4.R
Rscript nsclc/scripts/fig5.R
Rscript nsclc/scripts/Sfig1.R
Rscript nsclc/scripts/Sfig2.R
```

## Outputs
Write generated outputs to:
- `results/bcc/` for BCC results
- `results/nsclc/` for NSCLC results

## Data
Raw data are **not** included. Put inputs under `data/` (see `data/README.md`).
