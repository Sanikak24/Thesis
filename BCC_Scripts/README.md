# BCC Analysis Pipeline

This directory contains the scripts used to analyze the basal cell carcinoma (BCC) single-cell RNA sequencing dataset from Yost et al. (2019).

## Dataset

- **Study:** Yost et al., Nature Medicine (2019)
- **GEO:** GSE123813
- **Reference Atlas:** MD Anderson CD8 T-cell Atlas

## Pipeline

### 01_Atlas_Integration.R
Maps BCC CD8 T cells onto the reference CD8 T-cell atlas using Seurat label transfer.

### 02_Mapping.R
Generates predicted atlas labels and performs majority-vote mapping.

### 03_CD4_Analysis.R
Maps CD4 T cells onto the CD4 reference atlas and generates visualization outputs.

### 04_Alluvial_Plot.R
Creates alluvial plots comparing original annotations with predicted atlas states.

### 05_Original_Dataset_UMAP.R
Generates UMAP visualizations of the original BCC dataset.

### 06_Regression_Analysis_CD8.R
Builds logistic regression models to evaluate associations between predicted CD8 states and immunotherapy response.

## Running the pipeline

Run all scripts from the project root (`sigNATURE-Framework_backup`):

```bash
Rscript BCC_Scripts/01_Atlas_Integration.R
Rscript BCC_Scripts/02_Mapping.R
Rscript BCC_Scripts/03_CD4_Analysis.R
Rscript BCC_Scripts/04_Alluvial_Plot.R
Rscript BCC_Scripts/05_Original_Dataset_UMAP.R
Rscript BCC_Scripts/06_Regression_Analysis_CD8.R
```

Script 06 requires the majority-vote mapped object written by script 04.
Outputs are written under `results/bcc/`.
