# sigNATURE Framework

This repository contains R workflows for mapping T cells from basal cell
carcinoma (BCC) and non-small-cell lung cancer (NSCLC) datasets to shared
CD4/CD8 reference atlases.

Run every command below from the repository root.

## Datasets

- Reference T-cell atlases: <https://singlecell.mdanderson.org/TCM/>
- Liu et al. NSCLC dataset (GSE179994):
  <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE179994>
- Yost et al. BCC dataset (GSE123813):
  <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE123813>

Clinical annotations are obtained from the supplementary tables of:

- Liu, B., Hu, X., Feng, K. et al. *Temporal single-cell tracing reveals
  clonal revival and expansion of precursor exhausted T cells during
  anti-PD-1 therapy in lung cancer.* Nature Cancer 3, 108–121 (2022).
  <https://doi.org/10.1038/s43018-021-00292-8>
- Yost, K.E., Satpathy, A.T., Wells, D.K. et al. *Clonal replacement of
  tumor-specific T cells following PD-1 blockade.* Nature Medicine 25,
  1251–1259 (2019).
  <https://doi.org/10.1038/s41591-019-0522-3>

Raw data and generated results are intentionally excluded from Git.

## Repository layout

```text
BCC_Scripts/                  BCC setup and analysis scripts
NSCLC_Scripts/                NSCLC setup and analysis scripts
data/reference/               Shared CD4/CD8 reference atlases
data/bcc/                     BCC inputs
data/nsclc/                   NSCLC inputs
results/bcc/                  BCC outputs
results/nsclc/                NSCLC outputs
results/nsclc/intermediate/   Mapped NSCLC Seurat objects and caches
```

## NSCLC workflow

### 1. Prepare data

```bash
Rscript NSCLC_Scripts/00_Setup_Data.R
```

The setup script downloads or validates the Liu count matrix and metadata,
extracts the clinical response table from the official Liu supplementary
workbook, and validates the shared reference objects:

```text
data/reference/CD8_Obj_for_mapping.rds
data/reference/CD4_Obj_for_mapping.rds
data/nsclc/GSE179994_all.Tcell.rawCounts.rds
data/nsclc/GSE179994_Tcell.metadata.tsv
data/nsclc/response_info.xlsx
```

### 2. Run analyses

```bash
Rscript NSCLC_Scripts/Fig.1.R
Rscript NSCLC_Scripts/Fig2.R
Rscript NSCLC_Scripts/Fig3.R
Rscript NSCLC_Scripts/Fig4.R
Rscript NSCLC_Scripts/Fig.5.R
Rscript NSCLC_Scripts/SFig1.R
Rscript NSCLC_Scripts/SFig2.R
```

The dependency order is:

- `Fig2.R` creates `results/nsclc/intermediate/liu_cd8_mapped.rds`.
- `Fig3.R` uses the mapped CD8 object.
- `Fig4.R` creates `results/nsclc/intermediate/liu_cd4_mapped.rds`.
- `Fig.5.R` uses both mapped objects.
- `SFig1.R` uses the CD8 reference atlas.
- `SFig2.R` uses the mapped CD8 object from `Fig2.R`.

The mapping scripts can require substantial memory because the NSCLC count
matrix is approximately 2.3 GB.

## BCC workflow

Prepare the data and run the numbered scripts in order:

```bash
Rscript BCC_Scripts/00_Setup_Data.R
Rscript BCC_Scripts/01_Atlas_Integration.R
Rscript BCC_Scripts/02_Mapping.R
Rscript BCC_Scripts/03_CD4_Analysis.R
Rscript BCC_Scripts/04_Alluvial_Plot.R
Rscript BCC_Scripts/05_Original_Dataset_UMAP.R
Rscript BCC_Scripts/06_Regression_Analysis_CD8.R
```

See `BCC_Scripts/README.md` for BCC-specific details.

## Outputs

- BCC results: `results/bcc/`
- NSCLC results: `results/nsclc/`
- NSCLC intermediate objects: `results/nsclc/intermediate/`

## Validation

Parse-check all NSCLC scripts with:

```bash
for f in NSCLC_Scripts/*.R; do
  echo "Checking $f"
  Rscript -e "parse(file='$f')" || exit 1
done
```
