#!/usr/bin/env Rscript

# Lightweight validation: no downloads and no analysis execution.

expected_files <- c(
  "BCC_Scripts/00_Setup_Data.R",
  "BCC_Scripts/01_Atlas_Integration.R",
  "BCC_Scripts/02_Mapping.R",
  "BCC_Scripts/03_CD4_Analysis.R",
  "BCC_Scripts/04_Alluvial_Plot.R",
  "BCC_Scripts/05_Original_Dataset_UMAP.R",
  "BCC_Scripts/06_Regression_Analysis_CD8.R",
  "NSCLC_Scripts/00_Setup_Data.R",
  "NSCLC_Scripts/Fig.1.R",
  "NSCLC_Scripts/Fig2.R",
  "NSCLC_Scripts/Fig3.R",
  "NSCLC_Scripts/Fig4.R",
  "NSCLC_Scripts/Fig.5.R",
  "NSCLC_Scripts/SFig1.R",
  "NSCLC_Scripts/SFig2.R"
)

missing_files <- expected_files[!file.exists(expected_files)]
if (length(missing_files)) {
  stop("Missing required scripts: ", paste(missing_files, collapse = ", "))
}

parse_failures <- character()
for (script in expected_files) {
  tryCatch(
    {
      parse(file = script)
      message("Parsed: ", script)
    },
    error = function(error) {
      parse_failures <<- c(parse_failures, paste0(script, ": ", conditionMessage(error)))
    }
  )
}
if (length(parse_failures)) {
  stop("R parse failures:\n", paste(parse_failures, collapse = "\n"))
}

required_packages <- c(
  "FNN", "Matrix", "RColorBrewer", "Seurat", "SeuratObject", "cowplot",
  "data.table", "dplyr", "ggalluvial", "ggplot2", "ggpubr", "ggrepel",
  "ggsci", "gridExtra", "pROC", "patchwork", "readr", "readxl", "scales",
  "stringr", "tidyr", "tidytext", "tidyverse", "viridis"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages)) {
  stop(
    "Missing required R packages: ", paste(missing_packages, collapse = ", "),
    "\nRun renv::restore() from the repository root."
  )
}

message("Validation passed: repository structure, R syntax, and packages are ready.")
