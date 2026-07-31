library(Seurat)
library(ggplot2)
library(scales)
library(patchwork)
library(dplyr)
library(tidyr)

nsclc_data_dir <- file.path("data", "nsclc")
reference_dir <- file.path("data", "reference")
results_dir <- file.path("results", "nsclc")
intermediate_dir <- file.path(results_dir, "intermediate")

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(intermediate_dir, recursive = TRUE, showWarnings = FALSE)

require_input <- function(path) {
  if (
    !file.exists(path) ||
    is.na(file.info(path)$size) ||
    file.info(path)$size == 0
  ) {
    stop(
      "Required input file is missing or empty: ",
      path,
      "\nRun: Rscript NSCLC_Scripts/00_Setup_Data.R"
    )
  }
}

cd8_reference_file <- file.path(reference_dir, "CD8_Obj_for_mapping.rds")
require_input(cd8_reference_file)
CD8_Obj <- readRDS(cd8_reference_file)
genes_ordered <- c(
  # Immediate-early activation markers
  "FOS", "FOSB", "CD69",
  
  # Cytotoxic effectors
  "GZMA", "GZMK",
  
  # Chemokines
  "CCL4", "CCL4L2", "CXCL13",
  
  # Chemokine receptor
  "CXCR6",
  "EOMES", "TXNIP", "FKBP5")
gene_lists <- lapply(genes_ordered, function(g) g)
names(gene_lists) <- genes_ordered

# 5️⃣ Add module scores to CD8 object
CD8_Obj <- AddModuleScore(
  object   = CD8_Obj,
  features = gene_lists,
  name     = "Gene"
)

# 6️⃣ Define color scale
low_high_cols <- c("grey80", "red")

# 7️⃣ Create FeaturePlots for each gene
plot_list <- list()
for (i in seq_along(genes_ordered)) {
  feature_name <- paste0("Gene", i)
  
  p <- FeaturePlot(
    CD8_Obj,
    features   = feature_name,
    cols       = low_high_cols,
    min.cutoff = 0,
    max.cutoff = 2,
    raster     = TRUE
  ) +
    ggtitle(genes_ordered[i]) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 10)
    )
  
  plot_list[[i]] <- p
}

names(plot_list) <- genes_ordered

# 8️⃣ Combine all FeaturePlots into a grid (5 columns)
combined <- wrap_plots(plot_list, ncol = 5)
ggsave(
  file.path(results_dir, "CD8_gene_expression_rearranged.svg"),
  plot = combined,
  width = 16,
  height = 12,
  device = "svg"
)
