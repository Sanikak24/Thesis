library(Seurat)
library(ggplot2)
library(scales)
library(patchwork)
library(dplyr)
library(tidyr)
CD8_Obj <- readRDS("CD8.rds")
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
combined<-wrap_plots(plot_list, ncol = 5)
ggsave("CD8_gene_expression_rearranged.svg", combined_plot, width = 16, height = 12, device = "svg")
