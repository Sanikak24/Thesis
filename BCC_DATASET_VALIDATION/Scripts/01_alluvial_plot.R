suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(Matrix)
  library(FNN)
  library(ggplot2)
  library(scales)
  library(ggalluvial)
  library(gridExtra)
  library(cowplot)
})

# --- Custom helper: Add missing genes with zero expression for Seurat v5 ---
AddMissingGenes <- function(seurat_obj, gene_list, assay = "RNA") {
  current_genes <- rownames(seurat_obj)
  missing_genes <- setdiff(gene_list, current_genes)
  
  if (length(missing_genes) > 0) {
    zero_mat <- Matrix::Matrix(
      0,
      nrow = length(missing_genes),
      ncol = ncol(seurat_obj),
      dimnames = list(missing_genes, colnames(seurat_obj)),
      sparse = TRUE
    )
    
    counts <- GetAssayData(seurat_obj, assay = assay, layer = "counts")
    counts <- rbind(counts, zero_mat)
    seurat_obj <- SetAssayData(seurat_obj, assay = assay, layer = "counts", new.data = counts)
    
    data <- GetAssayData(seurat_obj, assay = assay, layer = "data")
    data <- rbind(data, zero_mat)
    seurat_obj <- SetAssayData(seurat_obj, assay = assay, layer = "data", new.data = data)
  }
  
  return(seurat_obj)
}

# ============================================================
# 0) INPUT OBJECTS YOU ALREADY HAVE FROM YOUR PIPELINE
# ============================================================
# CD8_ref     <- readRDS("CD8.rds")   # atlas reference
# yost_cd8    <- ...                 # Yost CD8-only Seurat object
# (or) yost_mapped <- ...            # if you already mapped, you can skip re-mapping below

# ============================================================
# 1) Preprocess reference atlas (needed for PCA + UMAP model)
# ============================================================
CD8_ref <- CD8_ref %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA()

# UMAP model for MapQuery projection
CD8_ref <- RunUMAP(CD8_ref, reduction = "pca", dims = 1:20, return.model = TRUE)

stopifnot("cell.type" %in% colnames(CD8_ref@meta.data))

# ============================================================
# 2) Preprocess Yost CD8 query
# ============================================================
yost_cd8 <- yost_cd8 %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData()

# ============================================================
# 3) Add missing atlas genes to Yost (so features = all_genes works)
# ============================================================
all_genes <- rownames(CD8_ref)
yost_cd8  <- AddMissingGenes(yost_cd8, all_genes)

# ============================================================
# 4) Find anchors using ALL atlas genes + MapQuery
# ============================================================
anchors_yost <- FindTransferAnchors(
  reference = CD8_ref,
  query = yost_cd8,
  dims = 1:20,
  features = all_genes,
  reference.reduction = "pca"
)

yost_mapped <- MapQuery(
  anchorset = anchors_yost,
  query = yost_cd8,
  reference = CD8_ref,
  refdata = list(predicted.celltype = "cell.type"),  # name can be anything
  reference.reduction = "pca",
  reduction.model = "umap"
)

# MapQuery will create:
# yost_mapped$predicted.celltype
# yost_mapped$predicted.celltype.score
# plus ref.pca + ref.umap reductions

# ============================================================
# 5) Majority vote label (exactly like Liu)
# ============================================================
ref_pca   <- Embeddings(CD8_ref, reduction = "pca")[, 1:20, drop = FALSE]
query_pca <- Embeddings(yost_mapped, reduction = "ref.pca")[, 1:20, drop = FALSE]

neighbors  <- get.knnx(data = ref_pca, query = query_pca, k = 10)$nn.index
ref_labels <- CD8_ref$cell.type

majority_vote <- apply(neighbors, 1, function(idx) {
  votes <- ref_labels[idx]
  tab <- sort(table(votes), decreasing = TRUE)
  if (length(tab) == 0 || max(tab) < 5) return(NA_character_)
  names(tab)[1]
})

yost_mapped$majority_vote_pca <- majority_vote

# ============================================================
# 6) Build label map + legend table (Predicted state -> number)
# ============================================================
predicted_labels <- unique(na.omit(yost_mapped$majority_vote_pca))
label_map <- setNames(seq_along(predicted_labels), predicted_labels)

predicted_label_legend <- data.frame(
  Predicted_Cluster_Num  = seq_along(predicted_labels),
  Predicted_Cluster_Name = predicted_labels
)

legend_table <- tableGrob(
  predicted_label_legend,
  rows = NULL,
  theme = ttheme_minimal(
    core = list(fg_params = list(fontsize = 8)),
    colhead = list(fg_params = list(fontsize = 9, fontface = "bold"))
  )
)

# ============================================================
# FIG C (like your Liu Fig C):
# "Proportion of original Yost clusters within each mapped atlas cluster"
# ============================================================
figC_data_yost <- yost_mapped@meta.data %>%
  mutate(
    majority_vote_pca = ifelse(is.na(majority_vote_pca), "NA", as.character(majority_vote_pca)),
    Yost_Cluster = as.character(cluster)   # original Yost cluster column
  ) %>%
  group_by(majority_vote_pca, Yost_Cluster) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(majority_vote_pca) %>%
  mutate(
    prop = n / sum(n),
    label_text = ifelse(prop > 0.05, paste0(round(prop * 100, 0), "%"), "")
  )

p_figC_yost <- ggplot(figC_data_yost, aes(x = majority_vote_pca, y = prop, fill = Yost_Cluster)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  geom_text(aes(label = label_text),
            position = position_stack(vjust = 0.5),
            size = 3, color = "white", fontface = "bold") +
  labs(
    title = "Proportion of original Yost phenotypes within each mapped atlas state",
    x = "Predicted atlas state (majority vote)",
    y = "Proportion of cells",
    fill = "Original Yost cluster"
  ) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid.major.x = element_blank(),
    legend.position = "right"
  )

print(p_figC_yost)

# ============================================================
# FIG D Alluvial (like your Liu Fig D):
# Yost cluster -> Predicted atlas cluster #
# ============================================================
yost_mapped$predicted_label_num <- as.character(label_map[as.character(yost_mapped$majority_vote_pca)])

summary_df_yost <- yost_mapped@meta.data %>%
  mutate(
    Yost_Cluster = as.character(cluster),
    Predicted_Label_Num = ifelse(is.na(predicted_label_num), "NA", predicted_label_num)
  ) %>%
  count(Yost_Cluster, Predicted_Label_Num, name = "Freq")

p_alluvial_yost <- ggplot(summary_df_yost,
                          aes(axis1 = Yost_Cluster, axis2 = Predicted_Label_Num, y = Freq)) +
  geom_alluvium(aes(fill = Yost_Cluster), width = 1/12, alpha = 0.8) +
  geom_stratum(width = 1/12, fill = "gray95", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_x_discrete(
    limits = c("Yost Cluster", "Predicted Cluster #"),
    expand = c(0.1, 0.1)
  ) +
  labs(
    title = "Alluvial Plot: Yost clusters → Predicted atlas clusters",
    x = NULL, y = "Number of Cells", fill = "Yost Cluster"
  ) +
  theme_minimal(base_size = 14)

final_plot_yost <- plot_grid(p_alluvial_yost, legend_table, rel_widths = c(2.5, 1), nrow = 1)
print(final_plot_yost)
#ggsave(filename = "outputs_yost_to_atlas/Alluvial_Plot.pdf",plot = p_alluvial_yost, width = 9, height = 12)
