library(Seurat)
library(ggplot2)
library(dplyr)
library(readr)
library(patchwork)

s.genes  <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

# Load raw data and subset to CD8 
liu_counts <- readRDS("GSE179994_all.Tcell.rawCounts.rds/GSE179994_all.Tcell.rawCounts.rds")
liu_meta <- readr::read_tsv("GSE179994_Tcell.metadata.tsv/GSE179994_Tcell.metadata.tsv", show_col_types = FALSE)

liu_meta_cd8 <- liu_meta %>% filter(celltype == "CD8")
counts_cd8   <- liu_counts[, liu_meta_cd8$cellid]

# Create the base Seurat object for Liu CD8 cells
liu_cd8_base <- CreateSeuratObject(counts = counts_cd8, meta.data = as.data.frame(liu_meta_cd8))
DefaultAssay(liu_cd8_base) <- "RNA"
liu_cd8_base <- NormalizeData(liu_cd8_base)
liu_cd8_base <- FindVariableFeatures(liu_cd8_base)
liu_cd8_base <- CellCycleScoring(liu_cd8_base, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)

# 1. BEFORE CELL-CYCLE NORMALIZATION (Uncorrected Data)
message("1. Processing data BEFORE cell-cycle normalization...")

# Copy the base object
liu_cd8_uncorrected <- liu_cd8_base

# Standard scaling and PCA (no cell-cycle regression/correction)
liu_cd8_uncorrected <- ScaleData(liu_cd8_uncorrected)
liu_cd8_uncorrected <- RunPCA(liu_cd8_uncorrected)

# Run UMAP on the uncorrected PCA
liu_cd8_uncorrected <- RunUMAP(
  liu_cd8_uncorrected,
  reduction = "pca",
  dims = 1:20
)

# Plot 1: UMAP Before Correction, colored by Cell Cycle Phase
p_before_phase <- DimPlot(liu_cd8_uncorrected, reduction = "umap", group.by = "Phase", pt.size = 0.5) +
  ggtitle("UMAP Before CC Normalization (Colored by Phase)")

# Plot 2: UMAP Before Correction, colored by Original Cluster
p_before_cluster <- DimPlot(liu_cd8_uncorrected, reduction = "umap", group.by = "cluster", pt.size = 0.5) +
  ggtitle("UMAP Before CC Normalization (Colored by Liu Cluster)")

# 2. AFTER CELL-CYCLE NORMALIZATION (Corrected Data)
# Copy the base object
liu_cd8_corrected <- liu_cd8_base

# Scaling and PCA WITH cell-cycle regression
# Here, the S.Score and G2M.Score are regressed out of the highly variable features
liu_cd8_corrected <- ScaleData(
  liu_cd8_corrected,
  vars.to.regress = c("S.Score", "G2M.Score")
)
liu_cd8_corrected <- RunPCA(liu_cd8_corrected, features = VariableFeatures(liu_cd8_corrected), npcs = 20)

# Run UMAP on the corrected PCA
liu_cd8_corrected <- RunUMAP(
  liu_cd8_corrected,
  reduction = "pca",
  dims = 1:20
)

#UMAP After Correction, colored by Cell Cycle Phase (to confirm correction)
p_after_phase <- DimPlot(liu_cd8_corrected, reduction = "umap", group.by = "Phase", pt.size = 0.5) +
  ggtitle("UMAP After CC Normalization (Colored by Phase)")

# Fig B: UMAP After Correction, colored by Original Cluster 
p_after_cluster <- DimPlot(liu_cd8_corrected, reduction = "umap", group.by = "cluster", pt.size = 0.5) +
  ggtitle("UMAP After CC Normalization (Colored by Liu Cluster)")

#FIG. C AND D
library(Seurat)
library(tidyverse)
library(FNN)
library(ggalluvial)
library(cowplot)
library(gridExtra)
library(grid)
library(readr)

#NEW
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
CD8_Obj <- readRDS("CD8.rds")
liu_counts <- readRDS("GSE179994_all.Tcell.rawCounts.rds/GSE179994_all.Tcell.rawCounts.rds")
liu_meta   <- read_tsv("GSE179994_Tcell.metadata.tsv/GSE179994_Tcell.metadata.tsv")

liu_meta_cd8 <- liu_meta %>% filter(celltype == "CD8")
counts_cd8   <- liu_counts[, liu_meta_cd8$cellid]

# 2. Create and Normalize Liu CD8 Seurat Object
liu_seurat <- CreateSeuratObject(counts_cd8, meta.data = as.data.frame(liu_meta_cd8)) %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%  # optional; keeps Seurat happy but not used
  ScaleData()
# 3. Preprocess Reference CD8 Seurat Object
CD8_Obj <- CD8_Obj %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  RunUMAP(reduction = "pca", dims = 1:20, return.model = TRUE)

# 4. Add missing genes (with zeros) to Liu object to match reference
all_genes <- rownames(CD8_Obj)
liu_seurat <- AddMissingGenes(liu_seurat, all_genes)

# 5. Find Anchors with All Genes
anchors <- FindTransferAnchors(
  reference = CD8_Obj,
  query = liu_seurat,
  dims = 1:20,
  features = all_genes,
  reference.reduction = "pca"
)

liu_seurat <- MapQuery(
  anchorset = anchors,
  query = liu_seurat,
  reference = CD8_Obj,
  refdata = list(predicted.celltype = "cell.type"),
  reference.reduction = "pca",
  reduction.model = "umap"
)

# 5. Majority Vote Labeling
ref_pca     <- Embeddings(CD8_Obj, reduction = "pca")[, 1:20]
query_pca   <- Embeddings(liu_seurat, reduction = "ref.pca")[, 1:20]
neighbors   <- get.knnx(data = ref_pca, query = query_pca, k = 10)$nn.index
ref_labels  <- CD8_Obj$cell.type
majority_vote <- apply(neighbors, 1, function(idx) {
  votes <- ref_labels[idx]
  tab <- sort(table(votes), decreasing = TRUE)
  if (length(tab) == 0 || max(tab) < 5) return(NA)
  names(tab)[1]
})
liu_seurat$majority_vote_pca <- majority_vote

# 6. Metadata and Label Mapping
liu_meta <- liu_seurat@meta.data %>%
  mutate(Liu_Cluster = cluster, Predicted_Label = majority_vote_pca)
predicted_labels <- unique(liu_meta$Predicted_Label)
label_map <- setNames(seq_along(predicted_labels), predicted_labels)
liu_meta$Predicted_Label_Num <- label_map[liu_meta$Predicted_Label]

# 7. Legend Table
predicted_label_legend <- data.frame(
  Predicted_Cluster_Num = seq_along(predicted_labels),
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


#FIG C
fig2_data <- liu_seurat@meta.data %>%
  # Convert logical NAs to the character string "NA" so ggplot sees it as a category
  mutate(majority_vote_pca = ifelse(is.na(majority_vote_pca), "NA", as.character(majority_vote_pca))) %>%
  # Group by the predicted cluster (including the new "NA" group)
  group_by(majority_vote_pca, cluster) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(majority_vote_pca) %>%
  mutate(
    prop = n / sum(n),
    # Create labels, but only if the slice is big enough to read
    label_text = ifelse(prop > 0.05, paste0(round(prop * 100, 0), "%"), "")
  )

p_fig2 <- ggplot(fig2_data, aes(x = majority_vote_pca, y = prop, fill = cluster)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  geom_text(aes(label = label_text), 
            position = position_stack(vjust = 0.5), 
            size = 3, color = "white", fontface = "bold") +
  labs(
    title = "Proportion of original Liu et al. phenotypes within each mapped TCellMap cluster",
    x = "Predicted TCellMap Cluster",
    y = "Proportion of Cells",
    fill = "Original Liu Cluster"
  ) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid.major.x = element_blank(),
    legend.position = "right"
  )

# 3. Display plot
print(p_fig2)

# FIG. D Alluvial Plot
liu_seurat$predicted_label_num <- as.character(label_map[liu_seurat$majority_vote_pca])
summary_df <- liu_seurat@meta.data %>%
  count(cluster, predicted_label_num, name = "Freq") %>%
  rename(Liu_Cluster = cluster, Predicted_Label_Num = predicted_label_num)

p_alluvial <- ggplot(summary_df, aes(axis1 = Liu_Cluster, axis2 = Predicted_Label_Num, y = Freq)) +
  geom_alluvium(aes(fill = Liu_Cluster), width = 1/12, alpha = 0.8) +
  geom_stratum(width = 1/12, fill = "gray95", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_x_discrete(
    limits = c("Liu Cluster", "Predicted Cluster #"),
    expand = c(0.1, 0.1)
  ) +
  labs(
    title = "Alluvial Plot: Liu Clusters → Predicted CD8 Clusters",
    x = NULL, y = "Number of Cells", fill = "Liu Cluster"
  ) +
  theme_minimal(base_size = 14)

final_plot <- plot_grid(p_alluvial, legend_table, rel_widths = c(2.5, 1), nrow = 1)
print(final_plot)

#FIG E
# ==============================================================================
# GENE SET EXPRESSION ON COMBINED FEATURE SPACE
# ==============================================================================
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(ggrastr)
library(readr)

CD8_Obj <- readRDS("CD8.rds")
liu_counts <- readRDS("GSE179994_all.Tcell.rawCounts.rds/GSE179994_all.Tcell.rawCounts.rds")
liu_meta <- readr::read_tsv("GSE179994_Tcell.metadata.tsv/GSE179994_Tcell.metadata.tsv")

liu_meta_cd8 <- liu_meta %>% filter(celltype == "CD8")
counts_cd8   <- liu_counts[, liu_meta_cd8$cellid]
liu_cd8 <- CreateSeuratObject(counts = counts_cd8, meta.data = as.data.frame(liu_meta_cd8))
liu_cd8 <- NormalizeData(liu_cd8)
liu_cd8 <- FindVariableFeatures(liu_cd8)
liu_cd8 <- ScaleData(liu_cd8)
liu_cd8 <- RunPCA(liu_cd8, features = VariableFeatures(liu_cd8))

CD8_Obj <- NormalizeData(CD8_Obj)        # if not already normalized
CD8_Obj <- FindVariableFeatures(CD8_Obj) # if not already done
CD8_Obj <- ScaleData(CD8_Obj)            # if not already done
CD8_Obj <- RunPCA(CD8_Obj)               # this creates the 'pca' slot
CD8_Obj <- RunUMAP(
  object = CD8_Obj,
  reduction = "pca",      # assuming PCA was used
  dims = 1:20,
  return.model = TRUE     # 🔑 this is the missing piece
)

# 3. Find anchors (reference.reduction must be set)
anchors <- FindTransferAnchors(
  reference = CD8_Obj,
  query = liu_cd8,
  dims = 1:20,
  reference.reduction = "pca"
)

# 4. Map Liu into atlas UMAP space
#    MapQuery will (by default) create a "ref.umap" (or "ref.*") embedding in the mapped object.
liu_mapped <- MapQuery(
  anchorset = anchors,
  query = liu_cd8,
  reference = CD8_Obj,
  refdata = list(Ref_cluster = "cell.type"),  # change if different column name in reference meta
  reduction.model = "umap"
)

# Confirm what embedding name MapQuery produced for projected UMAP:
print(Reductions(liu_mapped))
# usually it's "ref.umap" in the mapped object; check:
if("ref.umap" %in% Reductions(liu_mapped)) {
  liu_umap <- Embeddings(liu_mapped, "ref.umap")
} else if ("umap" %in% Reductions(liu_mapped)) {
  liu_umap <- Embeddings(liu_mapped, "umap")
} else {
  stop("No ref.umap/umap found in liu_mapped; inspect Reductions(liu_mapped).")
}

# Reference UMAP coords
ref_umap <- Embeddings(CD8_Obj, "umap")

# 5. Merge objects (union of genes automatically handled)
#    Use add.cell.ids so we can identify source and match cell name prefixes for embeddings
combined <- merge(
  x = CD8_Obj,
  y = liu_cd8,
  add.cell.ids = c("Reference", "Liu"),
  merge.data = TRUE
)

# 6. Prepare correctly-named UMAP embeddings for the merged object
#    After merge(..., add.cell.ids = c("Reference","Liu")), the merged cellnames are:
#      paste0("Reference_", colnames(CD8_Obj))
#      paste0("Liu_",       colnames(liu_cd8))
#    So we must prefix the ref_umap/liu_umap rownames the same way, then rbind in the same order
# create prefixed rownames
ref_umap_prefixed <- ref_umap
rownames(ref_umap_prefixed) <- paste0("Reference_", rownames(ref_umap_prefixed))

liu_umap_prefixed <- liu_umap
rownames(liu_umap_prefixed) <- paste0("Liu_", rownames(liu_umap_prefixed))

# combine embeddings
all_umap <- rbind(ref_umap_prefixed, liu_umap_prefixed)

# Ensure that the rownames of all_umap exactly match the cell order of combined
combined_cells <- Cells(combined)

# Reorder all_umap rows to match combined (this will error if any mismatch)
if(!all(combined_cells %in% rownames(all_umap))) {
  missing_cells <- setdiff(combined_cells, rownames(all_umap))
  stop("Some merged cells are missing embeddings. Missing example (first 10): ",
       paste(head(missing_cells, 10), collapse = ", "))
}
all_umap_ordered <- all_umap[combined_cells, , drop = FALSE]

# Create a DimReduc object and attach
umap_reduc <- CreateDimReducObject(
  embeddings = as.matrix(all_umap_ordered),
  key = "UMAP_",
  assay = DefaultAssay(combined)
)
combined[["umap"]] <- umap_reduc

# 7. Add a dataset label to metadata so plotting is easy
combined$dataset <- ifelse(startsWith(Cells(combined), "Reference_"), "Reference", "Liu")

# 8. Quick summaries: cells & genes per dataset, and how many cells are visualized in UMAP
get_counts <- function(obj) {
  data.frame(
    Num_Cells = ncol(obj),
    Num_Genes = nrow(obj)
  )
}
ref_counts <- get_counts(CD8_Obj) %>% mutate(Dataset = "Reference")
liu_counts <- get_counts(liu_cd8) %>% mutate(Dataset = "Liu")
combined_counts <- get_counts(combined) %>% mutate(Dataset = "Combined")
summary_counts <- bind_rows(ref_counts, liu_counts, combined_counts) %>% select(Dataset, everything())
print(summary_counts)

# How many cells are present in the UMAP embedding (should equal ncol(combined))
num_umap_cells <- nrow(Embeddings(combined, "umap"))
cat("Cells in UMAP embedding:", num_umap_cells, "\n")

# How many of those are Reference vs Liu (visualized)
table(combined$dataset)
p <- DimPlot(
  combined,
  reduction = "umap",
  group.by = "dataset",
  pt.size = 1
) +
  scale_color_manual(
    values = c("Reference" = "#E64B35",   # dark blue-red for Atlas
               "Liu"       = "#3C5488")   # blue for Liu
  ) +
  ggtitle("Liu CD8 overlaid on Reference CD8 UMAP)") +
  theme_minimal()

print(p)

