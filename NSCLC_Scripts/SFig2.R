library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(patchwork)

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
liu_cd8_file <- file.path(intermediate_dir, "liu_cd8_mapped.rds")
require_input(cd8_reference_file)
if (!file.exists(liu_cd8_file)) {
  stop(
    "Mapped CD8 object is missing: ",
    liu_cd8_file,
    "\nRun NSCLC_Scripts/Fig2.R first."
  )
}
require_input(liu_cd8_file)

s.genes  <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

# --- Load data objects ---
CD8_Obj <- readRDS(cd8_reference_file)
liu_mapped <- readRDS(liu_cd8_file)

# --- A. Query (Liu) Processing ---
liu_cd8 <- liu_mapped
if (!"Phase" %in% colnames(liu_cd8@meta.data)) {
  liu_cd8 <- CellCycleScoring(
    liu_cd8,
    s.features = s.genes,
    g2m.features = g2m.genes,
    set.ident = FALSE
  )
}

# Reference (CD8_Obj) Processing (Includes UMAP model fix) ---
CD8_Obj <- NormalizeData(CD8_Obj)
CD8_Obj <- FindVariableFeatures(CD8_Obj)
CD8_Obj <- ScaleData(CD8_Obj)
CD8_Obj <- RunPCA(CD8_Obj)

CD8_Obj <- RunUMAP(
  object = CD8_Obj,
  reduction = "pca",
  dims = 1:20,
  return.model = TRUE 
)
# Add Cell Cycle Scoring to the Reference object
CD8_Obj <- CellCycleScoring(CD8_Obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)

# CRITICAL FIX: Add the predicted labels to the original query object's metadata
if ("predicted.Ref_cluster" %in% colnames(liu_mapped@meta.data)) {
  liu_cd8$predicted.cluster <- liu_mapped$predicted.Ref_cluster
} else if ("predicted.predicted.celltype" %in% colnames(liu_mapped@meta.data)) {
  liu_cd8$predicted.cluster <- liu_mapped$predicted.predicted.celltype
} else {
  stop("Mapped CD8 object has no predicted cluster metadata column.")
}

# Extract UMAP embeddings
ref_umap <- Embeddings(CD8_Obj, "umap")
liu_umap <- Embeddings(liu_mapped, "ref.umap") # This reduction is created by MapQuery


# 3. MERGE OBJECTS AND CREATE UNIFIED REDUCTIONS/METADATA
# Merge objects (liu_cd8 now contains the predicted labels and Phase info)
combined <- merge(
  x = CD8_Obj,
  y = liu_cd8,
  add.cell.ids = c("Reference", "Liu"),
  merge.data = TRUE
)

# Add dataset label to metadata
combined$dataset <- ifelse(startsWith(Cells(combined), "Reference_"), "Reference", "Liu")

# Create a single 'Cluster' column (coalesce takes the first non-NA value)
# Ref cells get 'cell.type', Liu cells get 'predicted.cluster'
combined$Cluster <- coalesce(combined$cell.type, combined$predicted.cluster)

# Prepare and combine UMAP embeddings 
ref_umap_prefixed <- ref_umap
rownames(ref_umap_prefixed) <- paste0("Reference_", rownames(ref_umap_prefixed))
liu_umap_prefixed <- liu_umap
rownames(liu_umap_prefixed) <- paste0("Liu_", rownames(liu_umap_prefixed))
all_umap <- rbind(ref_umap_prefixed, liu_umap_prefixed)

# Reorder UMAP rows to match combined object cell order
all_umap_ordered <- all_umap[Cells(combined), , drop = FALSE]

# Create a DimReduc object and attach
umap_reduc <- CreateDimReducObject(
  embeddings = as.matrix(all_umap_ordered),
  key = "UMAP_",
  assay = DefaultAssay(combined)
)
combined[["umap"]] <- umap_reduc


# Prepare and combine PCA embeddings (for the PCA plot) 
ref_pca <- Embeddings(CD8_Obj, "pca")[, 1:20]
liu_pca <- Embeddings(liu_mapped, "ref.pca")[, 1:20]
rownames(ref_pca) <- paste0("Reference_", rownames(ref_pca))
rownames(liu_pca) <- paste0("Liu_", rownames(liu_pca))
all_pca <- rbind(ref_pca, liu_pca)
all_pca_ordered <- all_pca[Cells(combined), , drop = FALSE]

pca_reduc <- CreateDimReducObject(
  embeddings = as.matrix(all_pca_ordered),
  key = "PC_",
  assay = DefaultAssay(combined)
)
combined[["pca"]] <- pca_reduc

df_cycle_by_dataset <- combined@meta.data %>%
  count(dataset, Phase, name = "n") %>%
  group_by(dataset) %>%
  mutate(freq = n / sum(n)) %>%
  ungroup()

liu_cycle <- liu_cd8@meta.data %>%
  count(cluster, Phase, name = "n") %>%
  group_by(cluster) %>%
  mutate(freq = n / sum(n)) %>%
  ungroup() %>%
  rename(liu_cluster = cluster)

#FIG A
p_cycle_dataset <- ggplot(df_cycle_by_dataset, aes(x = dataset, y = freq, fill = Phase)) +
  geom_bar(stat = "identity", position = "fill") +
  ylab("Proportion of cells") + xlab("Dataset") +
  ggtitle("Cell Cycle Composition by Dataset") +
  theme_bw() +
  scale_fill_manual(values = c("G1" = "#1f77b4", "S" = "#ff7f0e", "G2M" = "#2ca02c")) +
  theme(plot.title = element_text(hjust = 0.5))
p_cycle_dataset
ggsave(
  file.path(results_dir, "SFig2_cell_cycle_by_dataset.pdf"),
  plot = p_cycle_dataset,
  width = 7,
  height = 5
)

#FIG B
p_cycle_liu <- ggplot(
  liu_cycle,
  aes(x = liu_cluster, y = freq, fill = Phase)
) +
  geom_col(position = position_dodge(width = 0.8)) +
  labs(
    title = "Cell-cycle distribution by original Liu clusters",
    x = "Liu cluster",
    y = "Proportion of cells"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

p_cycle_liu
ggsave(
  file.path(results_dir, "SFig2_cell_cycle_by_Liu_cluster.pdf"),
  plot = p_cycle_liu,
  width = 9,
  height = 6
)

#FIG C
