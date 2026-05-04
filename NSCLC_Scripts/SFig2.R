s.genes  <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

# --- Load data objects ---
CD8_Obj <- readRDS("CD8.rds")
liu_counts <- readRDS("GSE179994_all.Tcell.rawCounts.rds/GSE179994_all.Tcell.rawCounts.rds")
liu_meta <- readr::read_tsv("GSE179994_Tcell.metadata.tsv/GSE179994_Tcell.metadata.tsv")

# --- A. Query (Liu) Processing ---
liu_meta_cd8 <- liu_meta %>% filter(celltype == "CD8")
counts_cd8   <- liu_counts[, liu_meta_cd8$cellid]
liu_cd8 <- CreateSeuratObject(counts = counts_cd8, meta.data = as.data.frame(liu_meta_cd8))

liu_cd8 <- NormalizeData(liu_cd8)
liu_cd8 <- FindVariableFeatures(liu_cd8)
liu_cd8 <- ScaleData(liu_cd8)
liu_cd8 <- RunPCA(liu_cd8, features = VariableFeatures(liu_cd8))
# Add Cell Cycle Scoring to the Query object (needed for the bar plot)
liu_cd8 <- CellCycleScoring(liu_cd8, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)

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

# 2. INTEGRATION (MAPPING)

# Find anchors (reference.reduction must be set)
anchors <- FindTransferAnchors(
  reference = CD8_Obj,
  query = liu_cd8,
  dims = 1:20,
  reference.reduction = "pca"
)

# Map Liu into atlas UMAP space
liu_mapped <- MapQuery(
  anchorset = anchors,
  query = liu_cd8,
  reference = CD8_Obj,
  # Use 'cell.type' from the reference, predicted column will be 'predicted.Ref_cluster'
  refdata = list(Ref_cluster = CD8_Obj$cell.type), 
  reduction.model = "umap"
)

# CRITICAL FIX: Add the predicted labels to the original query object's metadata
# The actual name of the predicted column in liu_mapped will be 'predicted.Ref_cluster'
liu_cd8$predicted.cluster <- liu_mapped$predicted.Ref_cluster

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

#FIG A
p_cycle_dataset <- ggplot(df_cycle_by_dataset, aes(x = dataset, y = freq, fill = Phase)) +
  geom_bar(stat = "identity", position = "fill") +
  ylab("Proportion of cells") + xlab("Dataset") +
  ggtitle("Cell Cycle Composition by Dataset") +
  theme_bw() +
  scale_fill_manual(values = c("G1" = "#1f77b4", "S" = "#ff7f0e", "G2M" = "#2ca02c")) +
  theme(plot.title = element_text(hjust = 0.5))
p_cycle_dataset

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

#FIG C
