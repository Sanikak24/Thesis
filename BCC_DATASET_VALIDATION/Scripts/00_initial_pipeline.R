# ============================================================
# FIXED, FULL END-TO-END SCRIPT (runs from DATA directory)
#
# What this script produces (in ONE folder):
#   1) Original Yost CD8 UMAP using AUTHOR coordinates (UMAP1/UMAP2)
#   2) Scree plots (atlas + Yost CD8) from PCA stdev
#   3) Yost CD8 projected onto atlas UMAP (overlay: atlas+Yost points)
#   4) Projected/predicted cluster UMAP (ref.umap colored by predicted atlas label)
#   5) Projection confidence plot (predicted.Ref_cluster.score) in ref.umap space
#
# Key fixes vs your error log:
#   - Creates outputs folder + RESULTS subfolder up front (no interactive prompt)
#   - No FeaturePlot(NULL) call
#   - Confidence plot uses correct column names from Embeddings(): UMAP_1 / UMAP_2
#     (we rename columns safely in case they are "UMAP1"/"UMAP2")
#   - Filters clinical data to Tumor Type == "BCC" (your requested fix)
#   - Subsets CD8 BEFORE combining/overlay
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(readxl)
  library(Seurat)
  library(ggplot2)
  library(Matrix)
})

# -----------------------------
# 0) Paths (you said getwd() is .../DATA)
# -----------------------------
counts_file <- "GSE123813_bcc_scRNA_counts.txt.gz"
meta_t_file <- "GSE123813_bcc_tcell_metadata.txt.gz"
clin_file   <- "41591_2019_522_MOESM2_ESM.xlsx"
atlas_rds   <- "CD8.rds"  # reference atlas

# output directories (FIX: create them now so ggsave never prompts)
outdir   <- "outputs_yost_to_atlas"
resdir   <- file.path(outdir, "RESULTS")
dir.create(resdir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Helper: print cells/genes at each step
# -----------------------------
report_obj <- function(obj, label) {
  cat("\n====================\n", label, "\n====================\n", sep = "")
  cat("Cells:", ncol(obj), "\n")
  cat("Genes:", nrow(obj), "\n")
  cat("Assays:", paste(Assays(obj), collapse = ", "), "\n")
  cat("Reductions:", paste(Reductions(obj), collapse = ", "), "\n")
  cat("Meta cols:", ncol(obj@meta.data), "\n")
  invisible(NULL)
}

# -----------------------------
# Helper: Scree plot from PCA stdev
# -----------------------------
plot_scree <- function(seu, title = "Scree plot", ndims = 50) {
  stdev <- seu[["pca"]]@stdev
  df <- data.frame(PC = seq_along(stdev), Stdev = stdev)
  df <- df[df$PC <= ndims, ]
  ggplot(df, aes(x = PC, y = Stdev)) +
    geom_point(size = 1.5) +
    labs(title = title, x = "PC", y = "Standard deviation") +
    theme_classic()
}

# -----------------------------
# Helper: safe UMAP column naming for ggplot
# -----------------------------
standardize_umap_cols <- function(df) {
  # Accepts colnames like: UMAP_1/UMAP_2, UMAP1/UMAP2, umap_1/umap_2, etc.
  cn <- colnames(df)
  # find first two columns that look like umap
  # common cases: "UMAP_1","UMAP_2" or "UMAP1","UMAP2"
  if (all(c("UMAP_1","UMAP_2") %in% cn)) return(df)
  if (all(c("UMAP1","UMAP2") %in% cn)) {
    df <- dplyr::rename(df, UMAP_1 = UMAP1, UMAP_2 = UMAP2)
    return(df)
  }
  if (all(c("umap_1","umap_2") %in% cn)) {
    df <- dplyr::rename(df, UMAP_1 = umap_1, UMAP_2 = umap_2)
    return(df)
  }
  # fallback: assume first 2 columns are UMAP coords
  colnames(df)[1:2] <- c("UMAP_1","UMAP_2")
  df
}

# ============================================================
# 1) LOAD COUNTS (53,030 cells)  +  T-CELL METADATA (33,106)
# ============================================================
cat("\n[1] Loading counts...\n")
counts_dt <- fread(counts_file)
# Warning about 53031 cols is expected: first column is gene names/rownames
gene_col   <- counts_dt[[1]]
counts_mat <- as.matrix(counts_dt[, -1])
rownames(counts_mat) <- gene_col
rm(counts_dt); gc()

cat("Counts matrix dim (genes x cells): ", nrow(counts_mat), " x ", ncol(counts_mat), "\n", sep="")

cat("\n[2] Loading T-cell metadata...\n")
meta_tcell <- fread(meta_t_file)
cat("T-cell metadata dim: ", nrow(meta_tcell), " x ", ncol(meta_tcell), "\n", sep="")

# ============================================================
# 2) LOAD CLINICAL (Supp Table 1) + FILTER BCC ONLY (FIX)
# ============================================================
cat("\n[3] Loading clinical table + filtering to BCC...\n")
supp <- read_excel(clin_file, skip = 2)

clin_bcc <- supp %>%
  filter(!is.na(Patient)) %>%
  filter(`Tumor Type` == "BCC") %>%   # <-- IMPORTANT
  transmute(
    patient = Patient,
    response_raw = Response,
    response_binary = case_when(
      Response %in% c("Yes", "Yes (CR)") ~ "Responder",
      Response == "No" ~ "NonResponder",
      TRUE ~ NA_character_
    )
  )

cat("Clinical rows (BCC-only): ", nrow(clin_bcc), "\n", sep="")
print(table(clin_bcc$response_binary, useNA = "ifany"))

# join to cell-level metadata by patient
meta_tcell_final <- meta_tcell %>%
  left_join(clin_bcc, by = "patient") %>%
  as.data.frame()
rownames(meta_tcell_final) <- meta_tcell_final$cell.id

cat("\nCell-level response counts after join:\n")
print(table(meta_tcell_final$response_binary, useNA = "ifany"))

# ============================================================
# 3) SUBSET COUNTS TO T-CELLS + CREATE SEURAT OBJECT
# ============================================================
cat("\n[4] Subsetting counts to T-cells and creating Seurat object...\n")
common_cells <- intersect(colnames(counts_mat), rownames(meta_tcell_final))
counts_tcell_only <- counts_mat[, common_cells, drop = FALSE]
meta_tcell_final  <- meta_tcell_final[common_cells, , drop = FALSE]

obj <- CreateSeuratObject(
  counts = counts_tcell_only,
  meta.data = meta_tcell_final,
  project = "Yost_BCC_Tcells"
)
report_obj(obj, "Yost T-cells (raw)")

# ============================================================
# 4) ORIGINAL YOST UMAP USING AUTHOR COORDINATES
# ============================================================
cat("\n[5] Attaching author UMAP coordinates (UMAP1/UMAP2) into Seurat 'umap'...\n")
stopifnot(all(c("UMAP1","UMAP2") %in% colnames(obj@meta.data)))

author_umap <- as.matrix(obj@meta.data[, c("UMAP1","UMAP2")])
colnames(author_umap) <- c("UMAP_1","UMAP_2")

obj[["umap"]] <- CreateDimReducObject(
  embeddings = author_umap,
  key = "UMAP_",
  assay = "RNA"
)

p_umap_all <- DimPlot(obj, reduction = "umap", group.by = "cluster", label = TRUE) +
  ggtitle("Yost BCC T-cells: Original author UMAP (by cluster)") +
  theme_classic()

ggsave(file.path(resdir, "Yost_original_author_umap_allTcells.pdf"),
       p_umap_all, width = 9, height = 6)

# ============================================================
# 5) SUBSET TO CD8 BEFORE ANY COMBINING (REQUIRED)
# ============================================================
cat("\n[6] Subsetting to CD8 clusters only...\n")
cd8_clusters <- c("CD8_act","CD8_eff","CD8_ex","CD8_ex_act","CD8_mem")
yost_cd8 <- subset(obj, subset = cluster %in% cd8_clusters)
report_obj(yost_cd8, "Yost CD8 only (raw counts + author UMAP)")

p_umap_yost_cd8 <- DimPlot(yost_cd8, reduction = "umap", group.by = "cluster", label = TRUE) +
  ggtitle("Yost CD8: Original author UMAP (CD8 clusters only)") +
  theme_classic()

ggsave(file.path(resdir, "Yost_original_author_umap_CD8only.pdf"),
       p_umap_yost_cd8, width = 9, height = 6)

# ============================================================
# 6) LOAD ATLAS + PREPROCESS BOTH DATASETS (PCA/UMAP MODEL)
# ============================================================
cat("\n[7] Loading CD8 atlas reference...\n")
CD8_ref <- readRDS(atlas_rds)
report_obj(CD8_ref, "CD8 atlas (loaded)")

stopifnot("cell.type" %in% colnames(CD8_ref@meta.data))

# ---- preprocess Yost CD8
cat("\n[8] Preprocessing Yost CD8 (Normalize → HVG → Scale → PCA)...\n")
yost_cd8 <- NormalizeData(yost_cd8)
yost_cd8 <- FindVariableFeatures(yost_cd8)
yost_cd8 <- ScaleData(yost_cd8)
yost_cd8 <- RunPCA(yost_cd8, features = VariableFeatures(yost_cd8))
report_obj(yost_cd8, "Yost CD8 (after PCA)")

# ---- preprocess atlas
cat("\n[9] Preprocessing atlas (Normalize → HVG → Scale → PCA → UMAP model)...\n")
CD8_ref <- NormalizeData(CD8_ref)
CD8_ref <- FindVariableFeatures(CD8_ref)
CD8_ref <- ScaleData(CD8_ref)
CD8_ref <- RunPCA(CD8_ref)

# Scree plots (FIX: always valid)
p_scree_atlas <- plot_scree(CD8_ref, "CD8 atlas: scree plot (PCA stdev)", ndims = 50)
#ggsave(file.path(resdir, "Atlas_scree_plot.pdf"), p_scree_atlas, width = 7, height = 4)

p_scree_yost <- plot_scree(yost_cd8, "Yost CD8: scree plot (PCA stdev)", ndims = 50)
#ggsave(file.path(resdir, "YostCD8_scree_plot.pdf"), p_scree_yost, width = 7, height = 4)

# UMAP model (needed for MapQuery)
CD8_ref <- RunUMAP(CD8_ref, reduction = "pca", dims = 1:20, return.model = TRUE)
report_obj(CD8_ref, "CD8 atlas (after UMAP model)")

# ============================================================
# 7) ANCHORS + MAPQUERY (PROJECT YOST CD8 INTO ATLAS SPACE)
# ============================================================
cat("\n[10] Finding transfer anchors + mapping Yost CD8 into atlas...\n")
anchors <- FindTransferAnchors(
  reference = CD8_ref,
  query = yost_cd8,
  dims = 1:20,
  reference.reduction = "pca"
)

yost_mapped <- MapQuery(
  anchorset = anchors,
  query = yost_cd8,
  reference = CD8_ref,
  refdata = list(Ref_cluster = "cell.type"),
  reduction.model = "umap"
)
report_obj(yost_mapped, "Yost CD8 mapped (ref.umap present)")

# Projected/predicted cluster UMAP in atlas coordinates
p_pred <- DimPlot(
  yost_mapped,
  reduction = "ref.umap",
  group.by = "predicted.Ref_cluster",
  label = FALSE
) +
  ggtitle("Yost CD8 projected into atlas UMAP (colored by predicted atlas state)") +
  theme_classic()

#ggsave(file.path(resdir, "Yost_projected_predicted_clusters_refUMAP.pdf"),
       #p_pred, width = 9, height = 6)

# ============================================================
# 8) PROJECTION CONFIDENCE PLOT (FIXED UMAP_1/UMAP_2 NAMES)
# ============================================================
cat("\n[11] Plotting projection confidence (predicted.Ref_cluster.score)...\n")
yost_embed <- as.data.frame(Embeddings(yost_mapped, "ref.umap"))
yost_embed <- standardize_umap_cols(yost_embed)
yost_embed$score <- yost_mapped@meta.data$predicted.Ref_cluster.score

p_score <- ggplot(yost_embed, aes(x = UMAP_1, y = UMAP_2, color = score)) +
  geom_point(size = 0.4) +
  labs(
    title = "Projection confidence (predicted.Ref_cluster.score)",
    x = "UMAP_1", y = "UMAP_2", color = "score"
  ) +
  theme_classic()

#ggsave(file.path(resdir, "Yost_projection_score_refUMAP.pdf"),
       #p_score, width = 9, height = 6)

# ============================================================
# 9) OVERLAY: YOST CD8 POINTS ON ATLAS UMAP (COMBINE CD8 ONLY)
# ============================================================
cat("\n[12] Creating overlay (atlas + Yost CD8) in ONE UMAP space...\n")
ref_umap  <- Embeddings(CD8_ref, "umap")
yost_umap <- Embeddings(yost_mapped, "ref.umap")  # yost coords in atlas UMAP space

combined <- merge(
  x = CD8_ref,
  y = yost_cd8,                 # <-- CD8 only (as you requested)
  add.cell.ids = c("Reference","Yost"),
  merge.data = TRUE
)
report_obj(combined, "Combined (atlas + Yost CD8)")

# prefix embedding rownames to match merged cell names
ref_umap_pref <- ref_umap
rownames(ref_umap_pref) <- paste0("Reference_", rownames(ref_umap_pref))

yost_umap_pref <- yost_umap
rownames(yost_umap_pref) <- paste0("Yost_", rownames(yost_umap_pref))

all_umap <- rbind(ref_umap_pref, yost_umap_pref)

# order to match combined cells
combined_cells <- Cells(combined)
stopifnot(all(combined_cells %in% rownames(all_umap)))
all_umap_ord <- all_umap[combined_cells, , drop = FALSE]

combined[["umap"]] <- CreateDimReducObject(
  embeddings = as.matrix(all_umap_ord),
  key = "UMAP_",
  assay = DefaultAssay(combined)
)

combined$dataset <- ifelse(startsWith(Cells(combined), "Reference_"), "Reference", "Yost")

p_overlay <- DimPlot(
  combined,
  reduction = "umap",
  group.by = "dataset",
  pt.size = 0.8,
  raster = TRUE
) +
  ggtitle("Overlay: Yost CD8 projected onto CD8 atlas UMAP") +
  theme_classic()

#ggsave(file.path(resdir, "Overlay_Yost_on_Atlas.pdf"),
       #p_overlay, width = 9, height = 6)

# ============================================================
# 10) SAVE RDS OBJECTS
# ============================================================
#saveRDS(obj,         file.path(resdir, "Yost_BCC_Tcells_all.rds"))
#saveRDS(yost_cd8,    file.path(resdir, "Yost_CD8_only.rds"))
#saveRDS(CD8_ref,     file.path(resdir, "CD8_atlas_processed.rds"))
#saveRDS(yost_mapped, file.path(resdir, "Yost_CD8_mapped_to_atlas.rds"))
#saveRDS(combined,    file.path(resdir, "Atlas_plus_Yost_overlay.rds"))

cat("\nDONE.\nAll figures + RDS saved in:\n  ", normalizePath(resdir), "\n", sep = "")


library(dplyr)

dataset_table <- tibble::tibble(
  Dataset = c(
    "Yost T-cells (all)",
    "Yost CD8 only",
    "CD8 Atlas reference",
    "Combined atlas + Yost CD8"
  ),
  Cells = c(33106, 12364, 110218, 122582),
  Genes = c(23309, 23309, 70777, 72911)
)

dataset_table

