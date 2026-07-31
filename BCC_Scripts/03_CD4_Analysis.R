# ============================================================
# Project : sigNATURE Framework
# Script  : # ============================================================
# Project : sigNATURE Framework
# Script  : 03_CD4_Analysis.R

# Purpose : Map Yost et al. BCC CD8 T cells to the reference
#           CD8 T-cell atlas using Seurat label transfer.
#
# Dataset :
#   - Yost et al. (2019), GSE123813
#   - MD Anderson CD8 T-cell Atlas
#
# Input :
#   - BCC Seurat object
#   - Reference CD8 atlas
#
# Output :
#   - Integrated Seurat object
#   - Predicted atlas labels
#
# Author : Sanika Kamath
# ============================================================
# ============================================================
# FULL END-TO-END SCRIPT: Yost CD4 (Tfh + Treg) -> CD4 Atlas
# (runs from DATA directory)
#
# Produces (in ONE folder):
#   1) Original Yost UMAP using AUTHOR coordinates (UMAP1/UMAP2) for CD4 subset
#   2) Scree plots (CD4 atlas + Yost CD4) from PCA stdev
#   3) Yost CD4 projected into CD4 atlas UMAP (ref.umap)
#   4) Projected/predicted atlas-state UMAP (colored by predicted atlas label)
#   5) Projection confidence plot (predicted.Ref_cluster.score) in ref.umap space
#   6) Overlay: atlas + Yost CD4 points in one UMAP space
#   7) Saves RDS objects
#
# Key points:
#   - Filters clinical to Tumor Type == "BCC"
#   - Subsets Yost to CD4 clusters (Tfh + Treg) BEFORE mapping/overlay
#   - Creates output dirs up front (no ggsave prompts)
#   - Confidence plot robust to UMAP colname variants
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(readxl)
  library(Seurat)
  library(ggplot2)
  library(Matrix)
})

# Run from the project root:
#   Rscript BCC_Scripts/03_CD4_Analysis.R
target_vsize_mb <- 65536
if (is.finite(mem.maxVSize()) && mem.maxVSize() < target_vsize_mb) {
  invisible(mem.maxVSize(target_vsize_mb))
}

# -----------------------------
# 0) Paths 
# -----------------------------
counts_file <- file.path("data", "bcc", "GSE123813_bcc_scRNA_counts.txt.gz")
meta_t_file <- file.path("data", "bcc", "GSE123813_bcc_tcell_metadata.txt.gz")
clin_file   <- file.path("data", "bcc", "41591_2019_522_MOESM2_ESM.xlsx")

cd4_atlas_rds <- file.path("data", "reference", "CD4_Obj_for_mapping.rds")
cd8_mapped_rds <- file.path(
  "results", "bcc", "atlas_integration", "Yost_to_CD8Atlas_mapped.rds"
)

# output directories
resdir <- file.path("results", "bcc", "cd4_analysis")
dir.create(resdir, showWarnings = FALSE, recursive = TRUE)

required_files <- c(
  counts_file, meta_t_file, clin_file, cd4_atlas_rds, cd8_mapped_rds
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required input file(s): ", paste(missing_files, collapse = ", "))
}

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
  df <- df[df$PC <= ndims, , drop = FALSE]
  ggplot(df, aes(x = PC, y = Stdev)) +
    geom_point(size = 1.5) +
    labs(title = title, x = "PC", y = "Standard deviation") +
    theme_classic()
}

# -----------------------------
# Helper: safe UMAP column naming for ggplot
# -----------------------------
standardize_umap_cols <- function(df) {
  cn <- colnames(df)
  
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
# 1) LOAD COUNTS (53,030 cells) + T-CELL METADATA (33,106)
# ============================================================
cat("\n[1] Loading counts...\n")
counts_dt <- fread(counts_file)
gene_col   <- counts_dt[[1]]
counts_mat <- as.matrix(counts_dt[, -1])
rownames(counts_mat) <- gene_col
rm(counts_dt); gc()

cat("Counts matrix dim (genes x cells): ", nrow(counts_mat), " x ", ncol(counts_mat), "\n", sep="")

cat("\n[2] Loading T-cell metadata...\n")
meta_tcell <- fread(meta_t_file)
cat("T-cell metadata dim: ", nrow(meta_tcell), " x ", ncol(meta_tcell), "\n", sep="")

# ============================================================
# 2) LOAD CLINICAL + FILTER BCC ONLY
# ============================================================
cat("\n[3] Loading clinical table + filtering to BCC...\n")
supp <- read_excel(clin_file, skip = 2)

clin_bcc <- supp %>%
  filter(!is.na(Patient)) %>%
  filter(`Tumor Type` == "BCC") %>%
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

# save all T-cells author UMAP (optional)
p_umap_all <- DimPlot(obj, reduction = "umap", group.by = "cluster", label = TRUE) +
  ggtitle("Yost BCC T-cells: Original author UMAP (by cluster)") +
  theme_classic()

#ggsave(file.path(resdir, "Yost_original_author_umap_allTcells.pdf"),
       #p_umap_all, width = 9, height = 6)

# ============================================================
# 5) SUBSET TO CD4 (Tfh + Treg) BEFORE ANY MAPPING/OVERLAY
# ============================================================
cat("\n[6] Subsetting to CD4 clusters: Tfh + Treg...\n")
cd4_clusters <- c("Tfh", "Tregs")
yost_cd4 <- subset(obj, subset = cluster %in% cd4_clusters)

report_obj(yost_cd4, "Yost CD4 only (raw counts + author UMAP)")

p_umap_yost_cd4 <- DimPlot(yost_cd4, reduction = "umap", group.by = "cluster", label = TRUE) +
  ggtitle("Yost CD4: Original author UMAP (Tfh + Treg only)") +
  theme_classic()

#ggsave(file.path(resdir, "Yost_original_author_umap_CD4_Tfh_Treg.pdf"),
       #p_umap_yost_cd4, width = 9, height = 6)

# ============================================================
# 6) LOAD CD4 ATLAS + PREPROCESS BOTH DATASETS
# ============================================================
cat("\n[7] Loading CD4 atlas reference...\n")
CD4_ref <- readRDS(cd4_atlas_rds)
report_obj(CD4_ref, "CD4 atlas (loaded)")

# Confirm the atlas label column exists
stopifnot("cell.type" %in% colnames(CD4_ref@meta.data))

# ---- preprocess Yost CD4
cat("\n[8] Preprocessing Yost CD4 (Normalize → HVG → Scale → PCA)...\n")
yost_cd4 <- NormalizeData(yost_cd4)
yost_cd4 <- FindVariableFeatures(yost_cd4)
yost_cd4 <- ScaleData(yost_cd4)
yost_cd4 <- RunPCA(yost_cd4, features = VariableFeatures(yost_cd4))
report_obj(yost_cd4, "Yost CD4 (after PCA)")

# ---- preprocess CD4 atlas
cat("\n[9] Preprocessing CD4 atlas (Normalize → HVG → Scale → PCA → UMAP model)...\n")
CD4_ref <- NormalizeData(CD4_ref)
CD4_ref <- FindVariableFeatures(CD4_ref)
CD4_ref <- ScaleData(CD4_ref)
CD4_ref <- RunPCA(CD4_ref)

# Scree plots
p_scree_atlas <- plot_scree(CD4_ref, "CD4 atlas: scree plot (PCA stdev)", ndims = 50)
#ggsave(file.path(resdir, "CD4Atlas_scree_plot.pdf"), p_scree_atlas, width = 7, height = 4)

p_scree_yost <- plot_scree(yost_cd4, "Yost CD4 (Tfh+Treg): scree plot (PCA stdev)", ndims = 50)
#ggsave(file.path(resdir, "YostCD4_scree_plot.pdf"), p_scree_yost, width = 7, height = 4)

# UMAP model needed for MapQuery
CD4_ref <- RunUMAP(CD4_ref, reduction = "pca", dims = 1:20, return.model = TRUE)
report_obj(CD4_ref, "CD4 atlas (after UMAP model)")

# ============================================================
# 7) ANCHORS + MAPQUERY (PROJECT YOST CD4 INTO ATLAS SPACE)
# ============================================================
cat("\n[10] Finding transfer anchors + mapping Yost CD4 into CD4 atlas...\n")

anchors <- FindTransferAnchors(
  reference = CD4_ref,
  query = yost_cd4,
  dims = 1:20,
  reference.reduction = "pca"
)

yost_cd4_mapped <- MapQuery(
  anchorset = anchors,
  query = yost_cd4,
  reference = CD4_ref,
  refdata = list(Ref_cluster = "cell.type"),
  reduction.model = "umap"
)

report_obj(yost_cd4_mapped, "Yost CD4 mapped (ref.umap present)")

# Predicted atlas state in atlas coordinates
p_pred <- DimPlot(
  yost_cd4_mapped,
  reduction = "ref.umap",
  group.by = "predicted.Ref_cluster",
  label = TRUE
) +
  ggtitle("Yost CD4 projected into CD4 atlas UMAP (colored by predicted atlas state)") +
  theme_classic()

#ggsave(file.path(resdir, "YostCD4_projected_predicted_clusters_refUMAP.pdf"),
       #p_pred, width = 9, height = 6)

# ============================================================
# 8) PROJECTION CONFIDENCE PLOT
# ============================================================
cat("\n[11] Plotting projection confidence (predicted.Ref_cluster.score)...\n")

yost_embed <- as.data.frame(Embeddings(yost_cd4_mapped, "ref.umap"))
yost_embed <- standardize_umap_cols(yost_embed)
yost_embed$score <- yost_cd4_mapped@meta.data$predicted.Ref_cluster.score

p_score <- ggplot(yost_embed, aes(x = UMAP_1, y = UMAP_2, color = score)) +
  geom_point(size = 0.4) +
  labs(
    title = "CD4 Projection confidence (predicted.Ref_cluster.score)",
    x = "UMAP_1", y = "UMAP_2", color = "score"
  ) +
  theme_classic()

#ggsave(file.path(resdir, "YostCD4_projection_score_refUMAP.pdf"),
       #p_score, width = 9, height = 6)

# ============================================================
# 9) OVERLAY: YOST CD4 POINTS ON CD4 ATLAS UMAP
# ============================================================
cat("\n[12] Creating overlay (CD4 atlas + Yost CD4) in ONE UMAP space...\n")

ref_umap  <- Embeddings(CD4_ref, "umap")
yost_umap <- Embeddings(yost_cd4_mapped, "ref.umap")

combined <- merge(
  x = CD4_ref,
  y = yost_cd4,
  add.cell.ids = c("Reference","Yost"),
  merge.data = TRUE
)
report_obj(combined, "Combined (CD4 atlas + Yost CD4)")

# prefix embedding rownames to match merged cell names
ref_umap_pref <- ref_umap
rownames(ref_umap_pref) <- paste0("Reference_", rownames(ref_umap_pref))

yost_umap_pref <- yost_umap
rownames(yost_umap_pref) <- paste0("Yost_", rownames(yost_umap_pref))

all_umap <- rbind(ref_umap_pref, yost_umap_pref)

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
  ggtitle("Overlay: Yost CD4 (Tfh+Treg) projected onto CD4 atlas UMAP") +
  theme_classic()

ggsave(file.path(resdir, "Overlay_YostCD4_on_CD4Atlas.pdf"),
       p_overlay, width = 9, height = 6)

# ============================================================
# 10) SAVE RDS OBJECTS
# ============================================================
saveRDS(yost_cd4,        file.path(resdir, "Yost_CD4_Tfh_Treg_only.rds"))
saveRDS(CD4_ref,         file.path(resdir, "CD4_atlas_processed.rds"))
saveRDS(yost_cd4_mapped, file.path(resdir, "Yost_CD4_mapped_to_CD4Atlas.rds"))

cat("\nDONE.\nAll figures + RDS saved in:\n  ", normalizePath(resdir), "\n", sep = "")

#ALLUVIAL PLOT
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(ggalluvial)
  library(gridExtra)
  library(cowplot)
})

# -----------------------------
# 0) Load objects if needed
# -----------------------------
# CD4_ref        <- readRDS("CD4.rds")
# yost_cd4       <- readRDS("outputs_yost_to_cd4_atlas/RESULTS/Yost_CD4_Tfh_Treg_only.rds")
# yost_cd4_mapped<- readRDS("outputs_yost_to_cd4_atlas/RESULTS/Yost_CD4_mapped_to_CD4Atlas.rds")

# Continue using the project-level output directory defined above.
dir.create(resdir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# 1) Pick which label column to use for alluvial
# -----------------------------
# Prefer majority vote if you created it; otherwise use MapQuery predicted labels
label_col <- NULL
if ("majority_vote_pca" %in% colnames(yost_cd4_mapped@meta.data)) {
  label_col <- "majority_vote_pca"
} else if ("predicted.Ref_cluster" %in% colnames(yost_cd4_mapped@meta.data)) {
  label_col <- "predicted.Ref_cluster"
} else if ("predicted.celltype" %in% colnames(yost_cd4_mapped@meta.data)) {
  label_col <- "predicted.celltype"
} else {
  stop("No predicted label column found. Available cols:\n",
       paste(colnames(yost_cd4_mapped@meta.data), collapse = ", "))
}

cat("Using label column:", label_col, "\n")

# -----------------------------
# 2) Build numeric label map + legend (Predicted state -> number)
# -----------------------------
pred_levels <- sort(unique(na.omit(as.character(yost_cd4_mapped[[label_col]][,1]))))
label_map <- setNames(seq_along(pred_levels), pred_levels)

legend_df <- data.frame(
  Predicted_Label_Num  = seq_along(pred_levels),
  Predicted_Label_Name = pred_levels
)

legend_table <- tableGrob(
  legend_df,
  rows = NULL,
  theme = ttheme_minimal(
    core    = list(fg_params = list(fontsize = 8)),
    colhead = list(fg_params = list(fontsize = 9, fontface = "bold"))
  )
)

# -----------------------------
# 3) FIG C-style stacked proportions:
# "Proportion of original Yost clusters within each predicted atlas state"
# -----------------------------
figC_df <- yost_cd4_mapped@meta.data %>%
  mutate(
    Yost_Cluster = as.character(cluster),
    Predicted_State = as.character(.data[[label_col]]),
    Predicted_State = ifelse(is.na(Predicted_State), "NA", Predicted_State)
  ) %>%
  count(Predicted_State, Yost_Cluster, name = "n") %>%
  group_by(Predicted_State) %>%
  mutate(
    prop = n / sum(n),
    label_text = ifelse(prop > 0.05, paste0(round(prop * 100, 0), "%"), "")
  ) %>%
  ungroup()

p_figC <- ggplot(figC_df, aes(x = Predicted_State, y = prop, fill = Yost_Cluster)) +
  geom_bar(stat = "identity", width = 0.8) +
  geom_text(
    aes(label = label_text),
    position = position_stack(vjust = 0.5),
    size = 3, color = "white", fontface = "bold"
  ) +
  labs(
    title = "CD4: Proportion of Yost phenotypes within each predicted atlas state",
    x = "Predicted atlas state",
    y = "Proportion of cells",
    fill = "Yost cluster"
  ) +
  scale_y_continuous(labels = percent_format()) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid.major.x = element_blank()
  )

ggsave(file.path(resdir, "CD4_FigC_stacked_proportions.pdf"),
       p_figC, width = 10, height = 6)

# -----------------------------
# 4) FIG D alluvial:
# Yost cluster -> Predicted atlas cluster #
# -----------------------------
alluvial_df <- yost_cd4_mapped@meta.data %>%
  mutate(
    Yost_Cluster = as.character(cluster),
    Predicted_State = as.character(.data[[label_col]]),
    Predicted_State = ifelse(is.na(Predicted_State), "NA", Predicted_State),
    Predicted_Label_Num = ifelse(Predicted_State == "NA", "NA", as.character(label_map[Predicted_State]))
  ) %>%
  count(Yost_Cluster, Predicted_Label_Num, name = "Freq")

p_alluvial <- ggplot(alluvial_df,
                     aes(axis1 = Yost_Cluster, axis2 = Predicted_Label_Num, y = Freq)) +
  geom_alluvium(aes(fill = Yost_Cluster), width = 1/12, alpha = 0.8) +
  geom_stratum(width = 1/12, fill = "gray95", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_x_discrete(
    limits = c("Yost Cluster", "Predicted Cluster #"),
    expand = c(0.1, 0.1)
  ) +
  labs(
    title = "CD4 Alluvial: Yost clusters → predicted CD4 atlas states (#)",
    x = NULL, y = "Number of cells",
    fill = "Yost Cluster"
  ) +
  theme_minimal(base_size = 14)

final_plot <- plot_grid(p_alluvial, legend_table, rel_widths = c(2.6, 1), nrow = 1)

ggsave(file.path(resdir, "CD4_Alluvial_with_legend.pdf"),
       final_plot, width = 12, height = 7)

# -----------------------------
# 5) Save counts table + matrix view
# -----------------------------
#write.csv(alluvial_df,
          #file.path(resdir, "CD4_alluvial_counts.csv"),
          #row.names = FALSE)

print(alluvial_df)
print(xtabs(Freq ~ Yost_Cluster + Predicted_Label_Num, data = alluvial_df))

# Save legend too
#write.csv(legend_df,
          #file.path(resdir, "CD4_predicted_label_num_legend.csv"),
          #row.names = FALSE)

cat("\nSaved:\n",
    "- CD4_FigC_stacked_proportions.pdf\n",
    "- CD4_Alluvial_with_legend.pdf\n",
    "- CD4_alluvial_counts.csv\n",
    "- CD4_predicted_label_num_legend.csv\n",
    "in: ", normalizePath(resdir), "\n", sep = "")


##REGRESSION ANALYSIS

# ============================================================
# FIG.A — Yost BCC CD4: 3 MODELS (Yost / Atlas / Combined)
#   - KEEP "NA" atlas label (do not drop any labels)
#   - DO NOT include n as a predictor
#   - ggroc aesthetics + fixed colors (purple/orange/yellow)
#   - Print predictors in each model
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(pROC)
  library(ggplot2)
})

# -----------------------------
# 0) Choose atlas label column for CD4
#    Prefer majority vote if you created it; otherwise use MapQuery labels
# -----------------------------
if ("majority_vote_pca" %in% colnames(yost_cd4_mapped@meta.data)) {
  atlas_col <- "majority_vote_pca"
} else if ("predicted.Ref_cluster" %in% colnames(yost_cd4_mapped@meta.data)) {
  atlas_col <- "predicted.Ref_cluster"
} else if ("predicted.predicted.Ref_cluster" %in% colnames(yost_cd4_mapped@meta.data)) {
  atlas_col <- "predicted.predicted.Ref_cluster"
} else {
  stop("No predicted atlas label column found in yost_cd4_mapped@meta.data.\nCols:\n",
       paste(colnames(yost_cd4_mapped@meta.data), collapse = ", "))
}
cat("Using atlas label column:", atlas_col, "\n")

# -----------------------------
# 1) Build cell-level meta (CD4)
#    Response from response_binary already joined in your pipeline
# -----------------------------
meta_cd4 <- yost_cd4_mapped@meta.data %>%
  mutate(
    SampleID        = as.character(patient),
    Yost_Cluster    = as.character(cluster),                 # Tfh / Tregs
    Atlas_Cluster   = as.character(.data[[atlas_col]]),      # keep "NA" as label
    Response_Status = as.character(response_binary)
  ) %>%
  filter(!is.na(SampleID)) %>%
  filter(Response_Status %in% c("Responder", "NonResponder")) %>%
  mutate(Response_Bin = ifelse(Response_Status == "Responder", 1, 0))

sample_status_cd4 <- meta_cd4 %>%
  select(SampleID, Response_Status, Response_Bin) %>%
  distinct()

# -----------------------------
# 2) Helper: fit proportion-based glm + ROC (NO n predictor)
#    Fix levels/direction => stable AUC from pROC
# -----------------------------
fit_prop_model <- function(meta, sample_status, cluster_col) {
  
  props <- meta %>%
    count(SampleID, .data[[cluster_col]]) %>%
    group_by(SampleID) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup() %>%
    pivot_wider(
      names_from  = all_of(cluster_col),
      values_from = prop,
      values_fill = 0
    )
  
  df <- props %>% left_join(sample_status, by = "SampleID")
  
  drop_cols <- c("SampleID", "Response_Status", "Response_Bin", "n")
  features  <- setdiff(colnames(df), drop_cols)
  features  <- setdiff(features, "n")
  
  formula <- as.formula(
    paste("Response_Bin ~", paste0("`", features, "`", collapse = " + "))
  )
  
  model <- glm(formula, data = df, family = binomial)
  
  roc_obj <- roc(
    response  = df$Response_Bin,
    predictor = as.numeric(predict(model, type = "response")),
    levels    = c(0, 1),
    direction = "<",
    quiet     = TRUE
  )
  
  list(df=df, features=features, formula=formula, model=model, roc=roc_obj)
}

# -----------------------------
# 3) Model 1 — Yost clusters only
# -----------------------------
m1_cd4 <- fit_prop_model(meta_cd4, sample_status_cd4, "Yost_Cluster")
roc_m1 <- m1_cd4$roc

# -----------------------------
# 4) Model 2 — Atlas clusters only (NA INCLUDED)
# -----------------------------
m2_cd4 <- fit_prop_model(meta_cd4, sample_status_cd4, "Atlas_Cluster")
roc_m2 <- m2_cd4$roc

# -----------------------------
# 5) Model 3 — Combined Yost × Atlas (NA INCLUDED)
# -----------------------------
meta_comb <- meta_cd4 %>%
  mutate(combined_cluster = paste(Yost_Cluster, Atlas_Cluster, sep = "_"))

props_m3 <- meta_comb %>%
  count(SampleID, combined_cluster) %>%
  group_by(SampleID) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  pivot_wider(
    names_from  = combined_cluster,
    values_from = prop,
    values_fill = 0
  )

df_m3 <- props_m3 %>% left_join(sample_status_cd4, by = "SampleID")

drop_cols <- c("SampleID", "Response_Status", "Response_Bin", "n")
features_m3 <- setdiff(colnames(df_m3), drop_cols)
features_m3 <- setdiff(features_m3, "n")

formula_m3 <- as.formula(
  paste("Response_Bin ~", paste0("`", features_m3, "`", collapse = " + "))
)

model3 <- glm(formula_m3, data = df_m3, family = binomial)

roc_m3 <- roc(
  response  = df_m3$Response_Bin,
  predictor = as.numeric(predict(model3, type = "response")),
  levels    = c(0, 1),
  direction = "<",
  quiet     = TRUE
)

m3_cd4 <- list(df=df_m3, features=features_m3, formula=formula_m3, model=model3, roc=roc_m3)

# -----------------------------
# 6) Print predictors
# -----------------------------
cat("\n====================\nCD4 PREDICTORS\n====================\n")

cat("\nModel 1 (Yost) predictors:\n")
print(m1_cd4$features)

cat("\nModel 2 (Atlas) predictors:\n")
print(m2_cd4$features)

cat("\nModel 3 (Combined) predictors:\n")
print(m3_cd4$features)

# -----------------------------
# 7) ROC plot — ggroc aesthetics + fixed colors
#    purple = Yost, orange = Atlas, yellow = Combined
# -----------------------------
label_m1 <- paste0("Model 1 (Yost): AUC = ", round(auc(roc_m1), 3))
label_m2 <- paste0("Model 2 (Atlas): AUC = ", round(auc(roc_m2), 3))
label_m3 <- paste0("Model 3 (Combined): AUC = ", round(auc(roc_m3), 3))

roc_list <- list(
  "Model 1 (Yost)"     = roc_m1,
  "Model 2 (Atlas)"    = roc_m2,
  "Model 3 (Combined)" = roc_m3
)

cols <- c(
  "Model 1 (Yost)"     = "#7D3C98",  # purple
  "Model 2 (Atlas)"    = "#E67E22",  # orange
  "Model 3 (Combined)" = "#F1C40F"   # yellow
)

p_cd4_figA <- ggroc(roc_list, legacy.axes = TRUE, size = 1.1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey40") +
  scale_color_manual(
    name   = "Model Performance",
    values = cols,
    breaks = c("Model 1 (Yost)", "Model 2 (Atlas)", "Model 3 (Combined)"),
    labels = c(label_m1, label_m2, label_m3)
  ) +
  labs(
    title = "Yost BCC CD4 — Sample-level prediction of response",
    x = "Specificity (1 − FPR)",
    y = "Sensitivity (TPR)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.title    = element_text(face = "bold")
  )

print(p_cd4_figA)

ggsave(
  file.path(resdir, "Regression_plot_BCC_CD4.pdf"),
  p_cd4_figA,
  width = 12,
  height = 9,
  dpi = 300
)

## ALL PREDICTORS
# ============================================================
# FIG.A — COMBINED CD8 + CD4 RESPONSE PREDICTION (3 MODELS)
#   Model 1: Yost (CD8_Yost + CD4_Yost)
#   Model 2: Atlas (CD8_Atlas + CD4_Atlas)  [keeps "NA"]
#   Model 3: Combined (CD8_YostxAtlas + CD4_YostxAtlas) [keeps "NA"]
#   - NO "n" as predictor (ever)
#   - Guards against many-to-many joins (1 row per SampleID)
#   - ggroc aesthetics + fixed colors: purple/orange/yellow
#   - Prints predictors for each model
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(pROC)
  library(ggplot2)
})

# -----------------------------
# 0) Pick atlas label columns
# -----------------------------
yost_mapped <- readRDS(cd8_mapped_rds)
if ("majority_vote_pca" %in% colnames(yost_mapped@meta.data)) {
  cd8_atlas_col <- "majority_vote_pca"
} else if ("predicted.Ref_cluster" %in% colnames(yost_mapped@meta.data)) {
  cd8_atlas_col <- "predicted.Ref_cluster"
} else {
  stop("CD8 mapped object has no majority-vote or predicted atlas label column.")
}

if ("majority_vote_pca" %in% colnames(yost_cd4_mapped@meta.data)) {
  cd4_atlas_col <- "majority_vote_pca"
} else if ("predicted.Ref_cluster" %in% colnames(yost_cd4_mapped@meta.data)) {
  cd4_atlas_col <- "predicted.Ref_cluster"
} else if ("predicted.predicted.Ref_cluster" %in% colnames(yost_cd4_mapped@meta.data)) {
  cd4_atlas_col <- "predicted.predicted.Ref_cluster"
} else {
  stop("CD4 object yost_cd4_mapped has no predicted atlas label column.\nCols:\n",
       paste(colnames(yost_cd4_mapped@meta.data), collapse = ", "))
}

cat("CD8 atlas col:", cd8_atlas_col, "\n")
cat("CD4 atlas col:", cd4_atlas_col, "\n")

# -----------------------------
# 1) Build cell-level meta
# -----------------------------
meta_cd8 <- yost_mapped@meta.data %>%
  mutate(
    SampleID        = as.character(patient),
    Response_Status = as.character(response_binary),
    Response_Bin    = ifelse(Response_Status == "Responder", 1, 0),
    Yost_Cluster    = as.character(cluster),
    Atlas_Cluster   = as.character(.data[[cd8_atlas_col]])  # KEEP "NA" labels
  ) %>%
  filter(!is.na(SampleID)) %>%
  filter(Response_Status %in% c("Responder","NonResponder"))

meta_cd4 <- yost_cd4_mapped@meta.data %>%
  mutate(
    SampleID        = as.character(patient),
    Response_Status = as.character(response_binary),
    Response_Bin    = ifelse(Response_Status == "Responder", 1, 0),
    Yost_Cluster    = as.character(cluster),               # Tfh/Tregs
    Atlas_Cluster   = as.character(.data[[cd4_atlas_col]]) # KEEP "NA" labels
  ) %>%
  filter(!is.na(SampleID)) %>%
  filter(Response_Status %in% c("Responder","NonResponder"))

# Stable sample set: intersection CD8 & CD4
sample_status <- meta_cd8 %>%
  select(SampleID, Response_Status, Response_Bin) %>% distinct() %>%
  inner_join(
    meta_cd4 %>% select(SampleID, Response_Status, Response_Bin) %>% distinct(),
    by = "SampleID",
    suffix = c(".cd8", ".cd4")
  ) %>%
  transmute(
    SampleID,
    Response_Status = Response_Status.cd8,
    Response_Bin    = Response_Bin.cd8
  ) %>%
  distinct()

keep_ids <- sample_status$SampleID
meta_cd8 <- meta_cd8 %>% filter(SampleID %in% keep_ids)
meta_cd4 <- meta_cd4 %>% filter(SampleID %in% keep_ids)

cat("Samples used (intersection CD8 & CD4):", length(keep_ids), "\n")

# -----------------------------
# 2) Helper: make prop table (1 row per SampleID) + prefix feature names
#    IMPORTANT: never keep raw count column 'n' as a feature
# -----------------------------
make_prop_table <- function(meta, cluster_col, prefix) {
  
  props <- meta %>%
    count(SampleID, .data[[cluster_col]], name = ".cell_n") %>%
    group_by(SampleID) %>%
    mutate(prop = .cell_n / sum(.cell_n)) %>%
    ungroup() %>%
    select(-.cell_n) %>%
    pivot_wider(
      names_from  = all_of(cluster_col),
      values_from = prop,
      values_fill = 0
    )
  
  # prefix feature columns (avoid CD8/CD4 collisions)
  feat_cols <- setdiff(colnames(props), "SampleID")
  props <- props %>% rename_with(~ paste0(prefix, .x), all_of(feat_cols))
  
  # guarantee 1 row per SampleID (prevents many-to-many join warnings)
  props <- props %>%
    group_by(SampleID) %>%
    summarise(across(everything(), ~ ifelse(is.numeric(.x), mean(.x), dplyr::first(.x))), .groups = "drop")
  
  # just in case some label is literally called "n" in future: drop prefixed __n
  props <- props %>% select(-matches("__n$"))
  
  props
}

fit_glm_from_df <- function(df, sample_status) {
  
  # enforce identical sample set
  df2 <- sample_status %>%
    inner_join(df, by = "SampleID") %>%
    mutate(across(where(is.numeric), ~ replace_na(.x, 0)))
  
  # drop non-features
  drop_cols <- c("SampleID","Response_Status","Response_Bin")
  features  <- setdiff(colnames(df2), drop_cols)
  
  # extra safety: remove any feature ending in "__n"
  features  <- features[!grepl("__n$", features)]
  
  formula <- as.formula(
    paste("Response_Bin ~", paste0("`", features, "`", collapse = " + "))
  )
  
  model <- glm(formula, data = df2, family = binomial)
  
  roc_obj <- roc(
    response  = df2$Response_Bin,
    predictor = as.numeric(predict(model, type = "response")),
    levels    = c(0, 1),
    direction = "<",
    quiet     = TRUE
  )
  
  list(df=df2, features=features, formula=formula, model=model, roc=roc_obj)
}

# -----------------------------
# 3) Build the three combined designs
# -----------------------------
# Model 1: Yost only
df_m1 <- make_prop_table(meta_cd8, "Yost_Cluster", "CD8_Yost__") %>%
  inner_join(make_prop_table(meta_cd4, "Yost_Cluster", "CD4_Yost__"), by = "SampleID")
m1 <- fit_glm_from_df(df_m1, sample_status)

# Model 2: Atlas only (keep "NA")
df_m2 <- make_prop_table(meta_cd8, "Atlas_Cluster", "CD8_Atlas__") %>%
  inner_join(make_prop_table(meta_cd4, "Atlas_Cluster", "CD4_Atlas__"), by = "SampleID")
m2 <- fit_glm_from_df(df_m2, sample_status)

# Model 3: Combined (Yost×Atlas) (keep "NA")
meta_cd8_c <- meta_cd8 %>% mutate(combined_cluster = paste(Yost_Cluster, Atlas_Cluster, sep = "_"))
meta_cd4_c <- meta_cd4 %>% mutate(combined_cluster = paste(Yost_Cluster, Atlas_Cluster, sep = "_"))

df_m3 <- make_prop_table(meta_cd8_c, "combined_cluster", "CD8_Comb__") %>%
  inner_join(make_prop_table(meta_cd4_c, "combined_cluster", "CD4_Comb__"), by = "SampleID")
m3 <- fit_glm_from_df(df_m3, sample_status)

# -----------------------------
# 4) Print predictors (NO n)
# -----------------------------
cat("\n====================\nCOMBINED CD8+CD4 PREDICTORS (NO n)\n====================\n")
cat("\nModel 1 (Yost) predictors:\n"); print(m1$features)
cat("\nModel 2 (Atlas) predictors:\n"); print(m2$features)
cat("\nModel 3 (Combined) predictors:\n"); print(m3$features)

# -----------------------------
# 5) ggroc-style ROC plot with fixed colors
# -----------------------------
label_m1 <- paste0("Model 1 (Yost): AUC = ", round(auc(m1$roc), 3))
label_m2 <- paste0("Model 2 (Atlas): AUC = ", round(auc(m2$roc), 3))
label_m3 <- paste0("Model 3 (Combined): AUC = ", round(auc(m3$roc), 3))

roc_list <- list(
  "Model 1 (Yost)"     = m1$roc,
  "Model 2 (Atlas)"    = m2$roc,
  "Model 3 (Combined)" = m3$roc
)

cols <- c(
  "Model 1 (Yost)"     = "#7D3C98",  # purple
  "Model 2 (Atlas)"    = "#E67E22",  # orange
  "Model 3 (Combined)" = "#F1C40F"   # yellow
)

p_figA_cd8_cd4 <- ggroc(roc_list, legacy.axes = TRUE, size = 1.1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey40") +
  scale_color_manual(
    name   = "Model Performance",
    values = cols,
    breaks = c("Model 1 (Yost)", "Model 2 (Atlas)", "Model 3 (Combined)"),
    labels = c(label_m1, label_m2, label_m3)
  ) +
  labs(
    title = "Yost BCC — Sample-level prediction of response (CD8 + CD4 predictors)",
    x = "Specificity (1 − FPR)",
    y = "Sensitivity (TPR)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.title    = element_text(face = "bold")
  )

print(p_figA_cd8_cd4)

ggsave(
  file.path(resdir, "Regression_plot_CD8plusCD4_Yost_Atlas_Combined.pdf"),
  p_figA_cd8_cd4,
  width = 12,
  height = 9,
  dpi = 300
)
