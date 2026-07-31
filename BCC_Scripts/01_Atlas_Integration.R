# ============================================================
# Project : sigNATURE Framework
# Script  : 01_Atlas_Integration.R
# Purpose : Map Yost et al. BCC CD8 T cells to the reference
#           CD8 T-cell atlas using Seurat label transfer.
#
# Datasets:
#   - Yost et al. (2019), GSE123813
#   - MD Anderson CD8 T-cell Atlas
#
# Inputs:
#   - data/bcc/GSE123813_bcc_scRNA_counts.txt.gz
#       Raw single-cell RNA-seq count matrix from Yost et al.
#
#   - data/bcc/GSE123813_bcc_tcell_metadata.txt.gz
#       Cell-level metadata for the Yost BCC dataset.
#
#   - data/bcc/41591_2019_522_MOESM2_ESM.xlsx
#       Clinical response metadata from the Nature Medicine
#       supplementary files.
#
#   - data/reference/CD8_Obj_for_mapping.rds
#       MD Anderson CD8 T-cell reference atlas as a Seurat object.
#
# Outputs:
#   - results/bcc/atlas_integration/Yost_BCC_Tcells_Final.rds
#       Processed Seurat object for the Yost BCC T-cell dataset.
#
#   - results/bcc/atlas_integration/Yost_to_CD8Atlas_mapped.rds
#       CD8 cells after label transfer from the reference atlas.
#
#   - results/bcc/atlas_integration/Yost_CD8_ATLAS_PROJECTION.pdf
#       Overlay of mapped Yost CD8 cells on the MD Anderson
#       CD8 T-cell atlas.
#
# Author : Sanika Kamath
# ============================================================

# The reference/query merge near the end requires more than macOS R's
# default 16 GB vector-heap limit. Raise the ceiling before loading data so
# this script can be run directly with:
#   Rscript BCC_Scripts/01_Atlas_Integration.R
target_vsize_mb <- 65536
current_vsize_mb <- mem.maxVSize()
if (is.finite(current_vsize_mb) && current_vsize_mb < target_vsize_mb) {
  mem.maxVSize(target_vsize_mb)
}
message("R vector-memory limit: ", mem.maxVSize(), " MB")

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(data.table)
  library(dplyr)
  library(readxl)
  library(ggplot2)
  library(Matrix)
})

# 0) PATHS
bcc_data_dir <- file.path("data", "bcc")
ref_data_dir <- file.path("data", "reference")
out_dir      <- file.path("results", "bcc", "atlas_integration")

counts_file <- Sys.getenv(
  "BCC_COUNTS_FILE",
  unset = file.path(bcc_data_dir, "GSE123813_bcc_scRNA_counts.txt.gz")
)
meta_file <- Sys.getenv(
  "BCC_METADATA_FILE",
  unset = file.path(bcc_data_dir, "GSE123813_bcc_tcell_metadata.txt.gz")
)
clin_xlsx <- Sys.getenv(
  "BCC_CLINICAL_FILE",
  unset = file.path(bcc_data_dir, "41591_2019_522_MOESM2_ESM.xlsx")
)
ref_rds <- Sys.getenv(
  "CD8_REFERENCE_RDS",
  unset = file.path(ref_data_dir, "CD8_Obj_for_mapping.rds")
)

out_rds    <- file.path(out_dir, "Yost_BCC_Tcells_Final.rds")
out_mapped <- file.path(out_dir, "Yost_to_CD8Atlas_mapped.rds")
out_pdf    <- file.path(out_dir, "Yost_CD8_ATLAS_PROJECTION.pdf")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(out_dir)) {
  stop("Could not create output directory: ", out_dir)
}

input_files <- c(counts_file, meta_file, clin_xlsx, ref_rds)
missing_files <- input_files[!file.exists(input_files)]
if (length(missing_files) > 0) {
  stop("Missing required input file(s): ",
       paste(missing_files, collapse = ", "))
}

# Metadata column in the atlas that I want to transfer
ref_label_col <- "cell.type"

# Yost clusters that define CD8 states in the paper metadata
cd8_clusters <- c("CD8_act","CD8_eff","CD8_ex","CD8_ex_act","CD8_mem")

message("Loading BCC count matrix...")
# 1) LOAD COUNTS
counts_dt <- fread(counts_file)
# first column is gene names 
gene_col <- counts_dt[[1]]
counts_mat <- as.matrix(counts_dt[, -1, with = FALSE])
rownames(counts_mat) <- gene_col

rm(counts_dt); gc()

# OPTIONAL but recommended:
# Seurat replaces "_" with "-" in feature names. Do it explicitly here
# so intersections across objects behave predictably.
rownames(counts_mat) <- gsub("_", "-", rownames(counts_mat))

# Convert to sparse for memory
counts_mat <- Matrix(counts_mat, sparse = TRUE)

# 2) LOAD T-CELL METADATA + CLINICAL RESPONSE
message("Loading BCC metadata and clinical response data...")

meta_tcell <- fread(meta_file)
required_meta_cols <- c("cell.id", "patient", "cluster")
missing_meta_cols <- setdiff(required_meta_cols, colnames(meta_tcell))
if (length(missing_meta_cols) > 0) {
  stop("T-cell metadata is missing required column(s): ",
       paste(missing_meta_cols, collapse = ", "))
}

supp <- read_excel(clin_xlsx, skip = 2)
required_clin_cols <- c("Patient", "Response")
missing_clin_cols <- setdiff(required_clin_cols, colnames(supp))
if (length(missing_clin_cols) > 0) {
  stop("Clinical spreadsheet is missing required column(s): ",
       paste(missing_clin_cols, collapse = ", "))
}

clin <- supp %>%
  filter(!is.na(Patient)) %>%
  transmute(
    patient = Patient,
    response_binary = case_when(
      Response %in% c("Yes", "Yes (CR)") ~ "Responder",
      Response == "No" ~ "NonResponder",
      TRUE ~ NA_character_
    )
  )

meta_tcell_final <- meta_tcell %>%
  left_join(clin, by = "patient") %>%
  as.data.frame()

rownames(meta_tcell_final) <- meta_tcell_final$cell.id

# 3) ALIGN CELLS + CREATE SEURAT OBJECT
common_cells <- intersect(colnames(counts_mat), rownames(meta_tcell_final))
message("Common cells: ", length(common_cells))

if (length(common_cells) == 0) {
  stop("No common cells between counts and metadata. Check cell IDs.")
}

counts_tcell_only <- counts_mat[, common_cells]
meta_tcell_final  <- meta_tcell_final[common_cells, , drop = FALSE]

obj <- CreateSeuratObject(
  counts    = counts_tcell_only,
  meta.data = meta_tcell_final,
  project   = "Yost_BCC_Tcells"
)

print(obj)
print(table(obj$response_binary, useNA="ifany"))
print(table(obj$cluster, useNA="ifany"))

# 4) ADD Yost UMAP (from metadata columns UMAP1/UMAP2)
if (all(c("UMAP1","UMAP2") %in% colnames(obj@meta.data))) {
  author_umap <- as.matrix(obj@meta.data[, c("UMAP1","UMAP2")])
  colnames(author_umap) <- c("UMAP_1","UMAP_2")
  
  obj[["umap"]] <- CreateDimReducObject(
    embeddings = author_umap,
    key = "UMAP_",
    assay = DefaultAssay(obj)
  )
  
  p1 <- DimPlot(obj, reduction="umap", group.by="cluster", label=TRUE) +
    ggtitle("Yost BCC: Original Author UMAP")
  print(p1)
} else {
  message("UMAP1/UMAP2 not found in metadata; skipping author UMAP import.")
}

message("Saving processed BCC T-cell object...")
saveRDS(obj, out_rds)

# 5) LOAD CD8 REFERENCE ATLAS
message("Loading CD8 reference atlas...")
CD8_ref <- readRDS(ref_rds)

if (!(ref_label_col %in% colnames(CD8_ref@meta.data))) {
  stop("ref_label_col '", ref_label_col, "' not found in CD8_ref@meta.data. Available: ",
       paste(colnames(CD8_ref@meta.data), collapse=", "))
}

# Keep feature name convention consistent
rownames(CD8_ref) <- gsub("_", "-", rownames(CD8_ref))

# 6) SUBSET YOST TO CD8 CLUSTERS
message("Subsetting Yost to CD8 clusters...")
yost_cd8 <- subset(obj, subset = cluster %in% cd8_clusters)
print(yost_cd8)
if (ncol(yost_cd8) == 0) {
  stop("CD8 cluster subset is empty. Check cluster annotations.")
}

# 7) HARMONIZE GENES BETWEEN REFERENCE + QUERY

common_genes <- intersect(rownames(CD8_ref), rownames(yost_cd8))
message("Common genes: ", length(common_genes))

if (length(common_genes) == 0) {
  stop("No common genes between the CD8 reference and BCC query.")
}

if (length(common_genes) < 2000) {
  warning("Very few common genes found (", length(common_genes), "). Check gene naming conventions.")
}

CD8_ref2  <- subset(CD8_ref, features = common_genes)
yost_cd82 <- subset(yost_cd8, features = common_genes)

# 8) NORMALIZE / PCA / UMAP MODEL ON REFERENCE (needed for MapQuery projection)
message("Preprocessing query and reference...")

# Query processing
yost_cd82 <- NormalizeData(yost_cd82)
yost_cd82 <- FindVariableFeatures(yost_cd82)
yost_cd82 <- ScaleData(yost_cd82)
yost_cd82 <- RunPCA(yost_cd82, features = VariableFeatures(yost_cd82))

# Reference processing (ensure it has PCA + UMAP model)
CD8_ref2 <- NormalizeData(CD8_ref2)
CD8_ref2 <- FindVariableFeatures(CD8_ref2)
CD8_ref2 <- ScaleData(CD8_ref2)
CD8_ref2 <- RunPCA(CD8_ref2)

# IMPORTANT: return.model=TRUE so query can be projected into reference UMAP
CD8_ref2 <- RunUMAP(
  CD8_ref2,
  reduction = "pca",
  dims = 1:15,
  return.model = TRUE
)

# 9) TRANSFER LABELS + MAP QUERY
message("Finding transfer anchors and mapping the query...")

anchors <- FindTransferAnchors(
  reference = CD8_ref2,
  query = yost_cd82,
  dims = 1:15,
  reference.reduction = "pca"
)

yost_mapped <- MapQuery(
  anchorset = anchors,
  query = yost_cd82,
  reference = CD8_ref2,
  refdata = list(Ref_cluster = CD8_ref2@meta.data[[ref_label_col]]),
  reduction.model = "umap"
)

# 10) PLOTS
# Predicted labels in meta.data:
# predicted.Ref_cluster and predicted.Ref_cluster.score
grep("predicted|score", colnames(yost_mapped@meta.data), value = TRUE) |> print()

p2 <- DimPlot(
  yost_mapped,
  reduction = "ref.umap",
  group.by = "predicted.Ref_cluster",
  label = TRUE
) + ggtitle("Yost CD8 projected onto CD8 atlas (predicted states)")
print(p2)

# Plot prediction confidence: it's meta.data, not an expression feature.
yost_mapped$pred_score <- yost_mapped$predicted.Ref_cluster.score

p3 <- FeaturePlot(
  yost_mapped,
  reduction = "ref.umap",
  features = "pred_score"
) + ggtitle("Prediction confidence (predicted.Ref_cluster.score)")
print(p3)

# High-confidence subset
yost_hi <- subset(yost_mapped, subset = predicted.Ref_cluster.score >= 0.5)

p4 <- DimPlot(
  yost_hi,
  reduction = "ref.umap",
  group.by = "predicted.Ref_cluster",
  label = TRUE
) + ggtitle("High-confidence projected cells (score ≥ 0.5)")
print(p4)

# 11) OPTIONAL: OVERLAY REFERENCE + YOST ON SAME UMAP
message("Creating reference and mapped-query overlay...")

ref_umap <- Embeddings(CD8_ref2, "umap")
yost_umap <- Embeddings(yost_mapped, "ref.umap")  # projection coordinates

required_mapped_cols <- c(
  "predicted.Ref_cluster",
  "predicted.Ref_cluster.score"
)
missing_mapped_cols <- setdiff(
  required_mapped_cols,
  colnames(yost_mapped@meta.data)
)
if (length(missing_mapped_cols) > 0) {
  stop("Mapped query is missing required metadata column(s): ",
       paste(missing_mapped_cols, collapse = ", "))
}

if (!identical(rownames(yost_umap), Cells(yost_mapped))) {
  stop("Projected query embeddings do not align with mapped-query cell names.")
}

combined <- merge(
  x = CD8_ref2,
  y = yost_mapped,
  add.cell.ids = c("Reference", "Yost"),
  merge.data = FALSE
)

ref_umap_pref <- ref_umap
rownames(ref_umap_pref) <- paste0("Reference_", rownames(ref_umap_pref))

yost_umap_pref <- yost_umap
rownames(yost_umap_pref) <- paste0("Yost_", rownames(yost_umap_pref))

all_umap <- rbind(ref_umap_pref, yost_umap_pref)

if (!identical(sort(rownames(all_umap)), sort(Cells(combined)))) {
  stop("Projected embeddings do not align with merged cell names.")
}

all_umap <- all_umap[Cells(combined), , drop = FALSE]

if (!identical(rownames(all_umap), Cells(combined))) {
  stop("Projected embeddings could not be ordered to match merged cells.")
}

combined[["umap"]] <- CreateDimReducObject(
  embeddings = as.matrix(all_umap),
  key = "UMAP_",
  assay = DefaultAssay(combined)
)

combined$dataset <- ifelse(startsWith(Cells(combined), "Reference_"), "Reference", "Yost")

p5 <- DimPlot(combined, reduction="umap", group.by="dataset", pt.size=1) +
  ggtitle("Yost CD8 overlaid on CD8 atlas UMAP")
print(p5)

message("Saving atlas-integration outputs...")
ggsave(filename = out_pdf, plot = p5, width = 9, height = 12)
saveRDS(yost_mapped, out_mapped)
message("Atlas integration complete. Outputs written to: ", out_dir)
