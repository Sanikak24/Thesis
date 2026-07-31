
# 1. SETUP & LIBRARIES
library(Seurat)
library(data.table)
library(dplyr)
library(readxl)
library(ggplot2)

# Run from the project root:
#   Rscript BCC_Scripts/05_Original_Dataset_UMAP.R
target_vsize_mb <- 65536
if (is.finite(mem.maxVSize()) && mem.maxVSize() < target_vsize_mb) {
  invisible(mem.maxVSize(target_vsize_mb))
}

data_dir <- file.path("data", "bcc")
counts_file <- file.path(data_dir, "GSE123813_bcc_scRNA_counts.txt.gz")
meta_file <- file.path(data_dir, "GSE123813_bcc_tcell_metadata.txt.gz")
clin_file <- file.path(data_dir, "41591_2019_522_MOESM2_ESM.xlsx")
out_dir <- file.path("results", "bcc", "original_umap")
out_rds <- file.path(out_dir, "Yost_BCC_Tcells_with_author_UMAP.rds")
out_pdf <- file.path(out_dir, "Yost_BCC_original_author_UMAP.pdf")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

required_files <- c(counts_file, meta_file, clin_file)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required input file(s): ", paste(missing_files, collapse = ", "))
}

# 2. LOAD COUNTS 
counts <- fread(counts_file)

# Fix the column shift (Gene names in V1)
gene_col <- counts[[1]]
counts_mat <- as.matrix(counts[, -1])
rownames(counts_mat) <- gene_col

# Clean up memory
rm(counts)
gc()

# 3. LOAD T-CELL METADATA & CLINICAL RESPONSE
meta_tcell <- fread(meta_file)

# B. Load Clinical Data (Supplementary Table 1)
supp <- read_excel(clin_file, skip = 2)

# Clean Clinical Data
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

# C. Merge Clinical Data into T-Cell Metadata
meta_tcell_final <- meta_tcell %>%
  left_join(clin, by = "patient") %>%
  as.data.frame() # Convert to DF for Seurat

# Set row names for Seurat
rownames(meta_tcell_final) <- meta_tcell_final$cell.id

# 4. SUBSET COUNTS & CREATE SEURAT OBJECT
# Align the counts to the T-cell metadata (subsetting down from 53k to 33k cells)
common_cells <- intersect(colnames(counts_mat), rownames(meta_tcell_final))
counts_tcell_only <- counts_mat[, common_cells]
meta_tcell_final <- meta_tcell_final[common_cells, ]

# Create Object
obj <- CreateSeuratObject(
  counts = counts_tcell_only,
  meta.data = meta_tcell_final,
  project = "Yost_BCC_Tcells"
)


print(obj)
print(table(obj$response_binary, useNA="ifany"))
print(table(obj$cluster)) 

# Save
library(Seurat)

# 1. Check if the UMAP columns are in your metadata
# (They should be there from the previous 'merge' step)
head(obj@meta.data[, c("UMAP1", "UMAP2")])

# 2. Extract the coordinates into a matrix
# Seurat requires a matrix with 2 columns
author_umap <- as.matrix(obj@meta.data[, c("UMAP1", "UMAP2")])

# 3. Rename columns to what Seurat expects
colnames(author_umap) <- c("UMAP_1", "UMAP_2")

# 4. Create the Dimension Reduction Object
# This puts the author's coordinates into the "slot" where Seurat looks for UMAPs
obj[["umap"]] <- CreateDimReducObject(
  embeddings = author_umap,
  key = "UMAP_",
  assay = "RNA"
)

# Plot
p_author_umap <- DimPlot(obj, group.by = "cluster", label = TRUE) + 
  ggtitle("Yost BCC: Original Author UMAP")

ggsave(out_pdf, p_author_umap, width = 9, height = 6)
saveRDS(obj, out_rds)
