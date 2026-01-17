library(Seurat)
library(tidyverse)
library(FNN)
library(ggalluvial)
library(cowplot)
library(gridExtra)
library(grid)
library(readr)

CD4_Obj <- readRDS("CD4.rds")
liu_counts <- readRDS("GSE179994_all.Tcell.rawCounts.rds/GSE179994_all.Tcell.rawCounts.rds")
liu_meta   <- read_tsv("GSE179994_Tcell.metadata.tsv/GSE179994_Tcell.metadata.tsv")
liu_meta_cd4 <- liu_meta %>% filter(celltype == "CD4")
counts_cd4   <- liu_counts[, liu_meta_cd4$cellid]

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


# 2. Create and Normalize Liu CD4 Seurat Object
liu_cd4_seurat <- CreateSeuratObject(counts_cd4, meta.data = as.data.frame(liu_meta_cd4)) %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%  # optional; keeps Seurat happy but not used
  ScaleData()
# 3. Preprocess Reference CD4 Seurat Object
CD4_Obj <- CD4_Obj %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  RunUMAP(reduction = "pca", dims = 1:20, return.model = TRUE)

# 4. Add missing genes (with zeros) to Liu object to match reference
all_genes <- rownames(CD4_Obj)
liu_cd4_seurat <- AddMissingGenes(liu_cd4_seurat, all_genes)

# 5. Find Anchors with All Genes
anchors <- FindTransferAnchors(
  reference = CD4_Obj,
  query = liu_cd4_seurat,
  dims = 1:20,
  features = all_genes,
  reference.reduction = "pca"
)

liu_cd4_seurat <- MapQuery(
  anchorset = anchors,
  query = liu_cd4_seurat,
  reference = CD4_Obj,
  refdata = list(predicted.celltype = "cell.type"),
  reference.reduction = "pca",
  reduction.model = "umap"
)

# 5. Majority Vote Labeling
ref_pca     <- Embeddings(CD4_Obj, reduction = "pca")[, 1:20]
query_pca   <- Embeddings(liu_cd4_seurat, reduction = "ref.pca")[, 1:20]
neighbors   <- get.knnx(data = ref_pca, query = query_pca, k = 10)$nn.index
ref_labels  <- CD4_Obj$cell.type
majority_vote <- apply(neighbors, 1, function(idx) {
  votes <- ref_labels[idx]
  tab <- sort(table(votes), decreasing = TRUE)
  if (length(tab) == 0 || max(tab) < 5) return(NA)
  names(tab)[1]
})
liu_cd4_seurat$majority_vote_pca <- majority_vote

# 6. Metadata and Label Mapping
liu_meta <- liu_cd4_seurat@meta.data %>%
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

#bar plots

library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(Seurat)
library(readxl)
response_info <- read_excel("response_info.xlsx")
# 1. Clean sample names in Seurat object
liu_cd4_seurat$sample_clean <- liu_cd4_seurat$sample %>%
  str_replace("^P0*", "P") %>%
  str_replace_all("\\.pre\\.0*", ".pre.") %>%
  str_replace_all("\\.post\\.0*", ".post.")

# 2. Clean sample names in response info
response_info_clean <- response_info %>%
  mutate(
    sample_clean = `Sample Name` %>%
      str_replace("^P0*", "P") %>%
      str_replace_all("\\.pre\\.0*", ".pre.") %>%
      str_replace_all("\\.post\\.0*", ".post.")
  )

# 3. Join Response info by sample_clean
liu_cd4_seurat@meta.data <- liu_cd4_seurat@meta.data %>%
  left_join(response_info_clean %>% select(sample_clean, Response), by = "sample_clean")

# 4. Convert response values to labels
liu_cd4_seurat$Response_Status <- case_when(
  liu_cd4_seurat$Response == "Yes" ~ "Responder",
  liu_cd4_seurat$Response == "No"  ~ "Non-responder",
  TRUE ~ NA_character_
)

# 5. Handle NA cluster labels for predicted clusters
liu_cd4_seurat$majority_vote_pca <- as.character(liu_cd4_seurat$majority_vote_pca)
liu_cd4_seurat$majority_vote_pca[is.na(liu_cd4_seurat$majority_vote_pca)] <- "NA"

# 6. CD4 Reference Cluster Counts (with NA)
cd4_cluster_counts <- liu_cd4_seurat@meta.data %>%
  count(majority_vote_pca, Response_Status, name = "n") %>%
  arrange(majority_vote_pca)

cd4_cluster_counts$majority_vote_pca <- factor(cd4_cluster_counts$majority_vote_pca,
                                               levels = unique(cd4_cluster_counts$majority_vote_pca))

# 7. Liu Cluster Counts
liu_cd4_cluster_counts <- liu_cd4_seurat@meta.data %>%
  count(cluster, Response_Status, name = "n") %>%
  arrange(cluster)

liu_cd4_cluster_counts$cluster <- factor(liu_cd4_cluster_counts$cluster,
                                         levels = unique(liu_cd4_cluster_counts$cluster))

#FIG A
# Plot: Liu Original Clusters
ggplot(liu_cd4_cluster_counts, aes(x = cluster, y = n, fill = Response_Status)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = n), position = position_dodge(width = 0.9), vjust = -0.3, size = 3.5) +
  labs(
    title = "Liu CD4 Cells by Original Cluster",
    x = "Liu Cluster",
    y = "Cell Count",
    fill = "Responder Status"
  ) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#FIG B
ggplot(cd4_cluster_counts, aes(x = majority_vote_pca, y = n, fill = Response_Status)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = n), position = position_dodge(width = 0.9), vjust = -0.3, size = 3.5) +
  labs(
    title = "Liu CD4 Cells by Predicted CD4 Reference Cluster (incl. NA)",
    x = "Predicted CD4 Cluster (majority_vote_pca)",
    y = "Cell Count",
    fill = "Responder Status"
  ) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#FIG C
# 8. Alluvial Plot
liu_cd4_seurat$predicted_label_num <- as.character(label_map[liu_cd4_seurat$majority_vote_pca])
summary_df <- liu_cd4_seurat@meta.data %>%
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
    title = "Alluvial Plot: Liu Clusters → Predicted CD4 Clusters",
    x = NULL, y = "Number of Cells", fill = "Liu Cluster"
  ) +
  theme_minimal(base_size = 14)

# 9. Combine Plot and Legend
final_plot <- plot_grid(p_alluvial, legend_table, rel_widths = c(2.5, 1), nrow = 1)
print(final_plot)

#FIG D
library(dplyr)
library(tidyr)
library(pROC)
library(ggplot2)

# Assuming your CD4 mapped object is liu_cd4_seurat (or liu_cd4)
# RENAME this variable if your object is named differently (e.g., mapped_liu_cd4)
cd4_obj <- liu_cd4_seurat 


# 2. DATA PREPARATION: Create clean lookup tables for the sample level
#-------------------------------------------------
meta <- cd4_obj@meta.data %>%
  mutate(
    Liu_Cluster = as.character(cluster),
    TCellMap_Cluster = as.character(majority_vote_pca), # Use predicted cluster name
    Response_Status = as.character(Response_Status),
    SampleID = as.character(sample) # Using 'sample' as the sample identifier
  ) %>%
  filter(!is.na(Response_Status) & Response_Status != "NA")

# Create a lookup table for sample response status
sample_status <- meta %>%
  select(SampleID, Response_Status) %>%
  distinct() %>%
  mutate(Response_Bin = ifelse(Response_Status == "Responder", 1, 0))


# 3. MODELING: Replicating the original logic for each sample-level model
#-------------------------------------------------

### === Model 1: Liu Clusters (per Sample) - Use ALL Liu Clusters === ###
props_m1 <- meta %>%
  count(SampleID, Liu_Cluster) %>%
  group_by(SampleID) %>%
  mutate(prop = n / sum(n)) %>%
  pivot_wider(names_from = Liu_Cluster, values_from = prop, values_fill = 0)

df_m1 <- props_m1 %>%
  left_join(sample_status, by = "SampleID")

# Build formula using all Liu clusters dynamically
features_m1 <- setdiff(names(df_m1), c("SampleID", "Response_Status", "Response_Bin", "n"))
formula_m1 <- as.formula(paste("Response_Bin ~", paste0("`", features_m1, "`", collapse = " + ")))

model1 <- glm(formula_m1, data = df_m1, family = binomial)
roc_m1 <- roc(df_m1$Response_Bin, predict(model1, type = "response"), quiet = TRUE)


### === Model 2: TCellMap Clusters (per Sample) - Use ALL TCellMap Clusters === ###
props_m2 <- meta %>%
  count(SampleID, TCellMap_Cluster) %>%
  group_by(SampleID) %>%
  mutate(prop = n / sum(n)) %>%
  pivot_wider(names_from = TCellMap_Cluster, values_from = prop, values_fill = 0)

df_m2 <- props_m2 %>%
  left_join(sample_status, by = "SampleID")

features_m2 <- setdiff(names(df_m2), c("SampleID", "Response_Status", "Response_Bin","n"))
formula_m2 <- as.formula(paste("Response_Bin ~", paste0("`", features_m2, "`", collapse = " + ")))

model2 <- glm(formula_m2, data = df_m2, family = binomial)
roc_m2 <- roc(df_m2$Response_Bin, predict(model2, type = "response"), quiet = TRUE)


### === Model 3: TCellMap in Treg/Regulatory-like Subset (per Sample, Filtered) === ###
# Filter based on TCellMap predicted Treg cluster: "CD4_c1_Treg"
props_m3 <- meta %>%
  filter(TCellMap_Cluster == "CD4_c1_Treg") %>%
  count(SampleID, Liu_Cluster) %>% # Use Liu clusters within the Treg TCellMap subset
  group_by(SampleID) %>%
  mutate(prop = n / sum(n)) %>%
  pivot_wider(names_from = Liu_Cluster, values_from = prop, values_fill = 0)

df_m3 <- props_m3 %>%
  left_join(sample_status, by = "SampleID")

features_m3 <- setdiff(names(df_m3), c("SampleID", "Response_Status", "Response_Bin","n"))
formula_m3 <- as.formula(paste("Response_Bin ~", paste0("`", features_m3, "`", collapse = " + ")))

model3 <- glm(formula_m3, data = df_m3, family = binomial)
roc_m3 <- roc(df_m3$Response_Bin, predict(model3, type = "response"), quiet = TRUE)


### === Model 4: Combined Clusters (per Sample) === ###
props_m4 <- meta %>%
  count(SampleID, TCellMap_Cluster, Liu_Cluster) %>%
  group_by(SampleID) %>%
  mutate(prop = n / sum(n)) %>%
  unite("combined_cluster", TCellMap_Cluster, Liu_Cluster, sep = "_") %>%
  pivot_wider(names_from = combined_cluster, values_from = prop, values_fill = 0)

df_m4 <- props_m4 %>%
  left_join(sample_status, by = "SampleID")

# Note: The 'n' column (subgroup count) is excluded to match the original analysis
features_m4 <- setdiff(names(df_m4), c("SampleID", "Response_Status", "Response_Bin", "n"))
formula_m4 <- as.formula(paste("Response_Bin ~", paste0("`", features_m4, "`", collapse = " + ")))

model4 <- glm(formula_m4, data = df_m4, family = binomial)
roc_m4 <- roc(df_m4$Response_Bin, predict(model4, type = "response"), quiet = TRUE)


# 4. VISUALIZATION: Plot all ROC curves on a single graph
#-------------------------------------------------
label_m1 = paste0("Model 1 (Liu Clusters): AUC = ", round(auc(roc_m1), 3))
label_m2 = paste0("Model 2 (TCellMap Clusters): AUC = ", round(auc(roc_m2), 3))
label_m3 = paste0("Model 3 (Treg-focused): AUC = ", round(auc(roc_m3), 3)) # Using Model 3 label for TCellMap-in-Treg
label_m4 = paste0("Model 4 (Combined): AUC = ", round(auc(roc_m4), 3))

# Combine the ROC curves for plotting (using all four models for a comprehensive view)
combined_roc_plot_cd4 <- ggroc(
  list(Model_1 = roc_m1, Model_2 = roc_m2, Model_3_Filtered = roc_m3, Model_4_Combined = roc_m4),
  size = 1.1
) +
  geom_abline(intercept = 1, slope = 1, linetype = "dashed", color = "grey40") +
  labs(
    title = "CD4 Predictive Models (Sample Level)",
    x = "Specificity (1 - False Positive Rate)",
    y = "Sensitivity (True Positive Rate)"
  ) +
  scale_color_viridis_d(
    name = "Model Performance",
    labels = c(label_m1, label_m2, label_m3, label_m4)
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold"))

# Print the plot
print(combined_roc_plot_cd4)

# 5. Save the plot as an SVG file
#-------------------------------------------------
ggsave(
  filename = "CD4_ROC_Combined_Models_Sample_Level.svg",
  plot = combined_roc_plot_cd4,
  width = 12,
  height = 9,
  units = "in"
)

cat("\nAll four CD4 sample-level logistic regression models have been fit, and the combined ROC curves plot has been saved as CD4_ROC_Combined_Models_Sample_Level.svg. 📈\n")
# 4. VISUALIZATION: Plot all ROC curves on a single graph (excluding Model 3)
#-------------------------------------------------
label_m1 = paste0("Model 1 (Liu Clusters): AUC = ", round(auc(roc_m1), 3))
label_m2 = paste0("Model 2 (TCellMap Clusters): AUC = ", round(auc(roc_m2), 3))
# Model 3 is excluded
label_m4 = paste0("Model 3 (Combined): AUC = ", round(auc(roc_m4), 3)) # Relabeled to Model 3 for the plot

# Combine the ROC curves for plotting (only Models 1, 2, and 4)
combined_roc_plot_cd4 <- ggroc(
  list(Model_1 = roc_m1, Model_2 = roc_m2, Model_3_Combined = roc_m4),
  size = 1.1
) +
  geom_abline(intercept = 1, slope = 1, linetype = "dashed", color = "grey40") +
  labs(
    title = "CD4 Predictive Models (Sample Level)",
    x = "Specificity (1 - False Positive Rate)",
    y = "Sensitivity (True Positive Rate)"
  ) +
  scale_color_viridis_d(
    name = "Model Performance",
    # Use only labels for Models 1, 2, and 4 (relabeled as 3)
    labels = c(label_m1, label_m2, label_m4)
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold"))

# Print the plot
print(combined_roc_plot_cd4)