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

# 8. Alluvial Plot
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

# 9. Combine Plot and Legend
final_plot <- plot_grid(p_alluvial, legend_table, rel_widths = c(2.5, 1), nrow = 1)
print(final_plot)

count_table <- table(
  Liu_Cluster = liu_seurat$cluster,
  Predicted_Cluster = liu_seurat$majority_vote_pca,
  useNA = "ifany"
)

print(count_table)

# Create table with NA values included
plot_df <- as.data.frame(table(liu_seurat$majority_vote_pca, useNA = "ifany"))

# Rename columns
colnames(plot_df) <- c("Predicted_Cluster", "Cell_Count")

# Convert NA to a string label for plotting
plot_df$Predicted_Cluster <- as.character(plot_df$Predicted_Cluster)
plot_df$Predicted_Cluster[is.na(plot_df$Predicted_Cluster)] <- "NA"

# Plot with NA and count labels
ggplot(plot_df, aes(x = reorder(Predicted_Cluster, -Cell_Count), y = Cell_Count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = Cell_Count), vjust = -0.5, size = 4) +
  labs(
    title = "Majority Vote Assignment (Including NA)",
    x = "Predicted CD8 Cluster",
    y = "Number of Liu CD8 Cells"
  ) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(0, max(plot_df$Cell_Count) * 1.1)


# 1. Prepare Neighbor Label Matrix
# Each row represents a query cell; each column is the label of one of its 10 neighbors
neighbor_labels_mat <- matrix(
  ref_labels[as.vector(neighbors)], 
  nrow = nrow(neighbors), 
  ncol = ncol(neighbors)
)

# 2. Identify the Target Label for each cell
# For assigned cells, it's the majority_vote_pca. 
# For NA cells, we find the "most-identified with" (top vote) even if it didn't hit your threshold.
top_vote_label <- apply(neighbor_labels_mat, 1, function(row) {
  tab <- sort(table(row), decreasing = TRUE)
  if (length(tab) == 0) return(NA)
  names(tab)[1]
})

# 3. Calculate Counts for the two categories
# 'Match_Count': neighbors matching the assigned/top label
# 'Other_Count': neighbors that belong to any other cell type
match_counts <- sapply(1:nrow(neighbor_labels_mat), function(i) {
  sum(neighbor_labels_mat[i, ] == top_vote_label[i], na.rm = TRUE)
})

other_counts <- ncol(neighbors) - match_counts

# 4. Create a Long-Format Dataframe for ggplot
# We keep track of the final classification to group the boxplots
plot_data <- data.frame(
  CellID = rownames(liu_seurat@meta.data),
  Assigned_Cluster = liu_seurat$majority_vote_pca,
  Match_Count = match_counts,
  Other_Count = other_counts
) %>%
  mutate(Assigned_Cluster = ifelse(is.na(Assigned_Cluster), "NA / Low Confidence", Assigned_Cluster)) %>%
  pivot_longer(cols = c(Match_Count, Other_Count), 
               names_to = "Neighbor_Type", 
               values_to = "Count") %>%
  mutate(Neighbor_Type = factor(Neighbor_Type, 
                                levels = c("Match_Count", "Other_Count"),
                                labels = c("Neighbors in Assigned/Top Type", "Neighbors in Other Types")))

# 5. Generate the Box Plot
boxplot<- ggplot(plot_data, aes(x = Assigned_Cluster, y = Count, fill = Neighbor_Type)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
  scale_fill_manual(values = c("Neighbors in Assigned/Top Type" = "#2ecc71", 
                               "Neighbors in Other Types" = "#e74c3c")) +
  labs(
    title = "Neighbor Distribution by Predicted Cluster",
    subtitle = "Comparing neighbors matching the predicted type vs. all other types",
    x = "Predicted CD8 Cluster (Reference)",
    y = "Number of Neighbors (k=10)",
    fill = "Neighbor Category"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

#ggsave(filename = "kNN_neighbor_distribution_boxplot.pdf",plot = boxplot,device = cairo_pdf,   # ensures true vector output
#width = 12,height = 8, units = "in")

library(tidyverse)
library(ggplot2)

# --- 1. CD8: Calculate Differences and Reorder ---
# Assuming 'plot_data' (CD8) is already created from your previous script

# Calculate the median difference for each cluster
cd8_rank <- plot_data %>%
  group_by(Assigned_Cluster) %>%
  summarize(
    match_med = median(Count[Neighbor_Type == "Neighbors in Assigned/Top Type"]),
    other_med = median(Count[Neighbor_Type == "Neighbors in Other Types"]),
    diff = match_med - other_med
  ) %>%
  arrange(desc(diff))  # Rank from highest difference to lowest

# Reorder the factor levels based on this rank
plot_data$Assigned_Cluster <- factor(plot_data$Assigned_Cluster, levels = cd8_rank$Assigned_Cluster)

# Redraw Fig 5A
p_box_cd8_ranked <- ggplot(plot_data, aes(x = Assigned_Cluster, y = Count, fill = Neighbor_Type)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
  scale_fill_manual(values = c("Neighbors in Assigned/Top Type" = "#2ecc71", 
                               "Neighbors in Other Types" = "#e74c3c")) +
  labs(
    x = "Predicted CD8 Cluster",
    y = "Number of Neighbors (k=10)",
    fill = "Category"
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")

library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(Seurat)
library(readxl)
response_info <- read_excel("response_info.xlsx")
# 1. Clean sample names in Seurat object
liu_seurat$sample_clean <- liu_seurat$sample %>%
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
liu_seurat@meta.data <- liu_seurat@meta.data %>%
  left_join(response_info_clean %>% select(sample_clean, Response), by = "sample_clean")

# 4. Convert response values to labels
liu_seurat$Response_Status <- case_when(
  liu_seurat$Response == "Yes" ~ "Responder",
  liu_seurat$Response == "No"  ~ "Non-responder",
  TRUE ~ NA_character_
)

# 5. Handle NA cluster labels for predicted clusters
liu_seurat$majority_vote_pca <- as.character(liu_seurat$majority_vote_pca)
liu_seurat$majority_vote_pca[is.na(liu_seurat$majority_vote_pca)] <- "NA"

# 6. CD8 Reference Cluster Counts (with NA)
cd8_cluster_counts <- liu_seurat@meta.data %>%
  count(majority_vote_pca, Response_Status, name = "n") %>%
  arrange(majority_vote_pca)

cd8_cluster_counts$majority_vote_pca <- factor(cd8_cluster_counts$majority_vote_pca,
                                               levels = unique(cd8_cluster_counts$majority_vote_pca))

# 7. Liu Cluster Counts
liu_cluster_counts <- liu_seurat@meta.data %>%
  count(cluster, Response_Status, name = "n") %>%
  arrange(cluster)

liu_cluster_counts$cluster <- factor(liu_cluster_counts$cluster,
                                     levels = unique(liu_cluster_counts$cluster))

# 8. Plot: CD8 Predicted Clusters (incl. NA)
ggplot(cd8_cluster_counts, aes(x = majority_vote_pca, y = n, fill = Response_Status)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = n), position = position_dodge(width = 0.9), vjust = -0.3, size = 3.5) +
  labs(
    title = "Liu CD8 Cells by Predicted CD8 Reference Cluster (incl. NA)",
    x = "Predicted CD8 Cluster (majority_vote_pca)",
    y = "Cell Count",
    fill = "Responder Status"
  ) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 9. Plot: Liu Original Clusters
ggplot(liu_cluster_counts, aes(x = cluster, y = n, fill = Response_Status)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = n), position = position_dodge(width = 0.9), vjust = -0.3, size = 3.5) +
  labs(
    title = "Liu CD8 Cells by Original Cluster",
    x = "Liu Cluster",
    y = "Cell Count",
    fill = "Responder Status"
  ) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 10. Wide format tables
cd8_cluster_wide <- cd8_cluster_counts %>%
  pivot_wider(names_from = Response_Status, values_from = n, values_fill = 0)

liu_cluster_wide <- liu_cluster_counts %>%
  pivot_wider(names_from = Response_Status, values_from = n, values_fill = 0)

# 11. Print count tables
cd8_cluster_wide
liu_cluster_wide

library(ggplot2)
library(dplyr)
library(scales)

# 1. Calculate Proportions
cd8_cluster_prop <- cd8_cluster_counts %>%
  group_by(majority_vote_pca) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()

# 2. Determine Rank Order based on Responder Proportion
# We sum the proportion where status is 'Responder' for each cluster
cluster_order <- cd8_cluster_prop %>%
  group_by(majority_vote_pca) %>%
  summarize(responder_prop = sum(proportion[Response_Status == "Responder"], na.rm = TRUE)) %>%
  arrange(desc(responder_prop)) %>% # Sort High to Low
  pull(majority_vote_pca)

# 3. Apply the New Order to the Factor
cd8_cluster_prop$majority_vote_pca <- factor(
  cd8_cluster_prop$majority_vote_pca,
  levels = cluster_order
)
#FIG B
ggplot(cd8_cluster_prop,
       aes(x = majority_vote_pca, y = proportion, fill = Response_Status)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(
    aes(label = n),
    position = position_stack(vjust = 0.5),
    size = 3.5,
    color = "white"
  ) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  labs(
    title = "Responder Status by Predicted CD8 Cluster (Relative Proportion)",
    x = "Predicted CD8 Cluster (majority_vote_pca)",
    y = "Relative Proportion of Cells",
    fill = "Responder Status"
  ) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(liu_cluster_prop,
       aes(x = cluster, y = proportion, fill = Response_Status)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(
    aes(label = n),
    position = position_stack(vjust = 0.5),
    size = 3.5,
    color = "white"
  ) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  labs(
    title = "Responder Status by Liu CD8 Cluster (Relative Proportion)",
    x = "Liu Cluster",
    y = "Relative Proportion of Cells",
    fill = "Responder Status"
  ) +
  theme_minimal(base_size = 14)


#CD4
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

ggsave(
  filename = "Liu_CD4_Cluster_Prediction_Alluvial.svg",
  plot = final_plot,
  width = 12, # Adjust width as needed for better visualization of both the plot and the legend
  height = 8, # Adjust height as needed
  units = "in"
)

count_table <- table(
  Liu_Cluster = liu_cd4_seurat$cluster,
  Predicted_Cluster = liu_cd4_seurat$majority_vote_pca,
  useNA = "ifany"
)

print(count_table)

# Create table with NA values included
plot_df <- as.data.frame(table(liu_cd4_seurat$majority_vote_pca, useNA = "ifany"))

# Rename columns
colnames(plot_df) <- c("Predicted_Cluster", "Cell_Count")

# Convert NA to a string label for plotting
plot_df$Predicted_Cluster <- as.character(plot_df$Predicted_Cluster)
plot_df$Predicted_Cluster[is.na(plot_df$Predicted_Cluster)] <- "NA"

# Plot with NA and count labels
majority_count<-ggplot(plot_df, aes(x = reorder(Predicted_Cluster, -Cell_Count), y = Cell_Count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = Cell_Count), vjust = -0.5, size = 4) +
  labs(
    title = "Majority Vote Assignment (Including NA)",
    x = "Predicted CD4 Cluster",
    y = "Number of Liu CD4 Cells"
  ) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(0, max(plot_df$Cell_Count) * 1.1)


ggsave(
  filename = "CD4_majority_vote_count.svg",
  plot = majority_count,
  width = 10,
  height = 8,
  units = "in"
)

cat("\nThe CD4 UMAP overlay plot, 'Liu CD4 overlaid on Reference CD4 UMAP', has been successfully saved as Liu_CD4_UMAP_Overlay_Projection.svg\n")

# 1. Map labels to the neighbor indices
# 'neighbors' and 'ref_labels' should already exist from your CD4 script
neighbor_labels_mat <- matrix(
  ref_labels[as.vector(neighbors)], 
  nrow = nrow(neighbors), 
  ncol = ncol(neighbors)
)

# 2. Identify the target label for the boxes
# We use the majority vote, but for NA cells, we identify the 'top' vote 
# among the neighbors to see how close it was to the threshold.
top_vote_label <- apply(neighbor_labels_mat, 1, function(row) {
  tab <- sort(table(row), decreasing = TRUE)
  if (length(tab) == 0) return(NA)
  names(tab)[1]
})

# 3. Calculate match vs. other counts
match_counts <- sapply(1:nrow(neighbor_labels_mat), function(i) {
  sum(neighbor_labels_mat[i, ] == top_vote_label[i], na.rm = TRUE)
})

other_counts <- 10 - match_counts

# 4. Prepare data for plotting
plot_data_cd4 <- data.frame(
  Assigned_Cluster = liu_cd4_seurat$majority_vote_pca,
  Match_Count = match_counts,
  Other_Count = other_counts
) %>%
  mutate(Assigned_Cluster = ifelse(is.na(Assigned_Cluster), "NA / Low Confidence", Assigned_Cluster)) %>%
  pivot_longer(cols = c(Match_Count, Other_Count), 
               names_to = "Neighbor_Category", 
               values_to = "Neighbor_Count") %>%
  mutate(Neighbor_Category = factor(Neighbor_Category, 
                                    levels = c("Match_Count", "Other_Count"),
                                    labels = c("Neighbors in Assigned/Top Type", "Neighbors in Other Types")))

# 5. Generate Boxplot
ggplot(plot_data_cd4, aes(x = Assigned_Cluster, y = Neighbor_Count, fill = Neighbor_Category)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.8) +
  scale_fill_manual(values = c("Neighbors in Assigned/Top Type" = "#2ecc71", 
                               "Neighbors in Other Types" = "#e74c3c")) +
  labs(
    title = "Neighbor Distribution: CD4 Assignment Confidence",
    subtitle = "Comparing neighbors matching the top-voted type vs. all others",
    x = "Predicted CD4 Cluster (Reference)",
    y = "Count of Neighbors (k=10)",
    fill = "Neighbor Category"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )
# 1. Create the plot object
p_confidence_cd4 <- ggplot(plot_data_cd4, aes(x = Assigned_Cluster, y = Neighbor_Count, fill = Neighbor_Category)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.8) +
  scale_fill_manual(values = c("Neighbors in Assigned/Top Type" = "#2ecc71", 
                               "Neighbors in Other Types" = "#e74c3c")) +
  labs(
    title = "Neighbor Distribution: CD4 Assignment Confidence",
    subtitle = "Comparing neighbors matching the top-voted type vs. all others",
    x = "Predicted CD4 Cluster (Reference)",
    y = "Count of Neighbors (k=10)",
    fill = "Neighbor Category"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    plot.background = element_rect(fill = "white", color = NA)
  )

# 2. Save as SVG (Best for editing in Illustrator/Inkscape)
ggsave(
  filename = "CD4_Boxplot.svg",
  plot = p_confidence_cd4,
  width = 10,
  height = 7,
  units = "in",
  device = "svg"
)

# 3. Save as PDF (Best for document embedding)
ggsave(
  filename = "CD4_Boxplot.pdf",
  plot = p_confidence_cd4,
  width = 10,
  height = 7,
  units = "in",
  device = "pdf"
)

cd4_rank <- plot_data_cd4 %>%
  group_by(Assigned_Cluster) %>%
  summarize(
    # Use 'Neighbor_Count' for values and 'Neighbor_Category' for labels
    match = median(Neighbor_Count[Neighbor_Category == "Neighbors in Assigned/Top Type"], na.rm = TRUE),
    other = median(Neighbor_Count[Neighbor_Category == "Neighbors in Other Types"], na.rm = TRUE),
    diff = match - other
  ) %>%
  arrange(desc(diff))

print(cd4_rank)

plot_data_cd4$Assigned_Cluster <- factor(plot_data_cd4$Assigned_Cluster, levels = cd4_rank$Assigned_Cluster)

#FIG C
# Redraw Fig 5B
p_box_cd4_ranked <- ggplot(plot_data_cd4, aes(x = Assigned_Cluster, y = Neighbor_Count, fill = Neighbor_Category)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
  scale_fill_manual(values = c("Neighbors in Assigned/Top Type" = "#2ecc71", 
                               "Neighbors in Other Types" = "#e74c3c")) +
  labs(
    x = "Predicted CD4 Cluster",
    y = "Count of Neighbors (k=10)",
    fill = "Category"
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")


# Save as Vector Files 
ggsave("Fig_5A_CD8_Box_Plot_Ranked.pdf", p_box_cd8_ranked, width = 12, height = 8)
ggsave("Fig_5B_CD4_Box_Plot_Ranked.pdf", p_box_cd4_ranked, width = 12, height = 8)

#FIG D
#FINAL PLOT
library(dplyr)
library(tidyr)
library(scales)
library(cowplot)

# ==============================================================================
# 1. Common Data Prep
# ==============================================================================

base_data <- cd4_cluster_counts %>%
  mutate(
    Response_Status = replace_na(as.character(Response_Status), "NA"),
    majority_vote_pca = replace_na(as.character(majority_vote_pca), "NA")
  ) %>%
  group_by(majority_vote_pca) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()

# Common Ranking: Rank by proportion of "Responder" cells
cluster_order <- base_data %>%
  group_by(majority_vote_pca) %>%
  summarize(responder_prop = sum(proportion[Response_Status == "Responder"], na.rm = TRUE)) %>%
  arrange(desc(responder_prop)) %>%
  pull(majority_vote_pca)

base_data$majority_vote_pca <- factor(base_data$majority_vote_pca, levels = cluster_order)

# Define Colors Globally
my_colors <- c(
  "Non-responder" = "#F8766D",
  "Responder"     = "#00C0C0",
  "NA"            = "#808080"
)

# ==============================================================================
# Function to create plot with specific stack order
# ==============================================================================

create_stack_plot <- function(data, level_order, title_suffix) {
  
  # ggplot stacks from bottom → top, so reverse the desired top→bottom order
  data$Response_Status <- factor(data$Response_Status, levels = rev(level_order))
  
  ggplot(data, aes(x = majority_vote_pca, y = proportion, fill = Response_Status)) +
    geom_bar(stat = "identity", position = "stack", width = 0.85) +
    geom_text(
      aes(label = ifelse(proportion > 0.02, n, "")),
      position = position_stack(vjust = 0.5),
      size = 3,
      color = "white"
    ) +
    scale_y_continuous(expand = c(0, 0), labels = scales::number_format(accuracy = 0.01)) +
    scale_fill_manual(values = my_colors) +
    labs(
      y = "Relative Proportion of Cells",
      x = "Predicted CD4 Cluster(majority_vote_pca)",
      fill = "Responder Status"
    ) +
    theme_classic(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank(),
      plot.margin = margin(10, 10, 10, 10)
    )
}

# ==============================================================================
# 2. Generate Plot 1: Non-responder, Responder, NA
# ==============================================================================

order_1 <- c("NA", "Responder", "Non-responder")
p1 <- create_stack_plot(base_data, order_1, "Non-resp / Resp / NA")
print(p1)

# ==============================================================================
# 3. Generate Plot 2: Non-responder, NA, Responder
# ==============================================================================

order_2 <- c("Responder", "NA", "Non-responder")
p2 <- create_stack_plot(base_data, order_2, "Non-resp / NA / Resp")
print(p2)
#ggsave("Fig5b_CD4_relative_proportion.svg", plot = p1, width = 10, height = 6)
#ggsave("Fig5b_CD4_relative_proportion_NA_middle.svg", plot = p2, width = 10, height = 6)

ggplot(liu_cd4_cluster_prop,
       aes(x = cluster, y = proportion, fill = Response_Status)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(
    aes(label = n),
    position = position_stack(vjust = 0.5),
    size = 3.5,
    color = "white"
  ) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  labs(
    title = "Responder Status by Liu CD4 Cluster (Relative Proportion)",
    x = "Liu Cluster",
    y = "Relative Proportion of Cells",
    fill = "Responder Status"
  ) +
  theme_minimal(base_size = 14)
