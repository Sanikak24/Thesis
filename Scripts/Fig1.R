# Load required libraries
library(Seurat)
library(ggplot2)
library(tidyverse)
library(scales)
library(patchwork)
library(dplyr)
library(tidyr)
#INDIVIDUAL UMAPS
library(RColorBrewer)
#library(wesanderson)
library(viridis)
library(ggsci)
library(tidytext)
library(ggpubr)
library(cowplot)
# 1. Load data
CD8_Obj <- readRDS("CD8.rds")

#CD8 ATLAS
colorsForDataType <- c("#b5170f", "#f9844a", "#f9c74f", "#90be6d", "#43aa8b",
                       "#577590", "purple", "#4d9", "blue", "#f94144")

umapColor <- colorRampPalette(colorsForDataType)(length(unique(CD8_Obj$cell.type)))

# Assign identities
Idents(CD8_Obj) <- CD8_Obj$cell.type

# Extract UMAP coordinates
umap_df <- as.data.frame(Embeddings(CD8_Obj, reduction = "umap"))
umap_df$cell_type <- CD8_Obj$cell.type

# Compute centroid positions for each cell type
centroids <- umap_df %>%
  group_by(cell_type) %>%
  summarise(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2)) %>%
  mutate(cell_type_number = as.character(row_number()))

# Create a named vector for legend with numbers and names
legend_labels <- setNames(paste0(centroids$cell_type_number, ". ", centroids$cell_type), centroids$cell_type)

umap_plot <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = cell_type)) +
  geom_point(size = 0.4, alpha = 0.7) +  
  scale_color_manual(values = umapColor, labels = legend_labels) +
  theme_void() +
  theme(text = element_text(size = 14),
        legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 11),  # **Reduced legend text size from 13 to 11**
        legend.key.size = unit(1.8, "lines")) +  
  # Decreased centroid label size
  geom_text(data = centroids, aes(x = UMAP_1, y = UMAP_2, label = cell_type_number), 
            size = 3, color = "black", fontface = "bold") +  
  # Adjust legend point size
  guides(color = guide_legend(override.aes = list(size = 4)))


ggsave(file.path(figurePath, "UMAP_Figure_Optimized.pdf"), plot = umap_plot, 
       +        width = 160, height = 160, units = "mm", dpi = 3000, device = cairo_pdf)

#GENE EXPRESSION
GeneSet1 <- c("TCF7", "IL7R", "GPR183", "MGAT4A")
GeneSet2 <- c("GZMK")
GeneSet3 <- c("PDCD1", "CTLA4", "LAG3", "TIGIT")

gene_sets <- list(GeneSet1 = GeneSet1, GeneSet2 = GeneSet2, GeneSet3 = GeneSet3)
titles <- c("Sade-Feldman Response genes", 
            "Liu 2021 Precursor Exh\n(GZMK+)", 
            "Liu 2021 Coinhibitory\n(Terminal Exh)")

# Calculate Module Scores
# This creates columns: GeneSet1, GeneSet2, GeneSet3
CD8_Obj <- AddModuleScore(CD8_Obj, features = gene_sets, name = "GeneSet")

# 2. UMAPs: Gray -> Red (0 to 2.5)
# Any score < 0 (negative) will stay Gray.
# Any score > 2.5 will be Max Red.

my_palette <- colorRampPalette(c("gray", "red"))(50)
plot_list <- list()

for (i in 1:3) {
  col_name <- paste0("GeneSet", i) # Seurat's default naming
  
  p <- FeaturePlot(CD8_Obj, features = col_name, order = TRUE, raster = TRUE) +
    scale_color_gradientn(
      colors = my_palette,
      limits = c(0, 2.5),       # <--- FIXED LIMIT requested
      oob = scales::squish,     # Squish high values to max red, low values to gray
      name = "Module\nScore"
    ) +
    ggtitle(titles[i]) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
      legend.position = "right",
      legend.key.height = unit(0.5, "cm")
    )
  
  plot_list[[i]] <- p
}

# Combine and Save UMAPs
combined_umap <- wrap_plots(plot_list, ncol = 3) + plot_layout(guides = 'collect')
print(combined_umap)

ggsave("UMAP_gene_exp_Fig1.pdf", 
       plot = combined_umap, width = 12, height = 4.5)

# 3. Bar Plot: Average Module Score
# ==============================================================================
# Extract Data
meta_df <- FetchData(CD8_Obj, vars = c("cell.type", "GeneSet1", "GeneSet2", "GeneSet3"))

# Reshape
plot_data <- meta_df %>%
  pivot_longer(cols = starts_with("GeneSet"), names_to = "Module", values_to = "Score") %>%
  group_by(cell.type, Module) %>%
  summarise(Mean_Score = mean(Score), 
            SE = sd(Score) / sqrt(n()), 
            .groups = 'drop')

# Relabel for Legend
plot_data$Module <- factor(plot_data$Module, 
                           levels = c("GeneSet1", "GeneSet2", "GeneSet3"), 
                           labels = titles)

# Plot
p_bar <- ggplot(plot_data, aes(x = cell.type, y = Mean_Score, fill = Module)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = Mean_Score - SE, ymax = Mean_Score + SE),
                position = position_dodge(width = 0.8), width = 0.25) +
  
  # Colors
  scale_fill_manual(values = c("#4daf4a", "#e41a1c", "#377eb8")) + 
  
  # Styling
  labs(y = "Average Module Score", 
       x = "Clusters", 
       title = "Module Scores across Clusters") +
  
  # IMPORTANT: Add a dashed line at 0 so negative bars make sense
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(p_bar)

ggsave("BarPlot_ModuleScore.pdf", 
       plot = p_bar, width = 10, height = 6)