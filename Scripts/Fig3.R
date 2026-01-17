#FIG.A
#SAMPLE LEVEL
# 1. SETUP: Load all required libraries
#-------------------------------------------------
library(dplyr)
library(tidyr)
library(pROC)
library(ggplot2)

# 2. DATA PREPARATION: Create clean lookup tables for the sample level
# Create the main metadata table using 'sample_clean' as the identifier
meta <- liu_seurat@meta.data %>%
  mutate(
    Liu_Cluster = as.character(cluster),
    CD8_Cluster = as.character(majority_vote_pca),
    Response_Status = as.character(Response_Status),
    SampleID = as.character(sample_clean)
  ) %>%
  filter(!is.na(Response_Status) & Response_Status != "NA")

# Create a lookup table for sample response status
sample_status <- meta %>%
  select(SampleID, Response_Status) %>%
  distinct() %>%
  mutate(Response_Bin = ifelse(Response_Status == "Responder", 1, 0))


# 3. MODELING: Replicating the original logic for each sample-level model

### === Model 1: Liu Clusters (per Sample) === ###
props_m1 <- meta %>%
  count(SampleID, Liu_Cluster) %>%
  group_by(SampleID) %>%
  mutate(prop = n / sum(n)) %>%
  pivot_wider(names_from = Liu_Cluster, values_from = prop, values_fill = 0)

df_m1 <- props_m1 %>%
  left_join(sample_status, by = "SampleID")

model1 <- glm(Response_Bin ~ `Non-exhausted` + Prolif. + Tex, data = df_m1, family = binomial)
roc_m1 <- roc(df_m1$Response_Bin, predict(model1, type = "response"), quiet = TRUE)

### === Model 4: Combined Clusters (per Sample) === ###
props_m4 <- meta %>%
  count(SampleID, CD8_Cluster, Liu_Cluster) %>%
  group_by(SampleID) %>%
  mutate(prop = n / sum(n)) %>%
  unite("combined_cluster", CD8_Cluster, Liu_Cluster, sep = "_") %>%
  pivot_wider(names_from = combined_cluster, values_from = prop, values_fill = 0)

df_m4 <- props_m4 %>%
  left_join(sample_status, by = "SampleID")

# Note: The 'n' column (subgroup count) is excluded to match the original analysis
features_m4 <- setdiff(names(df_m4), c("SampleID", "Response_Status", "Response_Bin","n"))
formula_m4 <- as.formula(paste("Response_Bin ~", paste0("`", features_m4, "`", collapse = " + ")))

model4 <- glm(formula_m4, data = df_m4, family = binomial)
roc_m4 <- roc(df_m4$Response_Bin, predict(model4, type = "response"), quiet = TRUE)


# 4. VISUALIZATION: Plot all ROC curves on a single graph
#-------------------------------------------------
label_m1 = paste0("Model 1 (Liu): AUC = ", round(auc(roc_m1), 3))
#label_m2 = paste0("Model 2 (TCellMap): AUC = ", round(auc(roc_m2), 3))
label_m4 = paste0("Model 3 (Combined): AUC = ", round(auc(roc_m4), 3))

ggroc(list(roc_m1, roc_m2, roc_m4), size = 1.1) +
  geom_abline(intercept = 1, slope = 1, linetype = "dashed", color = "grey40") +
  labs(
    title = "Predictive Models (Sample Level)",
    x = "Specificity (1 - False Positive Rate)",
    y = "Sensitivity (True Positive Rate)"
  ) +
  scale_color_viridis_d(
    name = "Model Performance",
    labels = c(label_m1, label_m4)
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold"))
