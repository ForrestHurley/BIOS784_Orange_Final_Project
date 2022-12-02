## Use first 30 principal components for Gaussian Mixture Modeling to speed up computation
library(irlba)
library(mclust)
library(SummarizedExperiment)
library(tidyverse)
se <- readRDS("data/processed_data/se.rds")
express <- log2(t(assay(se)) + 1)
pc_express <- prcomp_irlba(express, n = 30, scale = T)

## 4 models
## G = 3, each tissue site
## G = 14 patients
## G = 12 Cell types
## G = 16 Publication
gmm <- Mclust(pc_express$x, G = c(3,12,14,16))

saveRDS(gmm, "data/processed_data/pc_gauss_mixture.rds")

se$gmm_cluster <- gmm$classification
gmm_bic <- gmm$BIC[,"VEV"] %>% as.data.frame()
colnames(gmm_bic) <- c("BIC")

png("figures/gmm_bic.png")
gmm_bic %>% rownames_to_column("Components") %>%
  ggplot(aes(x = as.numeric(Components), y = BIC)) +
  geom_point() +
  geom_line() +
  labs(x = "Components") +
  theme_bw()
dev.off()

## Overlay regular pca with GMM clusters
se$PC1 <- pc_express$x[,1]
se$PC2 <- pc_express$x[,2]

png("Figures/gmm_log2_30pc.png", width = 8.5, height = 6.8, units = "in",res = 200)
ggplot(as.data.frame(colData(se)), aes(x = PC1, y = PC2, color = as.factor(gmm_cluster))) +
  geom_point(alpha = 0.5) +
  stat_ellipse() +
  scale_color_discrete(name = "GMM Cluster") +
  theme_bw()
dev.off()

png("Figures/gmm_sampletype.png",width = 8.5, height = 6.8, units = "in",res = 200)
ggplot(as.data.frame(colData(se)), aes(x = PC1, y = PC2)) +
  geom_point(alpha = 0.5, aes(color = sampleType)) +
  stat_ellipse(aes(group = as.factor(gmm_cluster))) +
  scale_color_discrete(name = "Sample Type") +
  theme_bw()
dev.off()

png("Figures/gmm_patient.png",width = 8.5, height = 6.8, units = "in",res = 200 )
ggplot(as.data.frame(colData(se)), aes(x = PC1, y = PC2)) +
  geom_point(alpha = 0.5, aes(color = Patient)) +
  stat_ellipse(aes(group = as.factor(gmm_cluster))) +
  scale_color_discrete(name = "Patient") +
  theme_bw()
dev.off()

## Overlay GLM-PCA with GMM clusters
# glm_nb <- readRDS("data/processed_data/glmpca_nb_30pcs.rds")
# 
# nb <- glm_nb$factors %>%
#   rownames_to_column("UniqueCell_ID") %>%
#   inner_join(x = .,
#              y = as.data.frame(colData(se)),
#              by = c("UniqueCell_ID"))
# ggplot(nb, aes(x = -dim1, y = -dim2, color = as.factor(gmm_cluster))) +
#   geom_point(alpha = 0.7) +
#   # scale_color_brewer(palette = "Paired") +
#   theme_classic()
# 
# glm_poi <- readRDS("data/processed_data/glmpca_poi_30pcs.rds")
# 
# poi <- glm_poi$factors %>%
#   rownames_to_column("UniqueCell_ID") %>%
#   inner_join(x = .,
#              y = as.data.frame(colData(se)),
#              by = c("UniqueCell_ID"))
# ggplot(poi, aes(x = -dim1, y = -dim2, color = as.factor(gmm_cluster))) +
#   geom_point(alpha = 0.5) +
#   # scale_color_brewer(palette = "Paired") +
#   theme_classic()
