## Use first 30 principal components for Gaussian Mixture Modeling to speed up computation
library(irlba)
library(mclust)
library(SummarizedExperiment)
library(tidyverse)
se <- readRDS("data/processed_data/se.rds")
express <- log2(t(assay(se)) + 1)
pc_express <- prcomp_irlba(express, n = 30, scale = T)

gmm <- Mclust(pc_express$x, G = 12)

saveRDS(gmm, "data/processed_data/pc_gauss_mixture.rds")

se$gmm_cluster <- gmm$classification

## Overlay regular pca with GMM clusters
se$PC1 <- pc_express$x[,1]
se$PC2 <- pc_express$x[,2]
ggplot(as.data.frame(colData(se)), aes(x = PC1, y = PC2, color = as.factor(gmm_cluster))) +
  geom_point(alpha = 0.5) +
  scale_color_brewer(palette = "Paired") +
  theme_classic()

## Overlay GLM-PCA with GMM clusters
glm_nb <- readRDS("data/processed_data/glmpca_nb.rds")

nb <- glm_nb$factors %>%
  rownames_to_column("UniqueCell_ID") %>%
  inner_join(x = .,
             y = as.data.frame(colData(se)),
             by = c("UniqueCell_ID"))
ggplot(nb, aes(x = dim1, y = dim2, color = as.factor(gmm_cluster))) +
  geom_point(alpha = 0.7) +
  scale_color_brewer(palette = "Paired") +
  theme_classic()

glm_poi <- readRDS("data/processed_data/glmpca_poi.rds")

poi <- glm_poi$factors %>%
  rownames_to_column("UniqueCell_ID") %>%
  inner_join(x = .,
             y = as.data.frame(colData(se)),
             by = c("UniqueCell_ID"))
ggplot(poi, aes(x = dim1, y = dim2, color = as.factor(gmm_cluster))) +
  geom_point(alpha = 0.5) +
  scale_color_brewer(palette = "Paired") +
  theme_classic()