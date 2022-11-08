library(glmpca)
library(ggplot2)

source("src/get_data.R")
rm(d, info, zeros)
count <- assay(se)
metadata <- as.data.frame(colData(se))
colnames(count) <- metadata$UniqueCell_ID
rownames(count) <- rowData(se)$geneID

# GLM-PCA

# Negative Binomial:
res <- glmpca(count, 2, fam = "nb")
saveRDS(res, "data/processed_data/glmpca_nb.rds")

factors <- res$factors
colnames(factors) <- c("PC_1", "PC_2")
ggplot(factors, aes(x = PC_1, y = PC_2, color = metadata$sampleType)) +
  geom_point() +
  scale_color_discrete(name = "Sample Type")
ggplot(factors, aes(x = PC_1, y = PC_2, color = metadata$Patient)) +
  geom_point() +
  scale_color_discrete(name = "Patient")

# Poisson:
res2 <- glmpca(count, 2, fam = "poi")
saveRDS(res2, "data/processed_data/glmpca_poi.rds")

factors2 <- res2$factors
colnames(factors2) <- c("PC_1", "PC_2")
ggplot(factors2, aes(x = PC_1, y = PC_2, color = metadata$sampleType)) +
  geom_point() +
  scale_color_discrete(name = "Sample Type")
ggplot(factors2, aes(x = PC_1, y = PC_2, color = metadata$Patient)) +
  geom_point() +
  scale_color_discrete(name = "Patient")

