library(glmpca)
library(ggplot2)
library(NatParksPalettes)
library(gridExtra)

source("src/get_data.R")
count <- as.matrix(assay(se))
count <- as.matrix(count)
metadata <- as.data.frame(colData(se))
colnames(count) <- metadata$UniqueCell_ID
rownames(count) <- rowData(se)$geneID

metadata$Tissue <- sapply(metadata$sampleType, function(x){
  if (x %in% c("NTC", "NTH", "NTR")) {
    return("Adjacent Normal Lung Tissues")
  } else if (x %in% c("PTC", "PTH", "PTR", "PTY")) {
    return("Peripheral Blood")
  } else{
    return("Tumor")
  }
})

# GLM-PCA

# Negative Binomial:
res_pc30 <- glmpca(count, 30, fam = "nb")
saveRDS(res, "data/processed_data/glmpca_nb_30pcs.rds")

factors <- res_pc30$factors
metadata["PC_1"] = factors[,1]
metadata["PC_2"] = factors[,2]

# color by sample type and tissue
p1 = ggplot(aes(x=PC_1,y=PC_2,color=sampleType),data=metadata) + 
  geom_point(alpha=0.7) +
  labs(subtitle = "Color by Sample Type") + 
  scale_color_brewer(palette="Paired") + 
  theme_bw()
p2 = ggplot(aes(x=PC_1,y=PC_2,color=Tissue),data=metadata) + 
  geom_point(alpha=0.7) +
  labs(subtitle = "Color by Tissue") +
  scale_color_brewer(palette="Paired") + 
  theme_bw() +
  theme(legend.position = "bottom")

grid.arrange(p1, p2, ncol = 2)

# color by patient
ggplot(aes(x=PC_1,y=PC_2,color=Patient),data=metadata) + 
  geom_point(alpha=0.7) +
  labs(subtitle = "Color by Patient") +
  theme_bw() 

# color by major cluster
ggplot(aes(x=PC_1,y=PC_2,color=majorCluster),data=metadata) + 
  geom_point(alpha=0.7) +
  labs(subtitle = "Color by Major Cluster") +
  theme_bw()


# Poisson:
res2_pc30 <- glmpca(count, 30, fam = "poi")
saveRDS(res2_pc30, "data/processed_data/glmpca_poi_30pcs.rds")

metadata2 = subset(metadata, select = c(UniqueCell_ID, Patient, majorCluster,
                                        sampleType, Tissue))
factors2 <- res2_pc30$factors
metadata2["PC_1"] = factors2[,1]
metadata2["PC_2"] = factors2[,2]

# color by sample type and tissue
p1 = ggplot(aes(x=PC_1,y=PC_2,color=sampleType),data=metadata2) + 
  geom_point(alpha=0.7) +
  labs(subtitle = "Color by Sample Type") + 
  scale_color_brewer(palette="Paired") + 
  theme_bw()
p2 = ggplot(aes(x=PC_1,y=PC_2,color=Tissue),data=metadata2) + 
  geom_point(alpha=0.7) +
  labs(subtitle = "Color by Tissue") +
  scale_color_brewer(palette="Paired") + 
  theme_bw() +
  theme(legend.position = "bottom")

grid.arrange(p1, p2, ncol = 2)

# color by patient
ggplot(aes(x=PC_1,y=PC_2,color=Patient),data=metadata2) + 
  geom_point(alpha=0.7) +
  labs(subtitle = "Color by Patient") +
  theme_bw() 

# color by major cluster
ggplot(aes(x=PC_1,y=PC_2,color=majorCluster),data=metadata2) + 
  geom_point(alpha=0.7) +
  labs(subtitle = "Color by Major Cluster") +
  theme_bw()
