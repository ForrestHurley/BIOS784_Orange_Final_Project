library(glmpca)
library(ggplot2)
library(NatParksPalettes)

source("src/get_data.R")
rm(d, info, zeros)
count <- as.matrix(assay(se))
count <- as.matrix(count)
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

metadata$Cell_type <- sapply(metadata$majorCluster, function(x){
  if (substr(x, 1, 3) == "CD4") {
    return("CD4")
  } else if (substr(x, 1, 3) == "CD8") {
    return("CD8")
  } else{
    return("Other")
  }
})
  
# GLM-PCA

# Negative Binomial:
#res <- glmpca(count, 2, fam = "nb")
#saveRDS(res, "data/processed_data/glmpca_nb.rds")

res_pc30 <- glmpca(count, 30, fam = "nb")
saveRDS(res, "data/processed_data/glmpca_nb_30pcs.rds")

factors <- res_pc30$factors
metadata["PC_1"] = factors[,1]
metadata["PC_2"] = factors[,2]
# metadata["PC_3"] = factors[,3]
# metadata["PC_3_cat"] = cut(factors[,3], breaks=c(-Inf, -10, -5, 0, 5, 10, Inf))

ggplot(aes(x=PC_1,y=PC_2,color=sampleType),data=metadata) + 
  geom_point(alpha=0.7) +
  #facet_wrap(~PC_3_cat) + 
  scale_color_brewer(palette="Paired") + 
  theme_bw()

ggplot(aes(x=PC_1,y=PC_2,color=Tissue),data=metadata) + 
  geom_point(alpha=0.7) +
  #facet_wrap(~PC_3_cat) + 
  scale_color_brewer(palette="Paired") + 
  theme_bw() +
  theme(legend.position = "bottom")

ggplot(aes(x=PC_1,y=PC_2,color=Patient),data=metadata) + 
  geom_point(alpha=0.7) +
  theme_bw() 

#ggplot(aes(x=PC_1,y=PC_2,color=Cell_type),data=metadata) + 
#  geom_point(alpha=0.7) +
  #facet_wrap(~PC_3_cat) + 
#  scale_color_brewer(palette="Paired", name = "Cell Type") + 
#  theme_bw()


# Poisson:
#res2 <- glmpca(count, 2, fam = "poi")
#saveRDS(res2, "data/processed_data/glmpca_poi.rds")
res2_pc30 <- glmpca(count, 30, fam = "poi")
saveRDS(res2_pc30, "data/processed_data/glmpca_poi_30pcs.rds")

metadata2 = subset(metadata, select = c(UniqueCell_ID, Patient, majorCluster,
                                        sampleType, Tissue))
factors2 <- res2_pc30$factors
metadata2["PC_1"] = factors2[,1]
metadata2["PC_2"] = factors2[,2]
# metadata2["PC_3"] = factors2[,3]
# metadata2["PC_3_cat"] = cut(factors2[,3], breaks=c(-Inf, -10, -5, 0, 5, 10, Inf))

ggplot(aes(x=PC_1,y=PC_2,color=sampleType),data=metadata2) + 
  geom_point(alpha=0.7) +
  #facet_wrap(~PC_3_cat) + 
  scale_color_brewer(palette="Paired") + 
  theme_bw()

ggplot(aes(x=PC_1,y=PC_2,color=Tissue),data=metadata2) + 
  geom_point(alpha=0.7) +
  #facet_wrap(~PC_3_cat) + 
  scale_color_brewer(palette="Paired") + 
  theme_bw() +
  theme(legend.position = "bottom")

ggplot(aes(x=PC_1,y=PC_2,color=Patient),data=metadata2) + 
  geom_point(alpha=0.7) +
  theme_bw() 

