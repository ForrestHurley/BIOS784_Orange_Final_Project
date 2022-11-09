
library(SummarizedExperiment)
library(irlba)
library(NatParksPalettes)
library(ggplot2)
library(tsne)

se <- readRDS("data/processed_data/se.rds")
assay(se) <- as.matrix(assay(se)) # reduces memory usage by 1gb for rest of script??? IDK, but remove later

pca <- prcomp_irlba(log(1+t(assay(se))), n=30, scale=TRUE)

metadata <- as.data.frame(colData(se))

metadata["PCA1"] = pca$x[,1]
metadata["PCA2"] = pca$x[,2]
metadata["PCA3"] = pca$x[,3]
metadata["PCA3_cat"] = cut(pca$x[,3], breaks=c(-Inf, -10, -5, 0, 5, 10, Inf))

ggplot(aes(x=PCA1,y=PCA2,color=sampleType),data=metadata) + 
	geom_point(alpha=0.7) +
	facet_wrap(~PCA3_cat) + 
	scale_color_brewer(palette="Paired") + 
	theme_bw()

