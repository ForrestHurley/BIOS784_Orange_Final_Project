
library(SummarizedExperiment)
library(irlba)
library(NatParksPalettes)
library(ggplot2)

se <- readRDS("data/processed_data/se.rds")
assay(se) <- as.matrix(assay(se)) # reduces memory usage by 1gb for rest of script??? IDK, but remove later

row_vars.se <- sort(rowVars(log(1+assay(se))))
sum(tail(row_vars.se,30)) / sum(row_vars.se)

pca <- prcomp_irlba(log(1+t(assay(se))), n=10, scale=TRUE)
pca_unscaled <- prcomp_irlba(log(1+t(assay(se))), n=10, scale=FALSE)

metadata <- as.data.frame(colData(se))

metadata["PCA1"] = pca$x[,1]
metadata["PCA2"] = pca$x[,2]
metadata["PCA3"] = pca$x[,3]
metadata["PCA3_cat"] = cut(pca$x[,3], breaks=c(-Inf, -10, -5, 0, 5, 10, Inf))

ggplot(aes(x=PCA1,y=PCA2,color=sampleType),data=metadata) + 
	geom_point(alpha=0.4) +
	facet_wrap(~PCA3_cat) + 
	scale_color_brewer(palette="Paired") + 
	theme_bw()

pca_var = as.data.frame(pca$sdev)
colnames(pca_var) <- c("sdev")

ggplot(aes(x=1:length(sdev),y=sdev^2/pca$totalvar),data=pca_var) + 
	geom_line(color="red") + 
	geom_line(aes(y=cumsum(sdev^2)/pca$totalvar),color="blue") +
	labs(x="PC",y="Proportion of Variance") +
	scale_x_continuous(breaks=c(0,100,200,300,400,500,600,700,800,900,1000))


pca_var_unscaled = as.data.frame(pca_unscaled$sdev)
colnames(pca_var_unscaled) <- c("sdev")

ggplot(aes(x=1:length(sdev),y=sdev^2/pca_unscaled$totalvar),data=pca_var_unscaled) + 
	geom_line(color="red") + 
	geom_line(aes(y=cumsum(sdev^2)/pca_unscaled$totalvar),color="blue") +
	labs(x="PC",y="Proportion of Variance") +
	scale_x_continuous(breaks=c(0,10,20,30,40,50,60,70,80,90,100))

