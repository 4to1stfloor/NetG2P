
data.folder <- "/mnt/gluster_server/data/network"
setwd(data.folder)


library(Rtsne)
library(ggplot2)
library(gridExtra)
# library(ggrepel)


cancer.type.list <- c("KIDNEY", "COADREAD", "UCEC", "BRCA", "LGG", "LUSC", "OV", "LUAD", "LIHC", "STAD", "BLCA", "CESC")

pca.data <- readRDS("total_cancer_pat_for_PCA.rds")
pca.data.meta <- readRDS("total_cancer_pat_for_PCA_meta.rds")
# pca.meta <- data.frame(pat = pca.data.meta$pat, cancertype = pca.data.meta$cancertype.x)

### all 
tsne_out <- Rtsne(as.matrix(t(pca.data)))
tsne_plot <- data.frame(TSNE1 = tsne_out$Y[,1], TSNE2= tsne_out$Y[,2])
tsne_plot$cancertype <- pca.data.meta$cancertype.x
tsne_plot$status <- pca.data.meta$vitalstatus
tsne_plot$prognosis <- pca.data.meta$cluster
g1<-ggplot(tsne_plot) + geom_point(aes(x=TSNE1, y=TSNE2, color = cancertype))
g2<-ggplot(tsne_plot) + geom_point(aes(x=TSNE1, y=TSNE2, color = status))
g3<-ggplot(tsne_plot) + geom_point(aes(x=TSNE1, y=TSNE2, color = prognosis))
grid.arrange(g1, g2, g3, ncol=2)

# ggplot(tsne_plot) + geom_point(aes(x=TSNE1, y=TSNE2, color = cancertype)) + geom_label_repel(aes(x=TSNE1, y=TSNE2,label = cancertype, fill = cancertype),
#                                                                                              color="white", fontface="bold",
#                                                                                              segment.color="grey30")




###### COADREAD, LUSC excluded ############### 
cancer.type.list <- c("KIDNEY", "UCEC", "BRCA", "LGG", "OV", "LUAD", "LIHC", "STAD", "BLCA", "CESC")

pca.data <- readRDS("total_cancer_pat_for_PCA.rds")
pca.data.meta <- readRDS("total_cancer_pat_for_PCA_meta.rds")
# pca.meta <- data.frame(pat = pca.data.meta$pat, cancertype = pca.data.meta$cancertype.x)

temp.pat.list <- dplyr::filter(pca.data.meta, 
                               cancertype.x %in% sapply(cancer.type.list, function(x) paste0("TCGA-",x)))$pat
temp.data <- pca.data[,temp.pat.list]
temp.meta.data <- dplyr::filter(pca.data.meta, pat %in% temp.pat.list)
sub.data<-temp.data[rowSums(temp.data) > 0,]


### 
tsne_out <- Rtsne(as.matrix(t(sub.data)))
tsne_plot <- data.frame(TSNE1 = tsne_out$Y[,1], TSNE2= tsne_out$Y[,2])
tsne_plot$cancertype <- temp.meta.data$cancertype.x
tsne_plot$status <- temp.meta.data$vitalstatus
tsne_plot$prognosis <- temp.meta.data$cluster
g1<-ggplot(tsne_plot) + geom_point(aes(x=TSNE1, y=TSNE2, color = cancertype))
g2<-ggplot(tsne_plot) + geom_point(aes(x=TSNE1, y=TSNE2, color = status))
g3<-ggplot(tsne_plot) + geom_point(aes(x=TSNE1, y=TSNE2, color = prognosis))
grid.arrange(g1, g2, g3, ncol=2)







####################  each cancer type
icancer <- 12

temp.pat.list <- dplyr::filter(pca.data.meta, cancertype.x == paste0("TCGA-",cancer.type.list[icancer]))$pat
temp.data <- pca.data[,temp.pat.list]
temp.meta.data <- dplyr::filter(pca.data.meta, pat %in% temp.pat.list)
sub.data<-temp.data[rowSums(temp.data) > 0,]

tsne_out <- Rtsne(as.matrix(t(sub.data)))
tsne_plot <- data.frame(TSNE1 = tsne_out$Y[,1], TSNE2= tsne_out$Y[,2])

tsne_plot$status <- temp.meta.data$vitalstatus
tsne_plot$prognosis <- temp.meta.data$cluster

g1<-ggplot(tsne_plot) + geom_point(aes(x=TSNE1, y=TSNE2, color = status)) + ggtitle(cancer.type.list[icancer])
g2<-ggplot(tsne_plot) + geom_point(aes(x=TSNE1, y=TSNE2, color = prognosis))+ ggtitle(cancer.type.list[icancer])

grid.arrange(g1, g2, ncol=2)




