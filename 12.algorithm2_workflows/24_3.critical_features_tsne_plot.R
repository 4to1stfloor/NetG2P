
data.folder <- "/mnt/gluster_server/data/network"
setwd(data.folder)


library(Rtsne)
library(ggplot2)
library(gridExtra)
# library(ggrepel)
library(RColorBrewer)


###### COADREAD, LUSC excluded ############### 
cancer.type.list <- c("KIDNEY", "UCEC", "BRCA", "LGG", "OV", "LUAD", "LIHC", "STAD", "BLCA", "CESC")
# group_color <- brewer.pal(length(cancer.type.list), 'Set3')

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



        
###
color_map <- c(
  "TCGA-KIDNEY" = "#E64B35FF",
  "TCGA-UCEC" =  "#4DBBD5FF",
  "TCGA-BRCA" = "#00A087FF",
  "TCGA-LGG" = "#3C5488FF",
  "TCGA-OV" = "#F39B7FFF",
  "TCGA-LUAD" = "#8491B4FF",
  "TCGA-LIHC" = "#91D1C2FF",
  "TCGA-STAD" = "#DC0000FF",
  "TCGA-BLCA" = "#7E6148FF",
  "TCGA-CESC" = "#B09C85FF"
)

setwd("~/nas/04.Results/bestfeatures/feature_network/")
png(filename = "Pan-cancer_tsne_plot2.png",
    width = 15, height = 15,  units = "cm" ,pointsize = 12,
    bg = "white", res = 1200, family = "")

tsne_total =ggplot(tsne_plot) +
  geom_point(aes(x=TSNE1, y=TSNE2, color = cancertype)) +
  scale_color_manual(values = color_map)+ 
  theme_bw()

print(tsne_total)
dev.off()


# g1<-ggplot(tsne_plot) + geom_point(aes(x=TSNE1, y=TSNE2, color = cancertype))
# g2<-ggplot(tsne_plot) + geom_point(aes(x=TSNE1, y=TSNE2, color = status))
# g3<-ggplot(tsne_plot) + geom_point(aes(x=TSNE1, y=TSNE2, color = prognosis))
# grid.arrange(g1, g2, g3, ncol=2)






####################  each cancer type
icancer <- 3

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




