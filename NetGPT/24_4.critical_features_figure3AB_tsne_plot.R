
data.folder <- "~/nas/04.Results/PCA_tsne/"
setwd(data.folder)

library(Rtsne)
library(ggplot2)
library(gridExtra)
# library(ggrepel)
library(RColorBrewer)
# library(umap)
# library(scatterplot3d)
# library(plotly) 

# group_color <- brewer.pal(length(cancer.type.list), 'Set3')

pca.data <- readRDS("total_cancer_pat_for_PCA.rds")
pca.data.meta <- readRDS("total_cancer_pat_for_PCA_meta.rds")

# unique(pca.data.meta$cancertype.x)
cancer.type.list <- c("UCEC", "BRCA", "LGG","LUSC" ,"OV", "LUAD", "LIHC", "STAD", "BLCA", "CESC")
temp.pat.list <- dplyr::filter(pca.data.meta,
                               cancertype.x %in% sapply(cancer.type.list, function(x) paste0("TCGA-",x)))$pat
temp.data <- pca.data[,temp.pat.list]
sub.data<-temp.data[rowSums(temp.data) > 0,]

### 

tsne_out <- Rtsne(as.matrix(t(sub.data)))

tsne_plot <- data.frame(TSNE1 = tsne_out$Y[,1], TSNE2= tsne_out$Y[,2])
tsne_plot$cancertype <- pca.data.meta$cancertype.x
tsne_plot$status <- pca.data.meta$vitalstatus
tsne_plot$prognosis <- pca.data.meta$cluster

#####
color_map <- c(
  "TCGA-UCEC" =  "#E64B35FF",
  "TCGA-BRCA" = "#4DBBD5FF",
  "TCGA-LGG" = "#00A087FF",
  "TCGA-LUSC" = "#B09C85FF",
  "TCGA-OV" = "#3C5488FF",
  "TCGA-LUAD" = "#F39B7FFF",
  "TCGA-LIHC" = "#8491B4FF",
  "TCGA-STAD" = "#91D1C2FF",
  "TCGA-BLCA" = "#7E6148FF",
  "TCGA-CESC" = "#DC0000FF"
)

pat_cluster_map <- c(
  "short" =  "#E41A1C",
  "long" = "#4DAF4A"
)

total_tsne = ggplot(tsne_plot) +
  geom_point(aes(x=TSNE1, y=TSNE2, color = cancertype)) +
  scale_color_manual(values = color_map)+
  theme_bw()

ggsave(file="figure3A.svg", plot=total_tsne, width=10, height=10)

total_sl_tsne =ggplot(tsne_plot) +
  geom_point(aes(x=TSNE1, y=TSNE2, color = prognosis)) +
  scale_color_manual(values = pat_cluster_map)+ 
  theme_bw()

ggsave(file="figure3B.svg", plot=total_sl_tsne, width=10, height=10)

# plot_ly(umap_out_df, x = umap_out_df$V1, y = umap_out_df$V2, z= umap_out_df$V3, color = umap_out_df$cancertype , colors = color_map) 

# unmap_out = umap(as.matrix(t(pca.data)), n_components = 2, random_state = 15) 
# umap_out_df = as.data.frame(unmap_out$layout)
# umap_out_df$cancertype <- pca.data.meta$cancertype.x

# ggplot(umap_out_df) +
#   geom_point(aes(x=V1, y=V2, color = cancertype)) +
#   scale_color_manual(values = color_map)+ 
#   theme_bw()


# g1<-ggplot(tsne_plot) + geom_point(aes(x=TSNE1, y=TSNE2, color = cancertype))
# g2<-ggplot(tsne_plot) + geom_point(aes(x=TSNE1, y=TSNE2, color = status))
# g3<-ggplot(tsne_plot) + geom_point(aes(x=TSNE1, y=TSNE2, color = prognosis))
# grid.arrange(g1, g2, g3, ncol=2)






# ####################  each cancer type
# icancer <- 3
# 
# temp.pat.list <- dplyr::filter(pca.data.meta, cancertype.x == paste0("TCGA-",cancer.type.list[icancer]))$pat
# temp.data <- pca.data[,temp.pat.list]
# temp.meta.data <- dplyr::filter(pca.data.meta, pat %in% temp.pat.list)
# sub.data<-temp.data[rowSums(temp.data) > 0,]
# 
# tsne_out <- Rtsne(as.matrix(t(sub.data)))
# tsne_plot <- data.frame(TSNE1 = tsne_out$Y[,1], TSNE2= tsne_out$Y[,2])
# 
# tsne_plot$status <- temp.meta.data$vitalstatus
# tsne_plot$prognosis <- temp.meta.data$cluster
# 
# g1<-ggplot(tsne_plot) + geom_point(aes(x=TSNE1, y=TSNE2, color = status)) + ggtitle(cancer.type.list[icancer])
# g2<-ggplot(tsne_plot) + geom_point(aes(x=TSNE1, y=TSNE2, color = prognosis))+ ggtitle(cancer.type.list[icancer])
# 
# grid.arrange(g1, g2, ncol=2)




