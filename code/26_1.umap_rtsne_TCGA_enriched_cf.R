library(ggbiplot)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggfortify)
library(pheatmap)
library(RColorBrewer)
library(Rtsne)
# library(tsne)
# library(umap)
library(tidyr)
library(dplyr)
library(reshape2)
library(readxl)
library(ggridges)
library(ggdist)
library(ggpubr)
library(ggrepel)
library(e1071)
library(umap)

filepath = "/home/seokwon/nas/"
ref_path = paste0(filepath, "99.reference/")

## Read Dataset 
# meta.tbl <- read.csv(file = paste0(ref_path,'DepMap/Model.csv'), 
#                      sep=',', header=T, fill=T)
# 
# ge.tbl <- read.csv(file = paste0(ref_path, 'DepMap/CRISPRGeneEffect.csv'), 
#                    sep=',', header=T, fill=T, row.names=1)
# 
# meta_for_TCGA = read.csv(paste0(ref_path, "/Depmap_meta_filt_to_TCGA.csv"))

Cancerlist = dir(paste0(filepath, "/00.data/filtered_TCGA/"))
setwd("~/nas/04.Results/drug/depmap/")
####
Cancerlist = Cancerlist[c(-11,-12)]
# num_CancerType =  "04.TCGA-CESC"

for (num_CancerType in Cancerlist) {
  
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  Cancername = gsub('TCGA-' , '', CancerType)
  
  # call input
  # gc_cellline = readRDS(paste0("~/nas/00.data/filtered_TCGA/", num_CancerType, "/",Cancername,"_cellline_dual_all_log.rds"))
  gc_TCGA = readRDS(paste0("~/nas/04.Results/short_long/", CancerType,"_critical_features_short_long.rds"))
  
  # short_cluster = gc_TCGA %>% filter(cluster == "short") %>% dplyr::select(-vitalstatus, -duration,-status, -cluster)
  # long_cluster = gc_TCGA %>% filter(cluster == "long") %>% dplyr::select(-vitalstatus, -duration,-status, -cluster)
  wo_info = gc_TCGA %>% dplyr::select(-vitalstatus, -duration, -status)
  
  # gc_cellline$cluster = "cell"
  # gc_cellline_filt_df = gc_cellline[,colnames(wo_info)]
  # 
  # gc_cellline_fixed_df = readRDS(paste0("~/nas/04.Results/drug/depmap/gdsc/",CancerType, "_DM_sl_cluster.rds"))
  # info_predict = read.csv(paste0("~/nas/04.Results/cell/",Cancername,"_cellline_predict_for_bestmodel.csv"))
  # ### Dimension reduction (PCA)
  # # info_predict$adjust_predic
  # if (all.equal(rownames(gc_cellline_fixed_df) , info_predict$cellline_id)) {
  #   gc_cellline_fixed_df$predict = info_predict$adjust_predic
  # }
  
  set.seed(1)
  
  # combine_df = rbind(wo_info, gc_cellline_filt_df)
  test = umap(wo_info %>% select(-cluster))
  test_rtsne = Rtsne::Rtsne(wo_info %>% select(-cluster), check_duplicates = F)
  
  test_layout = as.data.frame(test$layout)
  colnames(test_layout) = c("PC1","PC2")
  test_layout$cluster = wo_info$cluster
  test_layout$vital_status = gc_TCGA$vitalstatus
  
  test_rtsne_layout = as.data.frame(test_rtsne$Y)
  colnames(test_rtsne_layout) = c("PC1","PC2")
  rownames(test_rtsne_layout) = rownames(wo_info)
  test_rtsne_layout$cluster = wo_info$cluster
  test_rtsne_layout$vital_status = gc_TCGA$vitalstatus

  total_cell_for_fig = gc_TCGA
  total_cell_meta = data.frame(cluster = total_cell_for_fig[,c("cluster")])
  total_cell_for_tmp = total_cell_for_fig %>% select(-cluster, -status, -duration,-vitalstatus)
  
  # # Convert the matrix to a numeric matrix
  
  tmp_numeric <- matrix(as.numeric(unlist(total_cell_for_tmp)), nrow = nrow(total_cell_for_tmp))
  vec <- as.vector(tmp_numeric)
  # hist(vec)
  total_cell_for_tmp[total_cell_for_tmp > quantile(vec, probs = 0.75) + 1.5*IQR(vec)] <- round(quantile(vec, probs = 0.75) + 1.5*IQR(vec),1)
  
  total_cell_for_fig_final = cbind(total_cell_for_tmp,total_cell_meta)

  tmp_cell = total_cell_for_fig_final %>% select(-cluster)
  min_max_normalize <- function(x) {
    (x - min(x)) / (max(x) - min(x)) * (1 - min(x)) + min(x)
  }

  non_zero = tmp_cell[tmp_cell!=0]
  non_zero = min_max_normalize(non_zero)
  tmp_cell_nor = tmp_cell
  tmp_cell_nor[tmp_cell != 0] = non_zero

  test_layout = cbind(test_layout,tmp_cell_nor)
  test_rtsne_layout = cbind(test_rtsne_layout,tmp_cell_nor)
  
  fig_path = paste0(filepath,"/04.Results/PCA_tsne/", Cancername)
  if(!dir.exists(fig_path)){
    dir.create(fig_path)
    print(paste0("Created folder: ", fig_path))
  } else {
    print(paste0("Folder already exists: ", fig_path))
  }
  setwd(fig_path)
  
  p_umap_sl = ggplot(test_layout , aes(x = PC1, 
                           y = PC2, 
                           color = cluster))+
    geom_point() +
    scale_color_manual(values = c("long" = "#4DAF4A" ,
                                  "short" = "#E41A1C")) +
    theme_classic() +
    xlab(label = "UMAP1")+
    ylab(label = "UMAP2")
  ggsave(file=paste0(CancerType , "_umap_sl.svg"), plot=p_umap_sl, width=6, height=6)
  
  p_rtsne_sl = ggplot(test_rtsne_layout , aes(x = PC1, 
                           y = PC2, 
                           color = cluster))+
    geom_point() +
    scale_color_manual(values = c("long" = "#4DAF4A" ,
                                  "short" = "#E41A1C")) +
    theme_classic() +
    xlab(label = "tSNE1")+
    ylab(label = "tSNE2")
  ggsave(file=paste0(CancerType , "_rtsne_sl.svg"), plot=p_rtsne_sl, width=6, height=6)
  
  p_umap_vital = ggplot(test_layout , aes(x = PC1, 
                           y = PC2, 
                           color = vital_status))+
    geom_point() +
    scale_color_manual(values = c("Alive" = "#4dbbd5" ,
                                  "Dead" = "#990066")) +
    theme_classic() +
    xlab(label = "UMAP1")+
    ylab(label = "UMAP2")
  ggsave(file=paste0(CancerType , "_umap_vital.svg"), plot=p_umap_vital, width=6, height=6)
  
  p_rtsne_vital = ggplot(test_rtsne_layout , aes(x = PC1, 
                           y = PC2, 
                           color = vital_status))+
    geom_point() +
    scale_color_manual(values = c("Alive" = "#4dbbd5" ,
                                  "Dead" = "#990066")) +
    theme_classic() +
    xlab(label = "tSNE1")+
    ylab(label = "tSNE2")
  ggsave(file=paste0(CancerType , "_rtsne_vital.svg"), plot=p_rtsne_vital, width=6, height=6)
  
  for (cf_u in colnames(test_rtsne_layout %>% select(-PC1,-PC2,-cluster, -vital_status))) {
    fic_layout = test_layout %>% select(PC1, PC2, any_of(cf_u))
    colnames(fic_layout) = c("PC1","PC2", "cf")
    
    p_umap = ggplot(fic_layout , 
                    aes(x = PC1, 
                        y = PC2,
                        color = cf))+
      geom_point() +
      scale_colour_gradient(low = "#F9F9C7", high = "red") +
      theme_classic() +
      xlab(label = "UMAP1")+
      ylab(label = "UMAP2")+
      labs(color= paste0(cf_u))+
      theme(panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_rect(fill = "white",
                                            colour = "white"))
    
    ggsave(file=paste0(CancerType ,"_",cf_u, "_umap.svg"), plot=p_umap, width=6, height=6)
  }

  
  for (cf_t in colnames(test_rtsne_layout %>% select(-PC1,-PC2,-cluster, -vital_status))) {
    fic_rtsne_layout = test_rtsne_layout %>% select(PC1, PC2, any_of(cf_t))
    colnames(fic_rtsne_layout) = c("PC1","PC2", "cf")
    
    p_rtsne = ggplot(fic_rtsne_layout , aes(x = PC1, 
                                   y = PC2, 
                                   color = cf))+
      geom_point() +
      scale_colour_gradient(low = "#F9F9C7", high = "red") +
      theme_classic() +
      xlab(label = "tSNE1")+
      ylab(label = "tSNE2")+
      labs(color= paste0(cf_t))+
      theme(panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_rect(fill = "white",
                                            colour = "white"))
    
    ggsave(file=paste0(CancerType ,"_",cf_t, "_rtsne.svg"), plot=p_rtsne, width=6, height=6)
  }

}


