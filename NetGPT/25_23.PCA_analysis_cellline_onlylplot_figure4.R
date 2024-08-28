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
filepath = "/home/seokwon/nas/"
ref_path = paste0(filepath, "99.reference/")

## Read Dataset 
meta.tbl <- read.csv(file = paste0(ref_path,'DepMap/Model.csv'), 
                     sep=',', header=T, fill=T)

ge.tbl <- read.csv(file = paste0(ref_path, 'DepMap/CRISPRGeneEffect.csv'), 
                   sep=',', header=T, fill=T, row.names=1)

meta_for_TCGA = read.csv(paste0(ref_path, "/Depmap_meta_filt_to_TCGA.csv"))

Cancerlist = dir(paste0(filepath, "/00.data/filtered_TCGA/"))

####
Cancerlist = Cancerlist[c(-11,-12)]
num_CancerType = "11.TCGA-STAD"

for (num_CancerType in Cancerlist) {
  setwd("~/nas/04.Results/drug/depmap/")
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  Cancername = gsub('TCGA-' , '', CancerType)
  
  # call input
  gc_cellline = readRDS(paste0("~/nas/00.data/filtered_TCGA/", num_CancerType, "/",Cancername,"_cellline_dual_all_log.rds"))
  gc_TCGA = readRDS(paste0("~/nas/04.Results/short_long/", CancerType,"_critical_features_short_long.rds"))
  
  # short_cluster = gc_TCGA %>% filter(cluster == "short") %>% dplyr::select(-vitalstatus, -duration,-status, -cluster)
  # long_cluster = gc_TCGA %>% filter(cluster == "long") %>% dplyr::select(-vitalstatus, -duration,-status, -cluster)
  wo_info = gc_TCGA %>% dplyr::select(-vitalstatus, -duration, -status)
  
  gc_cellline$cluster = "cell"
  gc_cellline_filt_df = gc_cellline[,colnames(wo_info)]
  
  gc_cellline_fixed_df = readRDS(paste0("~/nas/04.Results/drug/depmap/gdsc/",CancerType, "_DM_sl_cluster.rds"))
  info_predict = read.csv(paste0("~/nas/04.Results/cell/",Cancername,"_cellline_predict_for_bestmodel.csv"))
  ### Dimension reduction (PCA)
  # info_predict$adjust_predic
  if (all.equal(rownames(gc_cellline_fixed_df) , info_predict$cellline_id)) {
    gc_cellline_fixed_df$predict = info_predict$adjust_predic
  }
  
  set.seed(1)
  
  # combine_df = rbind(wo_info, gc_cellline_filt_df)
  
  pca_res_tcga <- prcomp(wo_info %>% dplyr::select(-cluster), retx=TRUE, center=TRUE, scale=TRUE)
  # pca_res_combine <- prcomp(combine_df %>% dplyr::select(-cluster), retx=TRUE, center=TRUE, scale=TRUE)
  
  cellline_pca_res = as.data.frame(predict(pca_res_tcga, gc_cellline_filt_df %>% select(-cluster)))
  # cellline_pca_res$cluster = gc_cellline_filt_df$cluster
  # cellline_pca_res$adjust_predict = gc_cellline_fixed_df$predict
  gc_cellline_fixed_df = gc_cellline_fixed_df %>% mutate(cluster_rename = case_when(cluster == "short" ~ "short_labelled",
                                                                                    cluster == "long" ~ "long_labelled",
                                                                                    .default = NA))
  cellline_pca_res$cluster_adjust = gc_cellline_fixed_df$cluster_rename
  cellline_pca_res$alpha = T
  cellline_pca_res$predict_state = gc_cellline_fixed_df$predict
  
  library(broom)
  variance_exp <- pca_res_tcga %>%  
    tidy("pcs") %>% 
    pull(percent)
  pca_TCGA = as.data.frame(pca_res_tcga$x)
  pca_TCGA$cluster_adjust = gc_TCGA$cluster
  pca_TCGA$alpha = F
  pca_TCGA$predict_state = "TCGA"
  
  merge_pca_res = rbind(pca_TCGA %>% 
                          select(PC1,PC2,cluster_adjust,alpha,predict_state) ,
                        cellline_pca_res %>% 
                          select(PC1,PC2,cluster_adjust,alpha,predict_state))

  # merge_pca_res = merge_pca_res %>% mutate(alpha = ifelse(grepl('TCGA', rownames(merge_pca_res)) , F, T),
  #                                          cluster_adjust = ifelse(grepl('ACH', rownames(.) ) , gc_cellline_fixed_df$cluster_rename, cluster),
  #                                          predict_state = ifelse(grepl('ACH', rownames(.) ) , gc_cellline_fixed_df$predict, "TCGA"))
  # 
  mean_of_long = pca_TCGA %>% filter(cluster_adjust  == "long") %>% select(PC1, PC2) %>%
    colMeans() %>% t() %>% as.data.frame()
  rownames(mean_of_long) = "TCGA-Long"
  
  mean_of_short = pca_TCGA %>% filter(cluster_adjust  == "short") %>% select(PC1, PC2) %>%
    colMeans() %>% t() %>% as.data.frame()
  rownames(mean_of_short) = "TCGA-Short"
  
  mean_of_long$alpha = T
  mean_of_long$cluster_adjust = "TCGA-Long"
  mean_of_long$predict_state = "TCGA_mean"
  
  mean_of_short$alpha = T
  mean_of_short$cluster_adjust = "TCGA-Short"
  mean_of_short$predict_state = "TCGA_mean"
  
  merge_pca_res = rbind(merge_pca_res, mean_of_long, mean_of_short)

  # manual_size =  ifelse(merge_pca_res$predict_state == "TCGA" , 2, 4)
  manual_size = ifelse(merge_pca_res$predict_state == "TCGA_mean" , 5, 3)
  merge_pca_res$text = ifelse(merge_pca_res$predict_state == "TCGA_mean", rownames(merge_pca_res),NA)
  
  # merge_pca_res$text = ifelse(merge_pca_res$predict_state == "TCGA_mean" | rownames(merge_pca_res) %in% test, rownames(merge_pca_res),NA)
  
  merge_pca_res$binary = ifelse(grepl("short",merge_pca_res$cluster_adjust) ,1 ,0)
  
  test_set = merge_pca_res %>% filter(!grepl( "TCGA" , predict_state))
  train_set = merge_pca_res %>% filter(predict_state == "TCGA")
  
  svmfit = svm(binary ~ PC1 + PC2, 
               data = train_set, 
               kernel = "linear",
               cost = 10,
               scale = FALSE,
               type = "C-classification")
  preds = predict(svmfit, test_set %>% select(PC1,PC2), decision.values = T)
  
  x1 = seq(min(train_set$PC1) - 1, max(train_set$PC1) + 1, by = 0.01)
  y1 = seq(min(train_set$PC2) - 1, max(train_set$PC2) + 1, by = 0.01)
  
  grid <- expand.grid(x1, y1)
  colnames(grid) = c('PC1','PC2')
  
  grid$value = predict(svmfit, grid)
  grid$z = as.vector(attributes(predict(svmfit, grid, decision.values = TRUE))$decision)

  grid$cluster_adjust = "line"
  grid$alpha = F
  grid$text = ""
  
  setwd("~/nas/04.Results/PCA_tsne/cellline_pca/")
  pca_plot = ggplot(merge_pca_res , aes(x = PC1, 
                                        y = PC2, 
                                        color = cluster_adjust,
                                        alpha = alpha,
                                        label = text))+
    geom_point(size = manual_size)+
    # geom_point(data = subset(merge_pca_res, predict_state %in% c("Alive" , "Dead")),
    #            col = "black", 
    #            stroke = 1,
    #            size = 3 ,
    #            shape = 21) +
    geom_point(data = subset(merge_pca_res, predict_state == "TCGA_mean"),
               col = "black", 
               stroke = 1,
               size = 5 ,
               shape = 21) +
    scale_alpha_discrete(range = c(0.3, 1)) +
    # geom_text_repel() +
    geom_text_repel(color = "black",
                    point.padding = 0.5,    
                    nudge_x = .3,
                    nudge_y = .35) +
    stat_contour(data = grid, aes(x = PC1, y = PC2, z = z), breaks = c(0)) +
    scale_shape_manual(values = c(19)) + 
    labs(x = paste0("PC1: ",round(variance_exp[1]*100), "%"),
         y = paste0("PC2: ",round(variance_exp[2]*100), "%")) +
    theme_classic()+  

  scale_color_manual(values=c("short_labelled" = "black",
                              "long_labelled" = "black",
                              "long" = "#4DAF4A",
                              "short" = "#E41A1C",
                              "TCGA-Long" = "#4dbbd5",
                              "TCGA-Short" = "#990066",
                              "line" = "black")) +
    
    theme(legend.position = "none",
          text = element_text(face = "bold" ,  family = "Helvetica-Bold", size = 10), 
          # axis.text.x = element_text(face = "bold", family = "Helvetica-Bold", size = 10),
          # axis.text.y = element_text(face = "bold", family = "Helvetica-Bold", size = 10),
          # axis.title.x = element_text(face = "bold", family = "Helvetica-Bold", size = 10),
          axis.title.y = element_text(face = "bold", family = "Helvetica-Bold", size = 10)) +
    annotate(geom = "text" ,
             x = min(merge_pca_res$PC1) + 0.5 , 
             y = max(merge_pca_res$PC2) - 0.5 , 
             label = "Long Group", 
             color = "#4DAF4A" ,
             size = 10 , 
             # face = "bold", 
             family = "Helvetica") +
    annotate(geom = "text" ,
             x = max(merge_pca_res$PC1) - 0.5 , 
             y = min(merge_pca_res$PC2) + 0.5 , 
             label = "Short Group", 
             color = "#E41A1C" ,
             size = 10 , 
             # face = "bold", 
             family = "Helvetica")

  ggsave(file=paste0(CancerType , "_PCA_plot_wo_predict_wo_label.svg"), plot=pca_plot, width=5, height=5)
  # write.csv(merge_pca_res,paste0(CancerType , "_PCA_.csv" ))
}
