library(survival) 
library(dplyr)
library(survminer)
library(pheatmap)
library(RColorBrewer)
library(readxl)
library(openxlsx)

filepath = "/home/seokwon/nas/"
ref_path = paste0(filepath, "99.reference/")

Cancerlist = dir(paste0(filepath, "/00.data/filtered_TCGA/"))

setwd("~/nas/04.Results/short_long/ttest")
# num_CancerType = "19.TCGA-LIHC"
# for all
for (num_CancerType in Cancerlist) {
  
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  
  # call input
  
  # duration_log_df = readRDS(paste0(main.path_tc, "/", CancerType,"_dual_add_duration_log.rds"))
  cancer_short_long_features = readRDS(paste0("~/nas/04.Results/short_long/",CancerType,"_critical_features_short_long.rds"))
  
  tmp_meta = cancer_short_long_features[,which(colnames(cancer_short_long_features) %in% c("vitalstatus" , "duration", "status", "cluster"))]
  rownames(tmp_meta) <- sub("TCGA",CancerType , rownames(tmp_meta))
  tmp_meta$pat = rownames(tmp_meta)
  rownames(tmp_meta) = NULL 
  tmp_meta = tmp_meta[, c("pat", setdiff(names(tmp_meta), "pat"))]
  tmp_meta$cancertype = CancerType
  tmp_meta = tmp_meta[, c("cancertype", setdiff(names(tmp_meta), "cancertype"))]
  
  if (CancerType == "TCGA-CESC") {
    total_meta = tmp_meta
  } else {
    total_meta = rbind(total_meta,tmp_meta)
    
  }
  
  cancer_short_long_features = cancer_short_long_features[,which(!colnames(cancer_short_long_features) %in% c("vitalstatus" , "duration", "status", "cluster"))]
  cancer_short_long_features_t = as.data.frame(t(cancer_short_long_features))
  
  colnames(cancer_short_long_features_t) <- sub("TCGA",CancerType , colnames(cancer_short_long_features_t))
  cancer_short_long_features_t$features = rownames(cancer_short_long_features_t)
  rownames(cancer_short_long_features_t) = NULL
  cancer_short_long_features_t = cancer_short_long_features_t[, c("features", setdiff(names(cancer_short_long_features_t), "features"))]
  
  if (CancerType == "TCGA-CESC") {
    total_pca_df = cancer_short_long_features_t
  } else {
    total_pca_df = merge(total_pca_df,cancer_short_long_features_t, by = "features", all = TRUE)
    
  }
  
  remove(cancer_short_long_features_t)
  
}

total_pca_df[is.na(total_pca_df)] = 0

rownames(total_pca_df) = total_pca_df$features
total_pca_df$features = NULL 
total_pca_df_filt = total_pca_df[,which(colSums(total_pca_df) != 0)]

# pca <- prcomp(total_pca_df_filt , scale. = TRUE)
pca <- prcomp(total_pca_df_filt, scale. = TRUE, retx = TRUE, rank. = 2)

pca_rank2 = as.data.frame(pca$rotation)
pca_rank2$pat = rownames(pca_rank2)
rownames(pca_rank2) = NULL 
pca_rank2$cancertype = sub("^(TCGA-.+?)-.*$", "\\1", pca_rank2$pat)
# pca_rank2$cancerpathway = NA
pca_rank2[is.na(pca_rank2)] = 0

pca_rank2_add_meta = merge(pca_rank2, total_meta, by = "pat", all = TRUE)
pca_rank2_add_meta$cancertype.y = NULL
pca_rank2_add_meta = pca_rank2_add_meta[which(!is.na(pca_rank2_add_meta$PC1)),]

num_pat_n = brewer.pal(length(unique(pca_rank2_add_meta$cancertype.x)), "Paired")
col_pat_n = setNames(num_pat_n, unique(pca_rank2_add_meta$cancertype.x))

ggplot(pca_rank2_add_meta, aes(PC1, PC2 ,shape = factor(cluster), color = factor(cancertype.x) ))+
  geom_point() +
  scale_color_manual(values = col_pat_n)
# stat_ellipse(aes(color = factor(cancertype.x)), geom = "path", linewidth = 1, alpha = 0.5)
