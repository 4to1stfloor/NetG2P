library(readxl)
library(devtools)
library(CMSclassifier)
library(org.Hs.eg.db)
library(data.table)
library(scater)
library(scran)
library(pheatmap)
require(TCGAbiolinks)
require(SummarizedExperiment)
require(SingleCellExperiment)
library(dplyr)
library(RColorBrewer)
setwd("~/nas/04.Results/subtypes/")

filepath = "/home/seokwon/nas/"
ref_path = paste0(filepath, "99.reference/")
sce_path = "/mnt/gluster_server/data/raw/TCGA_data/00.data/"

num_CancerType = "30.TCGA-BRCA"

# for all

main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))

# call input

sce = readRDS(paste0(sce_path,num_CancerType,"/", CancerType,".sce_raw_snv_cnv.rds"))
sce_df = as.data.frame(colData(sce))

cancer_bf = read.csv(paste0(filepath,"04.Results/bestfeatures/",CancerType, "_best_features.csv"))
duration_log_df = readRDS(paste0(main.path_tc, "/", CancerType,"_dual_add_duration_log.rds"))

exp_each_log_df = readRDS(paste0(main.path_tc, "/", CancerType,"_exp_pathwayeach_all_log.rds"))
exp_link_log_df = readRDS(paste0(main.path_tc, "/", CancerType,"_exp_pathwaylink_all_log.rds"))

mut_each_log_df = readRDS(paste0(main.path_tc, "/", CancerType,"_mut_pathwayeach_all_log.rds"))
mut_link_log_df = readRDS(paste0(main.path_tc, "/", CancerType,"_mut_pathwaylink_all_log.rds"))

# exp
exp_link_log_df_filt_wo = exp_link_log_df[intersect(rownames(exp_link_log_df), rownames(exp_each_log_df)),]
exp_each_log_df_filt_wo = exp_each_log_df[intersect(rownames(exp_link_log_df), rownames(exp_each_log_df)),]

if (all.equal(rownames(exp_link_log_df_filt_wo), rownames(exp_each_log_df_filt_wo))) {
  exp_link_log_df_filt_wo$vitalstatus = NULL
  exp_dual_log_df_filt_wo = cbind(exp_link_log_df_filt_wo,exp_each_log_df_filt_wo)
}

# mut

mut_link_log_df_filt_wo = mut_link_log_df[intersect(rownames(mut_link_log_df), rownames(mut_each_log_df)),]
mut_each_log_df_filt_wo = mut_each_log_df[intersect(rownames(mut_link_log_df), rownames(mut_each_log_df)),]

if (all.equal(rownames(mut_link_log_df_filt_wo), rownames(mut_each_log_df_filt_wo))) {
  mut_link_log_df_filt_wo$vitalstatus = NULL
  mut_dual_log_df_filt_wo = cbind(mut_link_log_df_filt_wo,mut_each_log_df_filt_wo)
}

dual_bf_filt = duration_log_df[,cancer_bf$variable]
exp_bf_filt = exp_dual_log_df_filt_wo[,cancer_bf$variable]
mut_bf_filt = mut_dual_log_df_filt_wo[,cancer_bf$variable]

# dual
tmp_numeric <- matrix(as.numeric(unlist(dual_bf_filt)), nrow = nrow(dual_bf_filt))
vec <- as.vector(tmp_numeric)
dual_bf_filt[dual_bf_filt > quantile(vec, probs = 0.75) + 1.5*IQR(vec)] <- round(quantile(vec, probs = 0.75) + 1.5*IQR(vec),1)

# exp
tmp_numeric <- matrix(as.numeric(unlist(exp_bf_filt)), nrow = nrow(exp_bf_filt))
vec <- as.vector(tmp_numeric)
exp_bf_filt[exp_bf_filt > quantile(vec, probs = 0.75) + 1.5*IQR(vec)] <- round(quantile(vec, probs = 0.75) + 1.5*IQR(vec),1)

# mut
tmp_numeric <- matrix(as.numeric(unlist(mut_bf_filt)), nrow = nrow(mut_bf_filt))
vec <- as.vector(tmp_numeric)
mut_bf_filt[mut_bf_filt > quantile(vec, probs = 0.75) + 1.5*IQR(vec)] <- round(quantile(vec, probs = 0.75) + 1.5*IQR(vec),1)

# filter patient
pam50_filt = sce_df[grep("*01A-*" ,sce_df$barcode ),]
pam50_filt$patient = substr(pam50_filt$barcode, 1,12)
dup_pam50 = pam50_filt[which(duplicated(pam50_filt$patient)),]

for (dup_pat in dup_pam50$patient) {
  tmp_dup = pam50_filt[which(pam50_filt$patient == dup_pat),]
  if (length(unique(tmp_dup$paper_BRCA_Subtype_PAM50) == 1)) {
    tmp_dup[-nrow(tmp_dup),]
    pam50_filt = pam50_filt[which(!pam50_filt$barcode == tmp_dup$barcode[-nrow(tmp_dup)]),]
  }
}

rownames(pam50_filt) = pam50_filt$patient

# mapping pam50 type
# dual
dual_bf_filt_inter = dual_bf_filt[intersect( rownames(pam50_filt) , rownames(dual_bf_filt)),]
pam50_dual_filt = pam50_filt[intersect( rownames(pam50_filt) , rownames(dual_bf_filt)),]

if (all.equal(rownames(dual_bf_filt_inter) , rownames(pam50_dual_filt)) ) {
  dual_bf_filt_inter$pam50 = pam50_dual_filt$paper_BRCA_Subtype_PAM50
}

# exp
exp_bf_filt_inter = exp_bf_filt[intersect( rownames(pam50_filt) , rownames(exp_bf_filt)),]
pam50_exp_filt = pam50_filt[intersect( rownames(pam50_filt) , rownames(exp_bf_filt)),]

if (all.equal(rownames(exp_bf_filt_inter) , rownames(pam50_exp_filt)) ) {
  exp_bf_filt_inter$pam50 = pam50_exp_filt$paper_BRCA_Subtype_PAM50
}

# mut
mut_bf_filt_inter = mut_bf_filt[intersect( rownames(pam50_filt) , rownames(mut_bf_filt)),]
pam50_mut_filt = pam50_filt[intersect( rownames(pam50_filt) , rownames(mut_bf_filt)),]

if (all.equal(rownames(mut_bf_filt_inter) , rownames(pam50_mut_filt)) ) {
  mut_bf_filt_inter$pam50 = pam50_mut_filt$paper_BRCA_Subtype_PAM50
}

# mapping cluster 
short_long = readRDS(paste0(filepath,"04.Results/short_long/",CancerType,"_best_features_short_long.rds"))

# dual 
dual_bf_filt_inter_sl = dual_bf_filt_inter[intersect( rownames(short_long) , rownames(dual_bf_filt_inter)),]
short_long_dual_filt = short_long[intersect( rownames(short_long) , rownames(dual_bf_filt_inter)),]

if (all.equal(rownames(short_long_dual_filt), rownames(dual_bf_filt_inter_sl))) {
  dual_bf_filt_inter_sl$cluster = short_long_dual_filt$cluster
}

# exp
exp_bf_filt_inter_sl = exp_bf_filt_inter[intersect( rownames(short_long) , rownames(exp_bf_filt_inter)),]
short_long_exp_filt = short_long[intersect( rownames(short_long) , rownames(exp_bf_filt_inter)),]

if (all.equal(rownames(short_long_exp_filt), rownames(exp_bf_filt_inter_sl))) {
  exp_bf_filt_inter_sl$cluster = short_long_exp_filt$cluster
}

# mut

mut_bf_filt_inter_sl = mut_bf_filt_inter[intersect( rownames(short_long) , rownames(mut_bf_filt_inter)),]
short_long_mut_filt = short_long[intersect( rownames(short_long) , rownames(mut_bf_filt_inter)),]

if (all.equal(rownames(short_long_mut_filt), rownames(mut_bf_filt_inter_sl))) {
  mut_bf_filt_inter_sl$cluster = short_long_mut_filt$cluster
}

dual_bf_filt_inter_sl = dual_bf_filt_inter_sl[which(!is.na(dual_bf_filt_inter_sl$pam50)),]
exp_bf_filt_inter_sl = exp_bf_filt_inter_sl[which(!is.na(exp_bf_filt_inter_sl$pam50)),]
mut_bf_filt_inter_sl = mut_bf_filt_inter_sl[which(!is.na(mut_bf_filt_inter_sl$pam50)),]

dual_bf_filt_inter_sl = dual_bf_filt_inter_sl[order(dual_bf_filt_inter_sl$pam50),]
exp_bf_filt_inter_sl = exp_bf_filt_inter_sl[order(exp_bf_filt_inter_sl$pam50),]
mut_bf_filt_inter_sl = mut_bf_filt_inter_sl[order(mut_bf_filt_inter_sl$pam50),]

# pic
# dual
annotation_df <- data.frame(pam50 = dual_bf_filt_inter_sl$pam50,
                            cluster = dual_bf_filt_inter_sl$cluster)
rownames(annotation_df) <- rownames(dual_bf_filt_inter_sl)

# Create a named color vector for the unique values of vital_status
num_pam50 = brewer.pal(length(unique(dual_bf_filt_inter_sl$pam50)), "Spectral")
col_pam50 = setNames(num_pam50, unique(dual_bf_filt_inter_sl$pam50))

cluster_colors <- c("short" = "yellow", "long" = "black")
subtypes_colors = list(pam50 = col_pam50 , cluster = cluster_colors)

png(filename = paste0(CancerType,"_dual_features_w_short_long_clusterT_pam50.png"),
    width = 35, height = 35,  units = "cm" ,pointsize = 12,
    bg = "white", res = 1200, family = "")

tmp = pheatmap::pheatmap(as.matrix(t(dual_bf_filt_inter_sl[,which(!colnames(dual_bf_filt_inter_sl) %in% c("pam50","cluster"))])),
                         annotation_col = annotation_df,
                         annotation_colors = subtypes_colors,
                         cluster_cols = T)
print(tmp)
dev.off()

png(filename = paste0(CancerType,"_dual_features_w_short_long_clusterF_pam50.png"),
    width = 35, height = 35,  units = "cm" ,pointsize = 12,
    bg = "white", res = 1200, family = "")

tmp2 = pheatmap::pheatmap(as.matrix(t(dual_bf_filt_inter_sl[,which(!colnames(dual_bf_filt_inter_sl) %in% c("pam50","cluster"))])),
                          annotation_col = annotation_df,
                          annotation_colors = subtypes_colors,
                          cluster_cols = F)
print(tmp2)

dev.off()

png(filename = paste0(CancerType,"_dual_features_w_short_long_complex_pam50.png"),
    width = 35, height = 35,  units = "cm" ,pointsize = 12,
    bg = "white", res = 1200, family = "")
tmp3 = ComplexHeatmap::pheatmap(as.matrix(t(dual_bf_filt_inter_sl[,which(!colnames(dual_bf_filt_inter_sl) %in% c("pam50","cluster"))])),
                                column_split = annotation_df$pam50,
                                annotation_col = annotation_df,
                                annotation_colors = subtypes_colors,
                                cluster_cols = T)

print(tmp3)
dev.off()

# exp

if (sum(colSums(exp_bf_filt_inter_sl[,which(!colnames(exp_bf_filt_inter_sl) %in% c("pam50","cluster"))]) == 0) == 
    ncol(exp_bf_filt_inter_sl[,which(!colnames(exp_bf_filt_inter_sl) %in% c("pam50","cluster"))])) {
  
  print("There are not a value in dataframe.")
  } else {
    
    annotation_df <- data.frame(pam50 = exp_bf_filt_inter_sl$pam50,
                                cluster = exp_bf_filt_inter_sl$cluster)
    rownames(annotation_df) <- rownames(exp_bf_filt_inter_sl)
    
    # Create a named color vector for the unique values of vital_status
    num_pam50 = brewer.pal(length(unique(exp_bf_filt_inter_sl$pam50)), "Spectral")
    col_pam50 = setNames(num_pam50, unique(exp_bf_filt_inter_sl$pam50))
    
    cluster_colors <- c("short" = "yellow", "long" = "black")
    subtypes_colors = list(pam50 = col_pam50 , cluster = cluster_colors)
    
    png(filename = paste0(CancerType,"_exp_features_w_short_long_clusterT_pam50.png"),
        width = 25, height = 25,  units = "cm" ,pointsize = 12,
        bg = "white", res = 1200, family = "")
    
    tmp = pheatmap::pheatmap(as.matrix(t(exp_bf_filt_inter_sl[,which(!colnames(exp_bf_filt_inter_sl) %in% c("pam50","cluster"))])),
                             annotation_col = annotation_df,
                             annotation_colors = subtypes_colors,
                             cluster_cols = T)
    print(tmp)
    dev.off()
    
    png(filename = paste0(CancerType,"_exp_features_w_short_long_clusterF_pam50.png"),
        width = 35, height = 35,  units = "cm" ,pointsize = 12,
        bg = "white", res = 1200, family = "")
    
    tmp2 = pheatmap::pheatmap(as.matrix(t(exp_bf_filt_inter_sl[,which(!colnames(exp_bf_filt_inter_sl) %in% c("pam50","cluster"))])),
                              annotation_col = annotation_df,
                              annotation_colors = subtypes_colors,
                              cluster_cols = F)
    print(tmp2)
    
    dev.off()
    
    png(filename = paste0(CancerType,"_exp_features_w_short_long_complex_pam50.png"),
        width = 35, height = 35,  units = "cm" ,pointsize = 12,
        bg = "white", res = 1200, family = "")
    tmp3 = ComplexHeatmap::pheatmap(as.matrix(t(exp_bf_filt_inter_sl[,which(!colnames(exp_bf_filt_inter_sl) %in% c("pam50","cluster"))])),
                                    column_split = annotation_df$pam50,
                                    annotation_col = annotation_df,
                                    annotation_colors = subtypes_colors,
                                    cluster_cols = T)
    
    print(tmp3)
    dev.off()
  }

# mut

if (sum(colSums(mut_bf_filt_inter_sl[,which(!colnames(mut_bf_filt_inter_sl) %in% c("pam50","cluster"))]) == 0) == 
    ncol(mut_bf_filt_inter_sl[,which(!colnames(mut_bf_filt_inter_sl) %in% c("pam50","cluster"))])) {
  
  print("There are not a value in dataframe.")
  
  } else {
    
    annotation_df <- data.frame(pam50 = mut_bf_filt_inter_sl$pam50,
                                cluster = mut_bf_filt_inter_sl$cluster)
    rownames(annotation_df) <- rownames(mut_bf_filt_inter_sl)
    
    # Create a named color vector for the unique values of vital_status
    num_pam50 = brewer.pal(length(unique(mut_bf_filt_inter_sl$pam50)), "Spectral")
    col_pam50 = setNames(num_pam50, unique(mut_bf_filt_inter_sl$pam50))
    
    cluster_colors <- c("short" = "yellow", "long" = "black")
    subtypes_colors = list(pam50 = col_pam50 , cluster = cluster_colors)
    
    png(filename = paste0(CancerType,"_mut_features_w_short_long_clusterT_pam50.png"),
        width = 25, height = 25,  units = "cm" ,pointsize = 12,
        bg = "white", res = 1200, family = "")
    
    tmp = pheatmap::pheatmap(as.matrix(t(mut_bf_filt_inter_sl[,which(!colnames(mut_bf_filt_inter_sl) %in% c("pam50","cluster"))])),
                             annotation_col = annotation_df,
                             annotation_colors = subtypes_colors,
                             cluster_cols = T)
    print(tmp)
    dev.off()
    
    png(filename = paste0(CancerType,"_mut_features_w_short_long_clusterF_pam50.png"),
        width = 35, height = 35,  units = "cm" ,pointsize = 12,
        bg = "white", res = 1200, family = "")
    
    tmp2 = pheatmap::pheatmap(as.matrix(t(mut_bf_filt_inter_sl[,which(!colnames(mut_bf_filt_inter_sl) %in% c("pam50","cluster"))])),
                              annotation_col = annotation_df,
                              annotation_colors = subtypes_colors,
                              cluster_cols = F)
    print(tmp2)
    
    dev.off()
    
    png(filename = paste0(CancerType,"_mut_features_w_short_long_complex_pam50.png"),
        width = 35, height = 35,  units = "cm" ,pointsize = 12,
        bg = "white", res = 1200, family = "")
    tmp3 = ComplexHeatmap::pheatmap(as.matrix(t(mut_bf_filt_inter_sl[,which(!colnames(mut_bf_filt_inter_sl) %in% c("pam50","cluster"))])),
                                    column_split = annotation_df$pam50,
                                    annotation_col = annotation_df,
                                    annotation_colors = subtypes_colors,
                                    cluster_cols = T)
    
    print(tmp3)
    dev.off()
  }
