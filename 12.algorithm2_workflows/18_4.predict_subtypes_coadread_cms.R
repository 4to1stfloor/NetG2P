library(readxl)
library(devtools)
library(CMSclassifier)
library(org.Hs.eg.db)
library(data.table)
library(scater)
library(scran)
library(pheatmap)
library(dplyr)
library(stringr)
library(RColorBrewer)
setwd("~/nas/04.Results/subtypes/")

filepath = "/home/seokwon/nas/"
ref_path = paste0(filepath, "99.reference/")

cli = fread("~/nas/99.reference/all_clin_indexed.csv")

# cms_genes = read_xlsx("~/nas/00.data/filtered_TCGA/34.TCGA-COADREAD/CMS_genes_merge.xlsx")
# cms = readRDS("~/nas/00.data/filtered_TCGA/34.TCGA-COADREAD/TCGA-COADREAD_CMS.rds")
cms = readRDS('/mnt/gluster_server/data/raw/TCGA_data/00.data/34.TCGA-COADREAD/TCGA-COADREAD.cms.tbl.rds')

num_CancerType = "34.TCGA-COADREAD"
main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))

# call input
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

cancer_bf = read.csv(paste0(filepath,"04.Results/bestfeatures/",CancerType, "_best_features.csv"))

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

cms_filt = cms[grep("*01A-*" ,cms$barcode ),]
cms_filt$patient = substr(cms_filt$barcode, 1,12)
dup_cms = cms_filt[which(duplicated(cms_filt$patient)),]

for (dup_pat in dup_cms$patient) {
  tmp_dup = cms_filt[which(cms_filt$patient == dup_pat),]
  if (length(unique(tmp_dup$CMStype) == 1)) {
    tmp_dup[-nrow(tmp_dup),]
    cms_filt = cms_filt[which(!cms_filt$barcode == tmp_dup$barcode[-nrow(tmp_dup)]),]
  }
}

rownames(cms_filt) = cms_filt$patient

# mapping cms type
# dual
dual_bf_filt_inter = dual_bf_filt[intersect( rownames(cms_filt) , rownames(dual_bf_filt)),]
cms_dual_filt = cms_filt[intersect( rownames(cms_filt) , rownames(dual_bf_filt)),]

if (all.equal(rownames(dual_bf_filt_inter) , rownames(cms_dual_filt)) ) {
  dual_bf_filt_inter$CMS = cms_dual_filt$CMStype
}

# exp
exp_bf_filt_inter = exp_bf_filt[intersect( rownames(cms_filt) , rownames(exp_bf_filt)),]
cms_exp_filt = cms_filt[intersect( rownames(cms_filt) , rownames(exp_bf_filt)),]

if (all.equal(rownames(exp_bf_filt_inter) , rownames(cms_exp_filt)) ) {
  exp_bf_filt_inter$CMS = cms_exp_filt$CMStype
}

# mut
mut_bf_filt_inter = mut_bf_filt[intersect( rownames(cms_filt) , rownames(mut_bf_filt)),]
cms_mut_filt = cms_filt[intersect( rownames(cms_filt) , rownames(mut_bf_filt)),]

if (all.equal(rownames(mut_bf_filt_inter) , rownames(cms_mut_filt)) ) {
  mut_bf_filt_inter$CMS = cms_mut_filt$CMStype
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

dual_bf_filt_inter_sl = dual_bf_filt_inter_sl[which(!is.na(dual_bf_filt_inter_sl$CMS)),]
exp_bf_filt_inter_sl = exp_bf_filt_inter_sl[which(!is.na(exp_bf_filt_inter_sl$CMS)),]
mut_bf_filt_inter_sl = mut_bf_filt_inter_sl[which(!is.na(mut_bf_filt_inter_sl$CMS)),]

dual_bf_filt_inter_sl = dual_bf_filt_inter_sl[order(dual_bf_filt_inter_sl$CMS),]
exp_bf_filt_inter_sl = exp_bf_filt_inter_sl[order(exp_bf_filt_inter_sl$CMS),]
mut_bf_filt_inter_sl = mut_bf_filt_inter_sl[order(mut_bf_filt_inter_sl$CMS),]

# pic
# dual
annotation_df <- data.frame(CMS = dual_bf_filt_inter_sl$CMS,
                            cluster = dual_bf_filt_inter_sl$cluster)
rownames(annotation_df) <- rownames(dual_bf_filt_inter_sl)

# Create a named color vector for the unique values of vital_status
num_CMS = brewer.pal(length(unique(dual_bf_filt_inter_sl$CMS)), "Spectral")
col_CMS = setNames(num_CMS, unique(dual_bf_filt_inter_sl$CMS))

cluster_colors <- c("short" = "yellow", "long" = "black")
subtypes_colors = list(CMS = col_CMS , cluster = cluster_colors)

png(filename = paste0(CancerType,"_dual_features_w_short_long_clusterT_CMS.png"),
    width = 25, height = 25,  units = "cm" ,pointsize = 12,
    bg = "white", res = 1200, family = "")

tmp = pheatmap::pheatmap(as.matrix(t(dual_bf_filt_inter_sl[,which(!colnames(dual_bf_filt_inter_sl) %in% c("CMS","cluster"))])),
                         annotation_col = annotation_df,
                         annotation_colors = subtypes_colors,
                         cluster_cols = T)
print(tmp)
dev.off()

png(filename = paste0(CancerType,"_dual_features_w_short_long_clusterF_CMS.png"),
    width = 35, height = 35,  units = "cm" ,pointsize = 12,
    bg = "white", res = 1200, family = "")

tmp2 = pheatmap::pheatmap(as.matrix(t(dual_bf_filt_inter_sl[,which(!colnames(dual_bf_filt_inter_sl) %in% c("CMS","cluster"))])),
                          annotation_col = annotation_df,
                          annotation_colors = subtypes_colors,
                          cluster_cols = F)
print(tmp2)

dev.off()

png(filename = paste0(CancerType,"_dual_features_w_short_long_complex_CMS.png"),
    width = 35, height = 35,  units = "cm" ,pointsize = 12,
    bg = "white", res = 1200, family = "")
tmp3 = ComplexHeatmap::pheatmap(as.matrix(t(dual_bf_filt_inter_sl[,which(!colnames(dual_bf_filt_inter_sl) %in% c("CMS","cluster"))])),
                                column_split = annotation_df$CMS,
                                annotation_col = annotation_df,
                                annotation_colors = subtypes_colors,
                                cluster_cols = T)

print(tmp3)
dev.off()

# exp
if (sum(colSums(exp_bf_filt_inter_sl[,which(!colnames(exp_bf_filt_inter_sl) %in% c("CMS","cluster"))]) == 0) == 
    ncol(exp_bf_filt_inter_sl[,which(!colnames(exp_bf_filt_inter_sl) %in%  c("CMS","cluster"))])) {
  print("There are not a value in dataframe.")
  } else {
    annotation_df <- data.frame(CMS = exp_bf_filt_inter_sl$CMS,
                                cluster = exp_bf_filt_inter_sl$cluster)
    rownames(annotation_df) <- rownames(exp_bf_filt_inter_sl)
    
    # Create a named color vector for the unique values of vital_status
    num_CMS = brewer.pal(length(unique(exp_bf_filt_inter_sl$CMS)), "Spectral")
    col_CMS = setNames(num_CMS, unique(exp_bf_filt_inter_sl$CMS))
    
    cluster_colors <- c("short" = "yellow", "long" = "black")
    subtypes_colors = list(CMS = col_CMS , cluster = cluster_colors)
    
    png(filename = paste0(CancerType,"_exp_features_w_short_long_clusterT_CMS.png"),
        width = 25, height = 25,  units = "cm" ,pointsize = 12,
        bg = "white", res = 1200, family = "")
    
    tmp = pheatmap::pheatmap(as.matrix(t(exp_bf_filt_inter_sl[,which(!colnames(exp_bf_filt_inter_sl) %in% c("CMS","cluster"))])),
                             annotation_col = annotation_df,
                             annotation_colors = subtypes_colors,
                             cluster_cols = T)
    print(tmp)
    dev.off()
    
    png(filename = paste0(CancerType,"_exp_features_w_short_long_clusterF_CMS.png"),
        width = 35, height = 35,  units = "cm" ,pointsize = 12,
        bg = "white", res = 1200, family = "")
    
    tmp2 = pheatmap::pheatmap(as.matrix(t(exp_bf_filt_inter_sl[,which(!colnames(exp_bf_filt_inter_sl) %in% c("CMS","cluster"))])),
                              annotation_col = annotation_df,
                              annotation_colors = subtypes_colors,
                              cluster_cols = F)
    print(tmp2)
    
    dev.off()
    
    png(filename = paste0(CancerType,"_exp_features_w_short_long_complex_CMS.png"),
        width = 35, height = 35,  units = "cm" ,pointsize = 12,
        bg = "white", res = 1200, family = "")
    tmp3 = ComplexHeatmap::pheatmap(as.matrix(t(exp_bf_filt_inter_sl[,which(!colnames(exp_bf_filt_inter_sl) %in% c("CMS","cluster"))])),
                                    column_split = annotation_df$CMS,
                                    annotation_col = annotation_df,
                                    annotation_colors = subtypes_colors,
                                    cluster_cols = T)
    
    print(tmp3)
    dev.off()
}
# mut


if (sum(colSums(mut_bf_filt_inter_sl[,which(!colnames(mut_bf_filt_inter_sl) %in% c("CMS","cluster"))]) == 0) == 
    ncol(mut_bf_filt_inter_sl[,which(!colnames(mut_bf_filt_inter_sl) %in%  c("CMS","cluster"))])) {
  print("There are not a value in dataframe.")
  } else {
    
    annotation_df <- data.frame(CMS = mut_bf_filt_inter_sl$CMS,
                                cluster = mut_bf_filt_inter_sl$cluster)
    rownames(annotation_df) <- rownames(mut_bf_filt_inter_sl)
    
    # Create a named color vector for the unique values of vital_status
    num_CMS = brewer.pal(length(unique(mut_bf_filt_inter_sl$CMS)), "Spectral")
    col_CMS = setNames(num_CMS, unique(mut_bf_filt_inter_sl$CMS))
    
    cluster_colors <- c("short" = "yellow", "long" = "black")
    subtypes_colors = list(CMS = col_CMS , cluster = cluster_colors)
    
    png(filename = paste0(CancerType,"_mut_features_w_short_long_clusterT_CMS.png"),
        width = 25, height = 25,  units = "cm" ,pointsize = 12,
        bg = "white", res = 1200, family = "")
    
    tmp = pheatmap::pheatmap(as.matrix(t(mut_bf_filt_inter_sl[,which(!colnames(mut_bf_filt_inter_sl) %in% c("CMS","cluster"))])),
                             annotation_col = annotation_df,
                             annotation_colors = subtypes_colors,
                             cluster_cols = T)
    print(tmp)
    dev.off()
    
    png(filename = paste0(CancerType,"_mut_features_w_short_long_clusterF_CMS.png"),
        width = 35, height = 35,  units = "cm" ,pointsize = 12,
        bg = "white", res = 1200, family = "")
    
    tmp2 = pheatmap::pheatmap(as.matrix(t(mut_bf_filt_inter_sl[,which(!colnames(mut_bf_filt_inter_sl) %in% c("CMS","cluster"))])),
                              annotation_col = annotation_df,
                              annotation_colors = subtypes_colors,
                              cluster_cols = F)
    print(tmp2)
    
    dev.off()
    
    png(filename = paste0(CancerType,"_mut_features_w_short_long_complex_CMS.png"),
        width = 35, height = 35,  units = "cm" ,pointsize = 12,
        bg = "white", res = 1200, family = "")
    tmp3 = ComplexHeatmap::pheatmap(as.matrix(t(mut_bf_filt_inter_sl[,which(!colnames(mut_bf_filt_inter_sl) %in% c("CMS","cluster"))])),
                                    column_split = annotation_df$CMS,
                                    annotation_col = annotation_df,
                                    annotation_colors = subtypes_colors,
                                    cluster_cols = T)
    
    print(tmp3)
    dev.off()
}
# mut has data -> 0