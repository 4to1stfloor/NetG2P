library(org.Hs.eg.db)
library(data.table)
library(scater)
library(scran)
library(pheatmap)
library(devtools)
library(data.table)
library(scuttle)
library(scran)
library(RColorBrewer)
library(readxl)
filepath = "/home/seokwon/nas/"
ref_path = paste0(filepath, "99.reference/")

Cancerlist = dir(paste0(filepath, "/00.data/filtered_TCGA/"))
sce_path = "/mnt/gluster_server/data/raw/TCGA_data/00.data/"
cli = fread("~/nas/99.reference/all_clin_indexed.csv")
surv_total_results = read_xlsx("~/nas/04.Results/Total_results_survpval2.xlsx")

setwd("~/nas/04.Results/compare_surv/survival_analysis_features")

# for all
for (num_CancerType in Cancerlist) {
  
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  
  # call input
  # 
  # sce = readRDS(paste0(sce_path,num_CancerType,"/", CancerType,".sce_raw_snv_cnv.rds"))
  # sce_df = as.data.frame(colData(sce))
  
  cancer_bf = read.csv(paste0(filepath,"04.Results/bestfeatures/",CancerType, "_best_features.csv"))
  cancer_bf_cut = cancer_bf[1:surv_total_results[which(surv_total_results$CancerType == CancerType),]$num_of_features,]
  
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
  
  dual_bf_filt = duration_log_df[,cancer_bf_cut$variable]
  exp_bf_filt = exp_dual_log_df_filt_wo[,cancer_bf_cut$variable]
  mut_bf_filt = mut_dual_log_df_filt_wo[,cancer_bf_cut$variable]

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
  
  if (CancerType == "TCGA-COADREAD") {
    cli_surv = cli[cli$project %in% c("TCGA-COAD","TCGA-READ"),
                   c("submitter_id",
                     "vital_status",
                     "days_to_death",
                     "days_to_last_follow_up")]
  } else if (CancerType == "TCGA-KIDNEY") {
    cli_surv = cli[cli$project %in% c("TCGA-KIRP","TCGA-KIRC","TCGA-KICH"),
                   c("submitter_id",
                     "vital_status",
                     "days_to_death",
                     "days_to_last_follow_up")]
  } else {
    cli_surv = cli[cli$project == CancerType,
                   c("submitter_id",
                     "vital_status",
                     "days_to_death",
                     "days_to_last_follow_up")]
  }
  
  cli_surv$deceased = cli_surv$vital_status == "Dead"
  cli_surv$overall_survival = ifelse(cli_surv$deceased,
                                     cli_surv$days_to_death,
                                     cli_surv$days_to_last_follow_up)
  
  cli_surv = cli_surv[!is.na(cli_surv$vital_status),]
  cli_surv = cli_surv[!(cli_surv$vital_status == "Not Reported"),]
  
  cli_surv$status = NA
  cli_surv$status[which(cli_surv$vital_status == "Dead")] = 1
  cli_surv$status[which(cli_surv$vital_status == "Alive")] = 0
  
  cli_surv$status = as.numeric(cli_surv$status)
  
  # filter patient

  dup_cli = cli_surv[which(duplicated(cli_surv$submitter_id)),]
  
  for (dup_pat in dup_cli$submitter_id) {
    tmp_dup = cli_surv[which(cli_surv$submitter_id == dup_pat),]
    if (length(unique(tmp_dup$status) == 1)) {
      tmp_dup[-nrow(tmp_dup),]
      cli_surv = cli_surv[which(!cli_surv$submitter_id == tmp_dup$submitter_id[-nrow(tmp_dup)]),]
    }
  }
  cli_surv = as.data.frame(cli_surv)
  rownames(cli_surv) = cli_surv$submitter_id
  
  # mapping pam50 type
  # dual
  dual_bf_filt_inter = dual_bf_filt[intersect( rownames(cli_surv) , rownames(dual_bf_filt)),]
  cli_dual_filt = cli_surv[intersect( rownames(cli_surv) , rownames(dual_bf_filt)),]
  
  if (all.equal(rownames(dual_bf_filt_inter) , rownames(cli_dual_filt)) ) {
    dual_bf_filt_inter$vital_status = cli_dual_filt$vital_status
  }
  
  # exp
  exp_bf_filt_inter = exp_bf_filt[intersect( rownames(cli_surv) , rownames(exp_bf_filt)),]
  cli_exp_filt = cli_surv[intersect( rownames(cli_surv) , rownames(exp_bf_filt)),]
  
  if (all.equal(rownames(exp_bf_filt_inter) , rownames(cli_exp_filt)) ) {
    exp_bf_filt_inter$vital_status = cli_exp_filt$vital_status
  }
  
  # mut
  mut_bf_filt_inter = mut_bf_filt[intersect( rownames(cli_surv) , rownames(mut_bf_filt)),]
  cli_mut_filt = cli_surv[intersect( rownames(cli_surv) , rownames(mut_bf_filt)),]
  
  if (all.equal(rownames(mut_bf_filt_inter) , rownames(cli_mut_filt)) ) {
    mut_bf_filt_inter$vital_status = cli_mut_filt$vital_status
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
  
  dual_bf_filt_inter_sl = dual_bf_filt_inter_sl[which(!is.na(dual_bf_filt_inter_sl$vital_status)),]
  exp_bf_filt_inter_sl = exp_bf_filt_inter_sl[which(!is.na(exp_bf_filt_inter_sl$vital_status)),]
  mut_bf_filt_inter_sl = mut_bf_filt_inter_sl[which(!is.na(mut_bf_filt_inter_sl$vital_status)),]
  
  dual_bf_filt_inter_sl = dual_bf_filt_inter_sl[order(dual_bf_filt_inter_sl$vital_status),]
  exp_bf_filt_inter_sl = exp_bf_filt_inter_sl[order(exp_bf_filt_inter_sl$vital_status),]
  mut_bf_filt_inter_sl = mut_bf_filt_inter_sl[order(mut_bf_filt_inter_sl$vital_status),]
  
  # pic
  # dual
  annotation_df <- data.frame(vital_status = dual_bf_filt_inter_sl$vital_status,
                              cluster = dual_bf_filt_inter_sl$cluster)
  rownames(annotation_df) <- rownames(dual_bf_filt_inter_sl)
  
  # Create a named color vector for the unique values of vital_status
  num_vital_status = brewer.pal(3, "Spectral")[c(1,3)]
  col_vital_status = setNames(num_vital_status, unique(dual_bf_filt_inter_sl$vital_status))
  
  cluster_colors <- c("short" = "yellow", "long" = "black")
  subtypes_colors = list(vital_status = col_vital_status , cluster = cluster_colors)
  
  png(filename = paste0(CancerType,"_dual_sa_features_w_short_long_clusterT_surv.png"),
      width = 35, height = 35,  units = "cm" ,pointsize = 12,
      bg = "white", res = 1200, family = "")
  
  tmp = pheatmap::pheatmap(as.matrix(t(dual_bf_filt_inter_sl[,which(!colnames(dual_bf_filt_inter_sl) %in% c("vital_status","cluster"))])),
                           annotation_col = annotation_df,
                           annotation_colors = subtypes_colors,
                           cluster_cols = T)
  print(tmp)
  dev.off()
  
  png(filename = paste0(CancerType,"_dual_sa_features_w_short_long_clusterF_surv.png"),
      width = 35, height = 35,  units = "cm" ,pointsize = 12,
      bg = "white", res = 1200, family = "")
  
  tmp2 = pheatmap::pheatmap(as.matrix(t(dual_bf_filt_inter_sl[,which(!colnames(dual_bf_filt_inter_sl) %in% c("vital_status","cluster"))])),
                            annotation_col = annotation_df,
                            annotation_colors = subtypes_colors,
                            cluster_cols = F)
  print(tmp2)
  
  dev.off()
  
  png(filename = paste0(CancerType,"_dual_sa_features_w_short_long_complex_surv.png"),
      width = 35, height = 35,  units = "cm" ,pointsize = 12,
      bg = "white", res = 1200, family = "")
  tmp3 = ComplexHeatmap::pheatmap(as.matrix(t(dual_bf_filt_inter_sl[,which(!colnames(dual_bf_filt_inter_sl) %in% c("vital_status","cluster"))])),
                                  column_split = annotation_df$vital_status,
                                  annotation_col = annotation_df,
                                  annotation_colors = subtypes_colors,
                                  cluster_cols = T)
  
  print(tmp3)
  dev.off()
  
  # exp
  if (sum(colSums(exp_bf_filt_inter_sl[,which(!colnames(exp_bf_filt_inter_sl) %in% c("vital_status","cluster"))]) == 0) == 
      ncol(exp_bf_filt_inter_sl[,which(!colnames(exp_bf_filt_inter_sl) %in% c("vital_status","cluster"))])) {
    next
  } else {
    
    annotation_df <- data.frame(vital_status = exp_bf_filt_inter_sl$vital_status,
                                cluster = exp_bf_filt_inter_sl$cluster)
    rownames(annotation_df) <- rownames(exp_bf_filt_inter_sl)
    
    # Create a named color vector for the unique values of vital_status
    num_vital_status = brewer.pal(3, "Spectral")[c(1,3)]
    col_vital_status = setNames(num_vital_status, unique(exp_bf_filt_inter_sl$vital_status))
    
    cluster_colors <- c("short" = "yellow", "long" = "black")
    subtypes_colors = list(vital_status = col_vital_status , cluster = cluster_colors)
    
    png(filename = paste0(CancerType,"_exp_sa_features_w_short_long_clusterT_surv.png"),
        width = 35, height = 35,  units = "cm" ,pointsize = 12,
        bg = "white", res = 1200, family = "")
    
    tmp = pheatmap::pheatmap(as.matrix(t(exp_bf_filt_inter_sl[,which(!colnames(exp_bf_filt_inter_sl) %in% c("vital_status","cluster"))])),
                             annotation_col = annotation_df,
                             annotation_colors = subtypes_colors,
                             cluster_cols = T)
    print(tmp)
    dev.off()
    
    png(filename = paste0(CancerType,"_exp_sa_features_w_short_long_clusterF_surv.png"),
        width = 35, height = 35,  units = "cm" ,pointsize = 12,
        bg = "white", res = 1200, family = "")
    
    tmp2 = pheatmap::pheatmap(as.matrix(t(exp_bf_filt_inter_sl[,which(!colnames(exp_bf_filt_inter_sl) %in% c("vital_status","cluster"))])),
                              annotation_col = annotation_df,
                              annotation_colors = subtypes_colors,
                              cluster_cols = F)
    print(tmp2)
    
    dev.off()
    
    png(filename = paste0(CancerType,"_exp_sa_features_w_short_long_complex_surv.png"),
        width = 35, height = 35,  units = "cm" ,pointsize = 12,
        bg = "white", res = 1200, family = "")
    tmp3 = ComplexHeatmap::pheatmap(as.matrix(t(exp_bf_filt_inter_sl[,which(!colnames(exp_bf_filt_inter_sl) %in% c("vital_status","cluster"))])),
                                    column_split = annotation_df$vital_status,
                                    annotation_col = annotation_df,
                                    annotation_colors = subtypes_colors,
                                    cluster_cols = T)
    
    print(tmp3)
    dev.off()
  }

  # mut
  
  if (sum(colSums(mut_bf_filt_inter_sl[,which(!colnames(mut_bf_filt_inter_sl) %in% c("vital_status","cluster"))]) == 0) == 
      ncol(mut_bf_filt_inter_sl[,which(!colnames(mut_bf_filt_inter_sl) %in% c("vital_status","cluster"))])) {
    next
  } else {
    
    annotation_df <- data.frame(vital_status = mut_bf_filt_inter_sl$vital_status,
                                cluster = mut_bf_filt_inter_sl$cluster)
    rownames(annotation_df) <- rownames(mut_bf_filt_inter_sl)
    
    # Create a named color vector for the unique values of vital_status
    num_vital_status = brewer.pal(3, "Spectral")[c(1,3)]
    col_vital_status = setNames(num_vital_status, unique(mut_bf_filt_inter_sl$vital_status))
    
    cluster_colors <- c("short" = "yellow", "long" = "black")
    subtypes_colors = list(vital_status = col_vital_status , cluster = cluster_colors)
    
    png(filename = paste0(CancerType,"_mut_sa_features_w_short_long_clusterT_surv.png"),
        width = 25, height = 25,  units = "cm" ,pointsize = 12,
        bg = "white", res = 1200, family = "")
    
    tmp = pheatmap::pheatmap(as.matrix(t(mut_bf_filt_inter_sl[,which(!colnames(mut_bf_filt_inter_sl) %in% c("vital_status","cluster"))])),
                             annotation_col = annotation_df,
                             annotation_colors = subtypes_colors,
                             cluster_cols = T)
    print(tmp)
    dev.off()
    
    png(filename = paste0(CancerType,"_mut_sa_features_w_short_long_clusterF_surv.png"),
        width = 35, height = 35,  units = "cm" ,pointsize = 12,
        bg = "white", res = 1200, family = "")
    
    tmp2 = pheatmap::pheatmap(as.matrix(t(mut_bf_filt_inter_sl[,which(!colnames(mut_bf_filt_inter_sl) %in% c("vital_status","cluster"))])),
                              annotation_col = annotation_df,
                              annotation_colors = subtypes_colors,
                              cluster_cols = F)
    print(tmp2)
    
    dev.off()
    
    png(filename = paste0(CancerType,"_mut_sa_features_w_short_long_complex_surv.png"),
        width = 35, height = 35,  units = "cm" ,pointsize = 12,
        bg = "white", res = 1200, family = "")
    tmp3 = ComplexHeatmap::pheatmap(as.matrix(t(mut_bf_filt_inter_sl[,which(!colnames(mut_bf_filt_inter_sl) %in% c("vital_status","cluster"))])),
                                    column_split = annotation_df$vital_status,
                                    annotation_col = annotation_df,
                                    annotation_colors = subtypes_colors,
                                    cluster_cols = T)
    
    print(tmp3)
    dev.off()
  }
    
  
}

