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
filepath = "/home/seokwon/nas/"
ref_path = paste0(filepath, "99.reference/")

Cancerlist = dir(paste0(filepath, "/00.data/filtered_TCGA/"))
sce_path = "/mnt/gluster_server/data/raw/TCGA_data/00.data/"
cli = fread("~/nas/99.reference/all_clin_indexed.csv")
setwd("~/nas/04.Results/subtypes/tmn/")

# for all
for (num_CancerType in Cancerlist) {
  
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  
  # call input
  # 
  # sce = readRDS(paste0(sce_path,num_CancerType,"/", CancerType,".sce_raw_snv_cnv.rds"))
  # sce_df = as.data.frame(colData(sce))
  
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
  
  if (CancerType == "TCGA-COADREAD") {
    cli_surv = cli[cli$project %in% c("TCGA-COAD","TCGA-READ"),
                   c("submitter_id",
                     "vital_status",
                     "days_to_death",
                     "days_to_last_follow_up",
                     "ajcc_pathologic_t",
                     "ajcc_pathologic_n",
                     "ajcc_pathologic_m")]
    
  } else if (CancerType == "TCGA-KIDNEY") {
    cli_surv = cli[cli$project %in% c("TCGA-KIRP","TCGA-KIRC","TCGA-KICH"),
                   c("submitter_id",
                     "vital_status",
                     "days_to_death",
                     "days_to_last_follow_up",
                     "ajcc_pathologic_t",
                     "ajcc_pathologic_n",
                     "ajcc_pathologic_m")]
    
  } else {
    cli_surv = cli[cli$project == CancerType,
                   c("submitter_id",
                     "vital_status",
                     "days_to_death",
                     "days_to_last_follow_up",
                     "ajcc_pathologic_t",
                     "ajcc_pathologic_n",
                     "ajcc_pathologic_m")]
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
  
  if (sum(!is.na(cli_surv$ajcc_pathologic_t)) == 0 || sum(!is.na(cli_surv$ajcc_pathologic_n)) == 0 || sum(!is.na(cli_surv$ajcc_pathologic_m)) == 0) {
    print(paste0( "there are not tmn data in ",CancerType))
    next
  }
  
  cli_surv_filtered = cli_surv[which(!is.na(cli_surv$ajcc_pathologic_t)),]
  cli_surv_filtered = cli_surv_filtered[which(!is.na(cli_surv_filtered$ajcc_pathologic_n)),]
  cli_surv_filtered = cli_surv_filtered[which(!is.na(cli_surv_filtered$ajcc_pathologic_m)),]
  
  cli_surv_filtered$ajcc_pathologic_t <- ifelse(grepl("^T1", cli_surv_filtered$ajcc_pathologic_t), "T1", cli_surv_filtered$ajcc_pathologic_t)
  cli_surv_filtered$ajcc_pathologic_t <- ifelse(grepl("^T2", cli_surv_filtered$ajcc_pathologic_t), "T2", cli_surv_filtered$ajcc_pathologic_t)
  cli_surv_filtered$ajcc_pathologic_t <- ifelse(grepl("^T3", cli_surv_filtered$ajcc_pathologic_t), "T3", cli_surv_filtered$ajcc_pathologic_t)
  cli_surv_filtered$ajcc_pathologic_t <- ifelse(grepl("^T4", cli_surv_filtered$ajcc_pathologic_t), "T4", cli_surv_filtered$ajcc_pathologic_t)
  
  cli_surv_filtered$ajcc_pathologic_n <- ifelse(grepl("^N1", cli_surv_filtered$ajcc_pathologic_n), "N1", cli_surv_filtered$ajcc_pathologic_n)
  cli_surv_filtered$ajcc_pathologic_n <- ifelse(grepl("^N2", cli_surv_filtered$ajcc_pathologic_n), "N2", cli_surv_filtered$ajcc_pathologic_n)
  cli_surv_filtered$ajcc_pathologic_n <- ifelse(grepl("^N3", cli_surv_filtered$ajcc_pathologic_n), "N3", cli_surv_filtered$ajcc_pathologic_n)
  
  cli_surv_filtered$ajcc_pathologic_m <- ifelse(grepl("^M1", cli_surv_filtered$ajcc_pathologic_m), "M1", cli_surv_filtered$ajcc_pathologic_m)
  
  
  # filter patient
  
  dup_cli = cli_surv_filtered[which(duplicated(cli_surv_filtered$submitter_id)),]
  
  for (dup_pat in dup_cli$submitter_id) {
    tmp_dup = cli_surv_filtered[which(cli_surv_filtered$submitter_id == dup_pat),]
    if (length(unique(tmp_dup$status) == 1)) {
      tmp_dup[-nrow(tmp_dup),]
      cli_surv_filtered = cli_surv_filtered[which(!cli_surv_filtered$submitter_id == tmp_dup$submitter_id[-nrow(tmp_dup)]),]
    }
  }
  cli_surv_filtered = as.data.frame(cli_surv_filtered)
  rownames(cli_surv_filtered) = cli_surv_filtered$submitter_id
  
  # mapping tmn
  # dual
  dual_bf_filt_inter = dual_bf_filt[intersect( rownames(cli_surv_filtered) , rownames(dual_bf_filt)),]
  cli_dual_filt = cli_surv_filtered[intersect( rownames(cli_surv_filtered) , rownames(dual_bf_filt)),]
  
  if (all.equal(rownames(dual_bf_filt_inter) , rownames(cli_dual_filt)) ) {

    dual_bf_filt_inter$pat_t = cli_dual_filt$ajcc_pathologic_t 
    dual_bf_filt_inter$pat_n = cli_dual_filt$ajcc_pathologic_n 
    dual_bf_filt_inter$pat_m = cli_dual_filt$ajcc_pathologic_m 
  }
  
  # exp
  exp_bf_filt_inter = exp_bf_filt[intersect( rownames(cli_surv_filtered) , rownames(exp_bf_filt)),]
  cli_exp_filt = cli_surv_filtered[intersect( rownames(cli_surv_filtered) , rownames(exp_bf_filt)),]
  
  if (all.equal(rownames(exp_bf_filt_inter) , rownames(cli_exp_filt)) ) {
    exp_bf_filt_inter$pat_t = cli_exp_filt$ajcc_pathologic_t 
    exp_bf_filt_inter$pat_n = cli_exp_filt$ajcc_pathologic_n 
    exp_bf_filt_inter$pat_m = cli_exp_filt$ajcc_pathologic_m 
  }
  
  # mut
  mut_bf_filt_inter = mut_bf_filt[intersect( rownames(cli_surv_filtered) , rownames(mut_bf_filt)),]
  cli_mut_filt = cli_surv_filtered[intersect( rownames(cli_surv_filtered) , rownames(mut_bf_filt)),]
  
  if (all.equal(rownames(mut_bf_filt_inter) , rownames(cli_mut_filt)) ) {
    mut_bf_filt_inter$pat_t = cli_mut_filt$ajcc_pathologic_t 
    mut_bf_filt_inter$pat_n = cli_mut_filt$ajcc_pathologic_n 
    mut_bf_filt_inter$pat_m = cli_mut_filt$ajcc_pathologic_m 
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
  
  dual_bf_filt_inter_sl = dual_bf_filt_inter_sl[which(!is.na(dual_bf_filt_inter_sl$pat_t)),]
  dual_bf_filt_inter_sl = dual_bf_filt_inter_sl[which(!is.na(dual_bf_filt_inter_sl$pat_n)),]
  dual_bf_filt_inter_sl = dual_bf_filt_inter_sl[which(!is.na(dual_bf_filt_inter_sl$pat_m)),]
  
  exp_bf_filt_inter_sl = exp_bf_filt_inter_sl[which(!is.na(exp_bf_filt_inter_sl$pat_t)),]
  exp_bf_filt_inter_sl = exp_bf_filt_inter_sl[which(!is.na(exp_bf_filt_inter_sl$pat_n)),]
  exp_bf_filt_inter_sl = exp_bf_filt_inter_sl[which(!is.na(exp_bf_filt_inter_sl$pat_m)),]
  
  mut_bf_filt_inter_sl = mut_bf_filt_inter_sl[which(!is.na(mut_bf_filt_inter_sl$pat_t)),]
  mut_bf_filt_inter_sl = mut_bf_filt_inter_sl[which(!is.na(mut_bf_filt_inter_sl$pat_n)),]
  mut_bf_filt_inter_sl = mut_bf_filt_inter_sl[which(!is.na(mut_bf_filt_inter_sl$pat_m)),]
  
  dual_bf_filt_inter_sl = dual_bf_filt_inter_sl[order(dual_bf_filt_inter_sl$cluster),]
  exp_bf_filt_inter_sl = exp_bf_filt_inter_sl[order(exp_bf_filt_inter_sl$cluster),]
  mut_bf_filt_inter_sl = mut_bf_filt_inter_sl[order(mut_bf_filt_inter_sl$cluster),]
  
  # pic
  # dual
  annotation_df <- data.frame(pat_t = dual_bf_filt_inter_sl$pat_t,
                              pat_n = dual_bf_filt_inter_sl$pat_n,
                              pat_m = dual_bf_filt_inter_sl$pat_m,
                              cluster = dual_bf_filt_inter_sl$cluster)
  rownames(annotation_df) <- rownames(dual_bf_filt_inter_sl)
  
  # Create a named color vector for the unique values of vital_status
  
  num_pat_t = brewer.pal(length(unique(dual_bf_filt_inter_sl$pat_t)), "RdGy")
  col_pat_t = setNames(num_pat_t, unique(dual_bf_filt_inter_sl$pat_t))
  
  num_pat_n = brewer.pal(length(unique(dual_bf_filt_inter_sl$pat_n)), "Spectral")
  col_pat_n = setNames(num_pat_n, unique(dual_bf_filt_inter_sl$pat_n))
  
  num_pat_m = brewer.pal(length(unique(dual_bf_filt_inter_sl$pat_m)), "Set3")
  col_pat_m = setNames(num_pat_m, unique(dual_bf_filt_inter_sl$pat_m))
  
  cluster_colors <- c("short" = "yellow", "long" = "black")
  subtypes_colors = list(pat_t = col_pat_t , pat_n = col_pat_n, pat_m = col_pat_m , cluster = cluster_colors)
  
  png(filename = paste0(CancerType,"_dual_features_w_short_long_clusterT_tmn.png"),
      width = 35, height = 35,  units = "cm" ,pointsize = 12,
      bg = "white", res = 1200, family = "")
  
  tmp = pheatmap::pheatmap(as.matrix(t(dual_bf_filt_inter_sl[,which(!colnames(dual_bf_filt_inter_sl) %in% c("pat_t","pat_n","pat_m","cluster"))])),
                           annotation_col = annotation_df,
                           annotation_colors = subtypes_colors,
                           cluster_cols = T)
  print(tmp)
  dev.off()
  
  png(filename = paste0(CancerType,"_dual_features_w_short_long_clusterF_tmn.png"),
      width = 35, height = 35,  units = "cm" ,pointsize = 12,
      bg = "white", res = 1200, family = "")
  
  tmp2 = pheatmap::pheatmap(as.matrix(t(dual_bf_filt_inter_sl[,which(!colnames(dual_bf_filt_inter_sl) %in% c("pat_t","pat_n","pat_m","cluster"))])),
                            annotation_col = annotation_df,
                            annotation_colors = subtypes_colors,
                            cluster_cols = F)
  print(tmp2)
  
  dev.off()
  
  png(filename = paste0(CancerType,"_dual_features_w_short_long_complex_t.png"),
      width = 35, height = 35,  units = "cm" ,pointsize = 12,
      bg = "white", res = 1200, family = "")
  tmp3 = ComplexHeatmap::pheatmap(as.matrix(t(dual_bf_filt_inter_sl[,which(!colnames(dual_bf_filt_inter_sl) %in% c("pat_t","pat_n","pat_m","cluster"))])),
                                  column_split = annotation_df$pat_t,
                                  annotation_col = annotation_df,
                                  annotation_colors = subtypes_colors,
                                  cluster_cols = T)
  
  print(tmp3)
  dev.off()
  
  png(filename = paste0(CancerType,"_dual_features_w_short_long_complex_n.png"),
      width = 35, height = 35,  units = "cm" ,pointsize = 12,
      bg = "white", res = 1200, family = "")
  tmp4 = ComplexHeatmap::pheatmap(as.matrix(t(dual_bf_filt_inter_sl[,which(!colnames(dual_bf_filt_inter_sl) %in% c("pat_t","pat_n","pat_m","cluster"))])),
                                  column_split = annotation_df$pat_n,
                                  annotation_col = annotation_df,
                                  annotation_colors = subtypes_colors,
                                  cluster_cols = T)
  
  print(tmp4)
  dev.off()
  
  png(filename = paste0(CancerType,"_dual_features_w_short_long_complex_m.png"),
      width = 35, height = 35,  units = "cm" ,pointsize = 12,
      bg = "white", res = 1200, family = "")
  tmp5 = ComplexHeatmap::pheatmap(as.matrix(t(dual_bf_filt_inter_sl[,which(!colnames(dual_bf_filt_inter_sl) %in% c("pat_t","pat_n","pat_m","cluster"))])),
                                  column_split = annotation_df$pat_m,
                                  annotation_col = annotation_df,
                                  annotation_colors = subtypes_colors,
                                  cluster_cols = T)
  
  print(tmp5)
  dev.off()
  
  # exp
  
  if (sum(colSums(exp_bf_filt_inter_sl[,which(!colnames(exp_bf_filt_inter_sl) %in%c("pat_t","pat_n","pat_m","cluster"))]) == 0) == 
      ncol(exp_bf_filt_inter_sl[,which(!colnames(exp_bf_filt_inter_sl) %in% c("pat_t","pat_n","pat_m","cluster"))])) {
    next
    
    } else {
    
      annotation_df <- data.frame(pat_t = exp_bf_filt_inter_sl$pat_t,
                                  pat_n = exp_bf_filt_inter_sl$pat_n,
                                  pat_m = exp_bf_filt_inter_sl$pat_m,
                                  cluster = exp_bf_filt_inter_sl$cluster)
      rownames(annotation_df) <- rownames(exp_bf_filt_inter_sl)
      
      # Create a named color vector for the unique values of vital_status
      
      num_pat_t = brewer.pal(length(unique(exp_bf_filt_inter_sl$pat_t)), "RdGy")
      col_pat_t = setNames(num_pat_t, unique(exp_bf_filt_inter_sl$pat_t))
      
      num_pat_n = brewer.pal(length(unique(exp_bf_filt_inter_sl$pat_n)), "Spectral")
      col_pat_n = setNames(num_pat_n, unique(exp_bf_filt_inter_sl$pat_n))
      
      num_pat_m = brewer.pal(length(unique(exp_bf_filt_inter_sl$pat_m)), "Set3")
      col_pat_m = setNames(num_pat_m, unique(exp_bf_filt_inter_sl$pat_m))
      
      cluster_colors <- c("short" = "yellow", "long" = "black")
      subtypes_colors = list(pat_t = col_pat_t , pat_n = col_pat_n, pat_m = col_pat_m , cluster = cluster_colors)
      
      png(filename = paste0(CancerType,"_exp_features_w_short_long_clusterT_tmn.png"),
          width = 35, height = 35,  units = "cm" ,pointsize = 12,
          bg = "white", res = 1200, family = "")
      
      tmp = pheatmap::pheatmap(as.matrix(t(exp_bf_filt_inter_sl[,which(!colnames(exp_bf_filt_inter_sl) %in% c("pat_t","pat_n","pat_m","cluster"))])),
                               annotation_col = annotation_df,
                               annotation_colors = subtypes_colors,
                               cluster_cols = T)
      print(tmp)
      dev.off()
      
      png(filename = paste0(CancerType,"_exp_features_w_short_long_clusterF_tmn.png"),
          width = 35, height = 35,  units = "cm" ,pointsize = 12,
          bg = "white", res = 1200, family = "")
      
      tmp2 = pheatmap::pheatmap(as.matrix(t(exp_bf_filt_inter_sl[,which(!colnames(exp_bf_filt_inter_sl) %in% c("pat_t","pat_n","pat_m","cluster"))])),
                                annotation_col = annotation_df,
                                annotation_colors = subtypes_colors,
                                cluster_cols = F)
      print(tmp2)
      
      dev.off()
      
      png(filename = paste0(CancerType,"_exp_features_w_short_long_complex_t.png"),
          width = 35, height = 35,  units = "cm" ,pointsize = 12,
          bg = "white", res = 1200, family = "")
      tmp3 = ComplexHeatmap::pheatmap(as.matrix(t(exp_bf_filt_inter_sl[,which(!colnames(exp_bf_filt_inter_sl) %in% c("pat_t","pat_n","pat_m","cluster"))])),
                                      column_split = annotation_df$pat_t,
                                      annotation_col = annotation_df,
                                      annotation_colors = subtypes_colors,
                                      cluster_cols = T)
      
      print(tmp3)
      dev.off()
      
      png(filename = paste0(CancerType,"_exp_features_w_short_long_complex_n.png"),
          width = 35, height = 35,  units = "cm" ,pointsize = 12,
          bg = "white", res = 1200, family = "")
      tmp4 = ComplexHeatmap::pheatmap(as.matrix(t(exp_bf_filt_inter_sl[,which(!colnames(exp_bf_filt_inter_sl) %in% c("pat_t","pat_n","pat_m","cluster"))])),
                                      column_split = annotation_df$pat_n,
                                      annotation_col = annotation_df,
                                      annotation_colors = subtypes_colors,
                                      cluster_cols = T)
      
      print(tmp4)
      dev.off()
      
      png(filename = paste0(CancerType,"_exp_features_w_short_long_complex_m.png"),
          width = 35, height = 35,  units = "cm" ,pointsize = 12,
          bg = "white", res = 1200, family = "")
      tmp5 = ComplexHeatmap::pheatmap(as.matrix(t(exp_bf_filt_inter_sl[,which(!colnames(exp_bf_filt_inter_sl) %in% c("pat_t","pat_n","pat_m","cluster"))])),
                                      column_split = annotation_df$pat_m,
                                      annotation_col = annotation_df,
                                      annotation_colors = subtypes_colors,
                                      cluster_cols = T)
      
      print(tmp5)
      dev.off()
    }
  
  if (sum(colSums(mut_bf_filt_inter_sl[,which(!colnames(mut_bf_filt_inter_sl) %in%c("pat_t","pat_n","pat_m","cluster"))]) == 0) == 
      ncol(mut_bf_filt_inter_sl[,which(!colnames(mut_bf_filt_inter_sl) %in% c("pat_t","pat_n","pat_m","cluster"))])) {
    next
    
    } else {
      
      # mut
      annotation_df <- data.frame(pat_t = mut_bf_filt_inter_sl$pat_t,
                                  pat_n = mut_bf_filt_inter_sl$pat_n,
                                  pat_m = mut_bf_filt_inter_sl$pat_m,
                                  cluster = mut_bf_filt_inter_sl$cluster)
      rownames(annotation_df) <- rownames(mut_bf_filt_inter_sl)
      
      # Create a named color vector for the unique values of vital_status
      
      num_pat_t = brewer.pal(length(unique(mut_bf_filt_inter_sl$pat_t)), "RdGy")
      col_pat_t = setNames(num_pat_t, unique(mut_bf_filt_inter_sl$pat_t))
      
      num_pat_n = brewer.pal(length(unique(mut_bf_filt_inter_sl$pat_n)), "Spectral")
      col_pat_n = setNames(num_pat_n, unique(mut_bf_filt_inter_sl$pat_n))
      
      num_pat_m = brewer.pal(length(unique(mut_bf_filt_inter_sl$pat_m)), "Set3")
      col_pat_m = setNames(num_pat_m, unique(mut_bf_filt_inter_sl$pat_m))
      
      cluster_colors <- c("short" = "yellow", "long" = "black")
      subtypes_colors = list(pat_t = col_pat_t , pat_n = col_pat_n, pat_m = col_pat_m , cluster = cluster_colors)
      
      png(filename = paste0(CancerType,"_mut_features_w_short_long_clusterT_surv.png"),
          width = 25, height = 25,  units = "cm" ,pointsize = 12,
          bg = "white", res = 1200, family = "")
      
      tmp = pheatmap::pheatmap(as.matrix(t(mut_bf_filt_inter_sl[,which(!colnames(mut_bf_filt_inter_sl) %in% c("pat_t","pat_n","pat_m","cluster"))])),
                               annotation_col = annotation_df,
                               annotation_colors = subtypes_colors,
                               cluster_cols = T)
      print(tmp)
      dev.off()
      
      png(filename = paste0(CancerType,"_mut_features_w_short_long_clusterF_surv.png"),
          width = 35, height = 35,  units = "cm" ,pointsize = 12,
          bg = "white", res = 1200, family = "")
      
      tmp2 = pheatmap::pheatmap(as.matrix(t(mut_bf_filt_inter_sl[,which(!colnames(mut_bf_filt_inter_sl) %in% c("pat_t","pat_n","pat_m","cluster"))])),
                                annotation_col = annotation_df,
                                annotation_colors = subtypes_colors,
                                cluster_cols = F)
      print(tmp2)
      
      dev.off()
      
      png(filename = paste0(CancerType,"_mut_features_w_short_long_complex_surv.png"),
          width = 35, height = 35,  units = "cm" ,pointsize = 12,
          bg = "white", res = 1200, family = "")
      tmp3 = ComplexHeatmap::pheatmap(as.matrix(t(mut_bf_filt_inter_sl[,which(!colnames(mut_bf_filt_inter_sl) %in% c("pat_t","pat_n","pat_m","cluster"))])),
                                      column_split = annotation_df$vital_status,
                                      annotation_col = annotation_df,
                                      annotation_colors = subtypes_colors,
                                      cluster_cols = T)
      
      print(tmp3)
      dev.off()
    }
}

