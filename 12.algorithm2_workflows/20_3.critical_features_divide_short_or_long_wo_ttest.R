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
surv_total_results = read_xlsx("~/nas/04.Results/Total_results_survpval2.xlsx")

setwd("~/nas/04.Results/short_long/")
# num_CancerType = "19.TCGA-LIHC"
# for all
for (num_CancerType in Cancerlist) {
  
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  
  # call input
  
  cancer_bf = read.csv(paste0(filepath,"04.Results/bestfeatures/",CancerType, "_best_features.csv"))
  duration_log_df = readRDS(paste0(main.path_tc, "/", CancerType,"_dual_add_duration_log.rds"))
  
  # cut the number of best features 
  cancer_bf_cut = cancer_bf[1:surv_total_results[which(surv_total_results$CancerType == CancerType),]$num_of_features,]
  best_features_df = duration_log_df[,cancer_bf_cut$variable]
  
  if (all.equal(rownames(best_features_df), rownames(duration_log_df))) {
    
    best_features_df$vitalstatus = duration_log_df$vitalstatus
    best_features_df$duration = duration_log_df$duration
    best_features_df = best_features_df[which(!is.na(best_features_df$duration)),]
    
    best_features_df$status = NA
    best_features_df$status[which(best_features_df$vitalstatus == "Dead")] = 1
    best_features_df$status[which(best_features_df$vitalstatus == "Alive")] = 0
    
  }
  
  data_bf = best_features_df
  # annotation cluster 1 or 2 (long or short) by best pval score
  # 18 means best p-value when devide two cluster
  
  out = pheatmap::pheatmap((data_bf[,which(!colnames(data_bf) %in% c("vitalstatus","duration","status"))] > -log(0.05))*1 ,
                           cluster_cols = T,
                           cluster_rows = T,
                           labels_cols = "", 
                           show_rownames = T,
                           silent = T)
  
  tmp_pheat_cut = as.data.frame (cutree(out$tree_row, 2) , out[["tree_row"]][["labels"]])
  colnames(tmp_pheat_cut) = "cluster"
  
  if (all.equal(rownames(data_bf), rownames(tmp_pheat_cut))) {
    data_bf$cluster = tmp_pheat_cut$cluster 
  }
  
  if (mean(data_bf[which(data_bf$cluster == 1),]$duration) <
      mean(data_bf[which(data_bf$cluster == 2),]$duration)) {
    data_bf[which(data_bf$cluster == 1),]$cluster = 3
    data_bf[which(data_bf$cluster == 2),]$cluster = 1
    data_bf[which(data_bf$cluster == 3),]$cluster = 2
  } else {
    data_bf = data_bf
  }
  
  fit = survfit(Surv(duration, status) ~ cluster, data = data_bf)
  # ggsurvplot(fit, data = data_bf, risk.table = TRUE,
  #            palette = "jco", pval = TRUE, surv.median.line = "hv", xlab = "days")
  
  # hypothesis : cluster 1 = better prognosis
  if (mean(fit$surv[1:fit[['strata']][['cluster=1']]]) > mean(fit$surv[fit[['strata']][['cluster=1']] + 1: fit[['strata']][['cluster=2']]])) {
    short_group = data_bf[which(data_bf$cluster == 2),]
    long_group = data_bf[which(data_bf$cluster == 1),]
    
    # short_group_for_fig = data_bf[which(data_bf$cluster == 2),]
    # long_group_for_fig = data_bf[which(data_bf$cluster == 1),]
    
  } else if (mean(fit$surv[1:fit[['strata']][['cluster=1']]]) < mean(fit$surv[fit[['strata']][['cluster=1']] + 1: fit[['strata']][['cluster=2']]])) {
    short_group = data_bf[which(data_bf$cluster == 1),]
    long_group = data_bf[which(data_bf$cluster == 2),]
    
    # short_group_for_fig = data_bf[which(data_bf$cluster == 1),]
    # long_group_for_fig = data_bf[which(data_bf$cluster == 2),]
  } else {
    print("I don't know")
  }
  
  short_group = na.omit(short_group)
  long_group = na.omit(long_group)
  
  # short_group_for_fig = na.omit(short_group_for_fig)
  # long_group_for_fig = na.omit(long_group_for_fig)
  
  # for ttest method (short vs long)
  short_group$cluster = "short" 
  long_group$cluster = "long" 
  
  total_group = rbind(long_group,short_group)
  
  saveRDS(total_group, paste0("~/nas/04.Results/short_long/",CancerType,"_critical_features_short_long.rds"))
  
  short_cluster_path = c()
  long_cluster_path = c()
  
  deg_short_long = colnames(total_group)[!colnames(total_group) %in% c("vitalstatus","duration","status","cluster")]
  deg_group = total_group[,c(deg_short_long,"cluster")]
  
  for (deg_path in deg_short_long) {
    if (mean(deg_group[which(deg_group$cluster == "short"), deg_path]) > mean(deg_group[which(deg_group$cluster == "long"), deg_path])) {
      short_cluster_path <- c(short_cluster_path, deg_path)
    } else {
      long_cluster_path <- c(long_cluster_path, deg_path)
    }
  }
  
  if (length(short_cluster_path) != 0 && length(long_cluster_path) != 0) {
    print("There are well divided as short or long pathway")
  }
  
  merge_short_long = c(short_cluster_path,long_cluster_path)
  
  cancer_short_long = cancer_bf[which(cancer_bf$variable %in% merge_short_long),]
  rownames(cancer_short_long) = NULL
  cancer_short_long$relative_importance =NULL
  cancer_short_long$scaled_importance = NULL
  cancer_short_long$percentage = NULL
  cancer_short_long$classification = NA
  
  if (length(short_cluster_path) != 0 && length(long_cluster_path) != 0) {
    cancer_short_long[which(cancer_short_long$variable %in% short_cluster_path),]$classification = "short"
    cancer_short_long[which(cancer_short_long$variable %in% long_cluster_path),]$classification = "long"
  } else if (length(short_cluster_path) == 0 ) {
    cancer_short_long[which(cancer_short_long$variable %in% long_cluster_path),]$classification = "long"
  } else {
    cancer_short_long[which(cancer_short_long$variable %in% short_cluster_path),]$classification = "short"
  }
  
  write.xlsx(cancer_short_long , paste0("~/nas/04.Results/short_long/",CancerType,"_critical_features_short_long_wo_ttest.xlsx"))
  
  
}  
