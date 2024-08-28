library(survminer) 
library(ggplot2) 
library(pheatmap)
library(survival) 
library(data.table) 
library(RColorBrewer) 
library(NbClust)
library(cluster) 
library(factoextra) 
library(dplyr)
library(stringr)
library(pheatmap)
setwd("/home/seokwon/nas/04.Results/")
filepath = "/home/seokwon/nas/"
ref_path = paste0(filepath, "99.reference/")
Cancerlist = dir(paste0(filepath, "/00.data/filtered_TCGA/"))
set.seed(13524)
total_results_pval = data.frame()
# num_CancerType = "04.TCGA-CESC" 
setwd("~/nas/04.Results/surv_features_lineplot")

for (num_CancerType in Cancerlist) {
  
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  
  # call input
  annotate_best_features =  read.csv(paste0(filepath,"13.analysis/",CancerType,"_best_features.csv"))
  duration_log_df = readRDS(paste0(main.path_tc, "/", CancerType,"_dual_add_duration_log.rds"))
  
  # results_pval
  result_surv_pval = data.frame(matrix(ncol = 1))
  colnames(result_surv_pval) = "pval"
  
  best_features_df = duration_log_df[,annotate_best_features$variable]
  best_features_df$vitalstatus = duration_log_df$vitalstatus
  best_features_df$duration = duration_log_df$duration
  
  best_features_df$status = NA
  best_features_df$status[which(best_features_df$vitalstatus == "Dead")] = 1
  best_features_df$status[which(best_features_df$vitalstatus == "Alive")] = 0
  
  result_surv_pval$CancerType = NA
  result_surv_pval$num_of_features = NA
  
  for (last_num in 2:length(annotate_best_features$variable)) {
    
    out = pheatmap((best_features_df[,1:last_num] > -log(0.05))*1 , 
                   cluster_cols = T,
                   cluster_rows = T, 
                   labels_cols = "",
                   show_rownames = T,
                   silent = T)
    print(paste0(CancerType , "_start : ",last_num ))
    print("----------------------------------------")
    # print(table(cutree(out$tree_row, 2)))
    
    tmp_pheat_cut = as.data.frame (cutree(out$tree_row, 2) , out[["tree_row"]][["labels"]])
    colnames(tmp_pheat_cut) = "cluster"
    
    best_features_df$cluster = tmp_pheat_cut$cluster
    fit = survfit(Surv(duration, status) ~ cluster, data = best_features_df)
    # ggsurvplot(fit, data = best_features_df, risk.table = TRUE,
    #            palette = "jco", pval = TRUE, surv.median.line = "hv", xlab = "days")
    
    result_surv_pval[last_num,"num_of_features"] = last_num
    result_surv_pval[last_num,"pval"] = surv_pvalue(fit)$pval
    result_surv_pval[last_num,"CancerType"] = CancerType
    
    if (round(fit$n[1] / (fit$n[1] + fit$n[2]), 2) <= 0.7 && round(fit$n[1] / (fit$n[1] + fit$n[2]), 2) >= 0.3) {
      result_surv_pval[last_num,"bias"] = "balanced"
    } else {
      result_surv_pval[last_num,"bias"] = "not balanced"
    }
    
    remove(tmp_pheat_cut,fit,plot_roc)
  }
  
  result_surv_pval = na.omit(result_surv_pval)
  result_surv_pval_spe = result_surv_pval[which(result_surv_pval$pval < 0.05),]
  if (nrow(result_surv_pval_spe) < 2 ) {
    result_surv_pval_spe = result_surv_pval[order(result_surv_pval$pval),][1:5,]
  }
    
  png(filename = paste0(CancerType,"_pval_spe_features_lineplot.png"),
      width = 6000, height = 6000,  
      # units = "cm" ,pointsize = 12,
      # bg = "white",
      res = 1200, family = "")
  
  tmp_line = ggplot(data = result_surv_pval_spe, aes(x = num_of_features, y = pval)) + geom_line()
  
  print(tmp_line)
  
  dev.off()
}