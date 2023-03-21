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

filepath = "/home/seokwon/nas/"
ref_path = paste0(filepath, "99.reference/")
Cancerlist = dir(paste0(filepath, "/00.data/filtered_TCGA/"))
Cancerlist = Cancerlist[1:11]
for (num_CancerType in Cancerlist) {
  
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  
  # call input
  annotate_best_features =  read.csv(paste0(filepath,"13.analysis/",CancerType,"_best_features.csv"))
  duration_log_df = readRDS(paste0(main.path_tc, "/", CancerType,"_dual_add_duration_log.rds"))
  
  # fig_folder_create
  dir.create(paste0(main.path_tc,"/cluster_fig"))   
  fig_path = paste0(main.path_tc,"/cluster_fig")
  setwd(fig_path)
  
  # results_pval
  result_surv_pval = data.frame(matrix(ncol = 1))
  colnames(result_surv_pval) = "pval"
  
  best_features_df = duration_log_df[,annotate_best_features$variable]
  best_features_df$vitalstatus = duration_log_df$vitalstatus
  best_features_df$duration = duration_log_df$duration
  
  best_features_df$status = NA
  best_features_df$status[which(best_features_df$vitalstatus == "Dead")] = 1
  best_features_df$status[which(best_features_df$vitalstatus == "Alive")] = 0
  
  
  for (last_num in 2:length(annotate_best_features$variable)) {
    
    out = pheatmap((best_features_df[,1:last_num] > -log(0.05))*1 , cluster_cols = T,
                   cluster_rows = T, labels_cols = "",
                   show_rownames = T)
    print(paste0(CancerType , "_start : ",last_num ))
    print("----------------------------------------")
    # print(table(cutree(out$tree_row, 2)))
    
    tmp_pheat_cut = as.data.frame (cutree(out$tree_row, 2) , out[["tree_row"]][["labels"]])
    colnames(tmp_pheat_cut) = "cluster"
    best_features_df$cluster = tmp_pheat_cut$cluster
    fit = survfit(Surv(duration, status) ~ cluster, data = best_features_df)
    
    png(filename = paste0("cut",last_num ,"_cluster_test.png"),
        width = 1000, height = 1000, units = "px", pointsize = 12,
        bg = "white", res = NA, family = "")
    
    plot_roc = ggsurvplot(fit, data = best_features_df, risk.table = TRUE,
                          palette = "jco", pval = TRUE, surv.median.line = "hv", xlab = "month")
    print(plot_roc)
    dev.off()
    
    
    result_surv_pval[last_num,"pval"] = surv_pvalue(fit)$pval
    remove(tmp_pheat_cut,fit,plot_roc)
  }
  remove(best_features_df,annotate_best_features,duration_log_df)
  write.csv(result_surv_pval, paste0(CancerType, "_result_survpval.csv"))
  
}




other_features_df = test_df[,!names(df_test) %in% df_test$variable]


colSums(best_features_df[which(best_features_df$vitalstatus == "Alive"),][1:(length(best_features_df)-1)]) / nrow(best_features_df[which(best_features_df$vitalstatus == "Alive"),])
colSums(best_features_df[which(best_features_df$vitalstatus == "Dead"),][1:(length(best_features_df)-1)]) / nrow(best_features_df[which(best_features_df$vitalstatus == "Dead"),])

a = colSums(best_features_df[which(best_features_df$vitalstatus == "Alive"),][1:(length(best_features_df)-1)]) / nrow(best_features_df[which(best_features_df$vitalstatus == "Alive"),])
b = colSums(best_features_df[which(best_features_df$vitalstatus == "Dead"),][1:(length(best_features_df)-1)]) / nrow(best_features_df[which(best_features_df$vitalstatus == "Dead"),])
sum(a>b)
a>b
length(a)

colSums(other_features_df[which(other_features_df$vitalstatus == "Alive"),][1:(length(other_features_df)-4)]) >
  colSums(other_features_df[which(other_features_df$vitalstatus == "Dead"),][1:(length(other_features_df)-4)])
