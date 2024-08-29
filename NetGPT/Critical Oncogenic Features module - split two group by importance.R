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
library(readxl)

set.seed(13524)

total_results_pval = data.frame()

CancerType = "Intered your cancer type"

# call input
annotate_best_features =  "your file : best_features.csv" # come from the output of extract important score.R
duration_log_df = "your files" # you need the information of patients period data (suvial time). data/TCGA-XXX_dual_with_duration.rds

best_features_df = duration_log_df[,annotate_best_features$variable]
best_features_df$vitalstatus = duration_log_df$vitalstatus
best_features_df$duration = duration_log_df$duration
best_features_df = best_features_df[which(!is.na(best_features_df$duration)),]

best_features_df$status = NA
best_features_df$status[which(best_features_df$vitalstatus == "Dead")] = 1
best_features_df$status[which(best_features_df$vitalstatus == "Alive")] = 0

# results_pval
result_surv_pval = data.frame(matrix(ncol = 1))
colnames(result_surv_pval) = "pval"

result_surv_pval$CancerType = NA
result_surv_pval$num_of_features = NA
result_surv_pval$diff_balance_ratio = NA

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
  
  if (mean(best_features_df[which(best_features_df$cluster == 1),]$duration) <
      mean(best_features_df[which(best_features_df$cluster == 2),]$duration)) {
    best_features_df[which(best_features_df$cluster == 1),]$cluster = 3
    best_features_df[which(best_features_df$cluster == 2),]$cluster = 1
    best_features_df[which(best_features_df$cluster == 3),]$cluster = 2
  } else {
    best_features_df = best_features_df
  }
  
  fit = survfit(Surv(duration, status) ~ cluster, data = best_features_df)
  
  # ### if you wana print the plot, plz proceed below annotated code
  # png(filename = paste0(CancerType, "_",last_num,"_survival_plot_adjust.png"),
  #     width = 35, height = 35,  units = "cm" ,pointsize = 12,
  #     bg = "white", res = 1200, family = "")
  # 
  # tmp_surv = ggsurvplot(
  #   fit,    # survfit object with calculated statistics.
  #   data = best_features_df,
  #   pval = TRUE,             # show p-value of log-rank test.
  #   conf.int = TRUE,         # show confidence intervals for 
  #   ggtheme = RTCGA::theme_RTCGA(),
  #   risk.table = TRUE,
  #   xlab = "days",
  #   surv.median.line = "hv"
  # )
  # 
  # print(tmp_surv)
  # 
  # dev.off()
  # ###
  
  result_surv_pval[last_num,"num_of_features"] = last_num
  result_surv_pval[last_num,"pval"] = surv_pvalue(fit)$pval
  result_surv_pval[last_num,"CancerType"] = CancerType
  
  if (round(fit$n[1] / (fit$n[1] + fit$n[2]), 2) > round(fit$n[2] / (fit$n[1] + fit$n[2]), 2) ) {
    result_surv_pval[last_num,"diff_balance_ratio"] = round(fit$n[1] / (fit$n[1] + fit$n[2]), 2) - round(fit$n[2] / (fit$n[1] + fit$n[2]), 2)
  } else {
    result_surv_pval[last_num,"diff_balance_ratio"] = round(fit$n[2] / (fit$n[1] + fit$n[2]), 2) - round(fit$n[1] / (fit$n[1] + fit$n[2]), 2) 
  }
  
  if (round(fit$n[1] / (fit$n[1] + fit$n[2]), 2) <= 0.7 && round(fit$n[1] / (fit$n[1] + fit$n[2]), 2) >= 0.3) {
    result_surv_pval[last_num,"bias"] = "balanced"
  } else {
    result_surv_pval[last_num,"bias"] = "not balanced"
  }
  
  remove(tmp_pheat_cut,fit,plot_roc)
}

# save the result
write.csv(result_surv_pval, paste0(CancerType, "_rawdata_survpval.csv"))



