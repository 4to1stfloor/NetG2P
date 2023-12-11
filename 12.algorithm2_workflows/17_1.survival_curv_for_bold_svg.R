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

filepath = "/home/seokwon/nas/"
ref_path = paste0(filepath, "99.reference/")
Cancerlist = dir(paste0(filepath, "/00.data/filtered_TCGA/"))

# for fic
fig_path = "~/nas/04.Results/survival/"
if(!dir.exists(fig_path)){
  dir.create(fig_path)
  print(paste0("Created folder: ", fig_path))
} else {
  print(paste0("Folder already exists: ", fig_path))
}
setwd(fig_path)
set.seed(13524)

num_CancerType = "29.TCGA-LGG"
total_surv_result = read_xlsx("~/nas/04.Results/Total_results_survpval2.xlsx")
total_results_pval = data.frame()

for (num_CancerType in Cancerlist) {
  
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  name_CancerType = gsub('TCGA-','',CancerType)
  
  # call input
  annotate_best_features =  read.csv(paste0(filepath,"13.analysis/",CancerType,"_best_features.csv"))
  duration_log_df = readRDS(paste0(main.path_tc, "/", CancerType,"_dual_add_duration_log.rds"))
  
  # results_pval
  result_surv_pval = data.frame(matrix(ncol = 1))
  colnames(result_surv_pval) = "pval"
  
  best_features_df = duration_log_df[,annotate_best_features$variable]
  best_features_df$vitalstatus = duration_log_df$vitalstatus
  best_features_df$duration = duration_log_df$duration
  best_features_df = best_features_df[which(!is.na(best_features_df$duration)),]
  
  best_features_df$status = NA
  best_features_df$status[which(best_features_df$vitalstatus == "Dead")] = 1
  best_features_df$status[which(best_features_df$vitalstatus == "Alive")] = 0

  last_num = total_surv_result[which(total_surv_result$CancerType == CancerType),]$num_of_features
  
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
  
  median(best_features_df[which(best_features_df$cluster == 1),]$duration) <
    median(best_features_df[which(best_features_df$cluster == 2),]$duration)
  
  if (mean(best_features_df[which(best_features_df$cluster == 1),]$duration) <
      mean(best_features_df[which(best_features_df$cluster == 2),]$duration)) {
    if (CancerType == "TCGA-BRCA") {
      best_features_df = best_features_df
    } else {
      best_features_df[which(best_features_df$cluster == 1),]$cluster = 3
      best_features_df[which(best_features_df$cluster == 2),]$cluster = 1
      best_features_df[which(best_features_df$cluster == 3),]$cluster = 2
    }
  } else {
    best_features_df = best_features_df
  }
  
  fit = survfit(Surv(duration, status) ~ cluster, data = best_features_df)
  
  custom = theme_bw()+ theme(plot.title = element_text(hjust = 0.99, face = "bold"))
  
  tmp_surv = ggsurvplot(
    fit,    # survfit object with calculated statistics.
    data = best_features_df,
    pval = F,             # show p-value of log-rank test.
    # conf.int = TRUE,         # show confidence intervals for
    # ggtheme = RTCGA::theme_RTCGA(base_family = "Arial", base_size = 30),
    ggtheme = custom,
    palette = c("#86AA00","#FF9E29"),
    # title = CancerType,
    risk.table = FALSE,
    xlab = "days",
    surv.median.line = "hv",
    font.main = c(50, "bold"),
    font.x = c(50, "bold"),
    font.y = c(50, "bold"),
    font.tickslab = c(50, "bold","black"),
    font.caption = c(50, "bold"),
    font.legend = c(50, "bold"),
    legend.title = "",
    legend.labs = c("Long-term", "Short-term"),
    pval.size = 30
    
    # ncensor.plot = TRUE
    # pval.coord =c(x,y) 
  )
  # library(graphics)
  
  if (max(tmp_surv$data.survtable$time) == 6000) {
    fine_tune_x = 5500 
  } else if (max(tmp_surv$data.survtable$time) == 4000) {
    fine_tune_x = 3700 
  } else if (max(tmp_surv$data.survtable$time) == 8000) {
    fine_tune_x = 7200
  } else if (max(tmp_surv$data.survtable$time) == 5000) {
    fine_tune_x = 4500
  } else {
    fine_tune_x = max(tmp_surv$data.survtable$time) - 1000
  }
  
  if (CancerType == "TCGA-COADREAD") {
    fine_tune_x = 3300 
  }
  
  if (CancerType == "TCGA-LGG") {
    pval_label = "p < 0.0001"
  } else {
    pval_label = paste0("p = ",round(surv_pvalue(fit)$pval, 3))
  }
  
  tmp_surv$plot <- tmp_surv$plot +
    annotate("text", 
             x = 0.9, 
             y = 0.1,
             hjust = 0.05, 
             size = 20, 
             label = pval_label,
             fontface = "bold")+
    annotate("text", 
             x = fine_tune_x, 
             y = 1,
             # hjust = 0.01, 
             # vjust = 0.6,
             size = 20, 
             label = name_CancerType,
             fontface = "bold")
  
  ggsave(file = paste0(CancerType, "_",last_num,"_survival_plot_adjust_wo_conf.svg"), tmp_surv$plot, width=17, height=15, device = svg)
  
}

