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
                          palette = "jco", pval = TRUE, surv.median.line = "hv", xlab = "days")
    print(plot_roc)
    dev.off()
    
    
    result_surv_pval[last_num,"pval"] = surv_pvalue(fit)$pval
    remove(tmp_pheat_cut,fit,plot_roc)
  }
  remove(best_features_df,annotate_best_features,duration_log_df)
  write.csv(result_surv_pval, paste0(CancerType, "_result_survpval.csv"))
  
}


annotate_best_features
duration_log_df

# combine all of pathwaylink and pathway (698)

best_features = annotate_best_features$variable
other_features = colnames(duration_log_df[,!names(duration_log_df) %in% c(annotate_best_features$variable,"vitalstatus","duration")])

relative_weight = data.frame()

relative_weight[best_features,"relative_importance"] = annotate_best_features$relative_importance
relative_weight[other_features,"relative_importance"] = 10^-4
relative_weight$relative_importance = ifelse(relative_weight$relative_importance == 0, 10^-4,relative_weight$relative_importance )

just_matrix = duration_log_df[,!names(duration_log_df) %in% c("vitalstatus","duration")]

df_ordered = as.data.frame(relative_weight[colnames(just_matrix),],row.names =colnames(just_matrix),optional = T)

df_ordered_mat = as.matrix(df_ordered)
just_matrix_mat = as.matrix(just_matrix)
combine_pl = as.data.frame(just_matrix_mat%*%df_ordered_mat)
colnames(combine_pl) = "weight"

if (all.equal(rownames(combine_pl),rownames(duration_log_df))) {
  
  combine_pl$status = ifelse(duration_log_df$vitalstatus == "Alive", 0, 1)
  combine_pl$duration = duration_log_df$duration
}


combine_pl$cluster = ifelse(combine_pl$weight > mean(combine_pl$weight), 1 , 0)

fit = survfit(Surv(duration, status) ~ cluster, data = combine_pl)

ggsurvplot(fit, data = combine_pl, risk.table = TRUE,
           palette = "jco", pval = TRUE, surv.median.line = "hv", xlab = "days")

# It just bestfeatures matrix
best_features = annotate_best_features$variable
df_best_features = duration_log_df[,best_features]

df_best_features_mat = as.matrix(df_best_features)
relative_weight_df = as.data.frame(annotate_best_features$relative_importance,row.names = annotate_best_features$variable)
colnames(relative_weight_df) = "weight"
relative_weight_df$weight = ifelse(relative_weight_df$weight == 0, 10^-4,relative_weight_df$weight )
relative_weight_mat = as.matrix(relative_weight_df)

combine_pl_best_features = as.data.frame(df_best_features_mat%*%relative_weight_mat)

if (all.equal(rownames(combine_pl_best_features),rownames(duration_log_df))) {
  
  combine_pl_best_features$status = ifelse(duration_log_df$vitalstatus == "Alive", 0, 1)
  combine_pl_best_features$duration = duration_log_df$duration
}


combine_pl_best_features$cluster = ifelse(combine_pl_best_features$weight > mean(combine_pl_best_features$weight), 1 , 0)

fit2 = survfit(Surv(duration, status) ~ cluster, data = combine_pl_best_features)

ggsurvplot(fit2, data = combine_pl_best_features, risk.table = TRUE,
           palette = "jco", pval = TRUE, surv.median.line = "hv", xlab = "days")

combine_pl_best_features$weight_zscore <- scale(combine_pl_best_features$weight)

combine_pl_zscore = combine_pl_best_features[,c("weight_zscore","status")]

pheatmap(combine_pl_zscore , 
         cluster_rows = T, labels_cols = "",
         show_rownames = T)

pheatmap(combine_pl_zscore)
