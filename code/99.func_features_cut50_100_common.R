library(survival) 
library(dplyr)
library(survminer)

filepath = "/home/seokwon/nas/"
ref_path = paste0(filepath, "99.reference/")

Cancerlist = dir(paste0(filepath, "/00.data/filtered_TCGA/"))
surv_total_results = read_xlsx("~/nas/04.Results/Total_results_survpval2.xlsx")

# num_CancerType = "19.TCGA-LIHC"
# for all
for (num_CancerType in Cancerlist) {
  
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  
  # call input
  
  cancer_bf = read.csv(paste0(filepath,"04.Results/bestfeatures/",CancerType, "_critical_features.csv"))
  cancer_short_long_features = readRDS(paste0("~/nas/04.Results/short_long/",CancerType,"_critical_features_short_long.rds"))
  # cancer_bf_cut = cancer_bf[1:surv_total_results[which(surv_total_results$CancerType == CancerType),]$num_of_features,]
  duration_log_df = readRDS(paste0(main.path_tc, "/", CancerType,"_dual_add_duration_log.rds"))
  duration_log_df = duration_log_df[which(!is.na(duration_log_df$duration)),]
  
  cancer_bf = cancer_bf %>% mutate( minmax = (relative_importance - min(relative_importance)) / (max(relative_importance) - min(relative_importance)))
  
  cancer_bf_ml = duration_log_df[,cancer_bf$variable]
 
  # cancer_bf_cut_50
  # cancer_bf_cut_100
  # cancer_bf_cut
  # cancer_bf_ml
  
  cancer_bf_ml = cancer_bf_ml[intersect(rownames(cancer_bf_ml)  ,rownames(cancer_short_long_features)),]
  cancer_bf_ml = cancer_bf_ml[rownames(cancer_short_long_features),]
  
  if (all.equal(rownames(cancer_bf_ml) , rownames(cancer_short_long_features))) {
    cancer_bf_ml$cluster = cancer_short_long_features$cluster
  }
  
  # for 50 100 
  short_group_ml = cancer_bf_ml[which(cancer_bf_ml$cluster == "short"),]
  long_group_ml = cancer_bf_ml[which(cancer_bf_ml$cluster == "long"),]
  
  total_group_ml = rbind(long_group_ml,short_group_ml)
  pvals_ml <- as.data.frame(matrix(nrow = c(length(which(!colnames(total_group_ml) %in% "cluster"))))) # -1 means except cluster column
  rownames(pvals_ml) = colnames(total_group_ml)[1:length(which(!colnames(total_group_ml) %in% "cluster"))] # -1 means except cluster column
  colnames(pvals_ml) = "pval"
  
  # t_test for divide genes by long or short 
  for(i in 1:length(which(!colnames(total_group_ml) %in% "cluster"))) { # -1 means except cluster column
    path_short_ml <- short_group_ml[,i]
    path_long_ml <- long_group_ml[,i]
    
    ttest_result_ml <- tryCatch({
      t.test(path_short_ml, path_long_ml)
      
    }, error = function(e) {
      NA
    })
    if (sum(is.na(ttest_result_ml)) != 0) {
      pvals_ml[i,"pval"] <- NA
    }else {
      pvals_ml[i,"pval"] <- ttest_result_ml$p.value
    }
    
  }
  
  cancer_bf$pval = pvals_ml$pval
  
  deg_short_long = rownames(pvals_ml)[which(pvals_ml < 0.05)]
  common_features = rownames(pvals_ml)[!rownames(pvals_ml) %in% deg_short_long]
  
  # dep_group = total_group_ml[,c(deg_short_long,"cluster")]
  # common_group = total_group_ml[,c(common_features,"cluster")]
  
  dep_group = total_group_ml[,deg_short_long]
  common_group = total_group_ml[,common_features]
  
  # cancer_bf_cut_50$pval = pvals_ml$pval[1:50]
  # cancer_bf_cut_100 = pvals_ml$pval[1:100]
  # cancer_bf_cut$pval = pvals_ml$pval[1:nrow(cancer_bf_cut)]
  
  short_cluster_path = c()
  long_cluster_path = c()

  for (deg_path in colnames(dep_group)) {
    if (mean(cancer_bf_ml[which(cancer_bf_ml$cluster == "short"), deg_path]) > mean(cancer_bf_ml[which(cancer_bf_ml$cluster == "long"), deg_path])) {
      short_cluster_path <- c(short_cluster_path, deg_path)
    } else {
      long_cluster_path <- c(long_cluster_path, deg_path)
    }
  }

  rownames(cancer_bf) = NULL
  # cancer_short_long$relative_importance =NULL
  cancer_bf$scaled_importance = NULL
  cancer_bf$percentage = NULL
  cancer_bf$classification = NA
  
  if (length(short_cluster_path) != 0 & length(long_cluster_path) != 0) {
    cancer_bf[which(cancer_bf$variable %in% short_cluster_path),]$classification = "short"
    cancer_bf[which(cancer_bf$variable %in% long_cluster_path),]$classification = "long"
  } else if (length(short_cluster_path) == 0  & length(long_cluster_path) != 0) {
    cancer_bf[which(cancer_bf$variable %in% long_cluster_path),]$classification = "long"
  } else if (length(short_cluster_path) != 0  & length(long_cluster_path) == 0) {
    cancer_bf[which(cancer_bf$variable %in% short_cluster_path),]$classification = "short"
  } else {
    print("I don't know")
  }
 
  if (sum(is.na(cancer_bf$classification)) != 0) {
    cancer_bf[is.na(cancer_bf$classification),]$classification = "common"
  } 
  
  
  if (nrow(cancer_bf) < 50) {
    
    cancer_bf_cut_50 = cancer_bf
    cancer_bf_cut_100 = cancer_bf
    print(paste0(CancerType , "does not have more than 50 features."))
    
  } else if (nrow(cancer_bf) > 50 & nrow(cancer_bf) <100 ) {
    cancer_bf_cut_50 = cancer_bf[1:50,]
    cancer_bf_cut_100 = cancer_bf
    
    print(paste0(CancerType , "does not have more than 100 features."))
  } else {
    cancer_bf_cut_50 = cancer_bf[1:50,]
    cancer_bf_cut_100 = cancer_bf[1:100,]
    
  }
  
  cancer_bf_cut = cancer_bf[1:surv_total_results[which(surv_total_results$CancerType == CancerType),]$num_of_features,]
  
  write.csv(cancer_bf , paste0("/mnt/gluster_server/data/network/",CancerType,"_raw_from_ml_short_long_common.csv"))
  write.csv(cancer_bf_cut_50 , paste0("/mnt/gluster_server/data/network/",CancerType,"_cut50_short_long_common.csv"))
  write.csv(cancer_bf_cut_100 , paste0("/mnt/gluster_server/data/network/",CancerType,"_cut100_short_long_common.csv"))
  write.csv(cancer_bf_cut , paste0("/mnt/gluster_server/data/network/",CancerType,"_critical_features_short_long_common.csv"))
  
}


