library(survival) 
library(dplyr)
library(survminer)

filepath = "/home/seokwon/nas/"
ref_path = paste0(filepath, "99.reference/")

Cancerlist = dir(paste0(filepath, "/00.data/filtered_TCGA/"))
surv_total_results = read.csv("~/nas/04.Results/Total_results_survpval.csv")

setwd("~/nas/04.Results/short_long/")
num_CancerType = "24.TCGA-OV"
# for all
for (num_CancerType in Cancerlist) {
  
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  
  # call input
  
  cancer_bf = read.csv(paste0(filepath,"04.Results/bestfeatures/",CancerType, "_best_features.csv"))
  duration_log_df = readRDS(paste0(main.path_tc, "/", CancerType,"_dual_add_duration_log.rds"))
  
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
  
  data_tc_link = readRDS(file = paste0(main.path_tc,"/",CancerType,"_pathwaylink_all_log.rds"))
  data_tc_each = readRDS(file = paste0(main.path_tc,"/",CancerType,"_pathwayeach_all_log.rds"))
  
  data_tc_link_filt = data_tc_link[,which(colSums(data_tc_link[-ncol(data_tc_link)]) != 0)]
  data_tc_each_filt = data_tc_each[,which(colSums(data_tc_each[-ncol(data_tc_each)]) != 0)]
  
  data_tc_link_filt$vitalstatus = data_tc_link$vitalstatus
  data_tc_each_filt$vitalstatus = data_tc_each$vitalstatus
  
  data_tc_link_filt = data_tc_link_filt[intersect(rownames(data_tc_link_filt), rownames(data_tc_each_filt)),]
  data_tc_each_filt = data_tc_each_filt[intersect(rownames(data_tc_link_filt), rownames(data_tc_each_filt)),]
  
  if (all.equal(rownames(data_tc_link_filt), rownames(data_tc_each_filt))) {
    data_tc_link_filt$vitalstatus = NULL
    data_tc = cbind(data_tc_link_filt,data_tc_each_filt)
  }
  cancer_filtered = cancer_bf[1:surv_total_results[which(surv_total_results$CancerType == CancerType),]$num_of_features,]
  
  data_bf =  data_tc[,cancer_filtered$variable]
  # annotation cluster 1 or 2 (long or short) by best pval score
  # 18 means best p-value when devide two cluster
  out = pheatmap::pheatmap((data_bf > -log(0.05))*1 ,
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
  
  if (nrow(data_bf) == nrow(duration_log_df)) {
    if (all.equal(rownames(data_bf), rownames(duration_log_df))) {
      data_bf$vital_status = duration_log_df$vitalstatus
      data_bf$duration = duration_log_df$duration
    }
    
  } else {
    common_pat2 = intersect(rownames(data_bf), rownames(duration_log_df))
    data_bf = data_bf[common_pat2,]
    duration_log_df = duration_log_df[common_pat2,]
    if (all.equal(rownames(data_bf), rownames(duration_log_df))) {
      data_bf$vital_status = duration_log_df$vitalstatus
      data_bf$duration = duration_log_df$duration
    }
    
  }
  
  data_bf$status = NA
  data_bf$status = ifelse(data_bf$vital_status == "Alive", 0 , 1)
  data_bf$vital_status = NULL
  data_bf = na.omit(data_bf)
  
  data_bf = data_bf[which(data_bf$duration >= 0),]
  
  data_bf2 = data_bf[,which(!colnames(data_bf) %in% c("cluster","duration","status"))]
  data_bf2[data_bf2 > 5] = 5
  data_bf2[data_bf2 < -5] = -5
  
  data_bf_max_df = cbind( data_bf2, data_bf[,which(colnames(data_bf) %in% c("cluster","duration","status"))])
  
  fit = survfit(Surv(duration, status) ~ cluster, data = data_bf_max_df)
  # ggsurvplot(fit, data = data_bf_max_df, risk.table = TRUE,
  #            palette = "jco", pval = TRUE, surv.median.line = "hv", xlab = "days")
  
  name_of_50_100 = c("cancer_bf_cut_50","cancer_bf_cut_100")
  n = 0

  for (cut in name_of_50_100) {
    n = n + 1
    data_cut = data_tc[,get(cut)$variable]
    
    if (nrow(data_cut) == nrow(data_bf_max_df)) {
      
      if (all.equal( rownames(data_cut), rownames(data_bf_max_df))) {
        data_cut$cluster = data_bf_max_df$cluster
        
      }
      
    } else {
      common_pat = intersect(rownames(data_cut) ,rownames(data_bf_max_df) )
      data_cut = data_cut[common_pat,]
      data_bf_max_df = data_bf_max_df[common_pat,]
      if (all.equal( rownames(data_cut), rownames(data_bf_max_df))) {
        data_cut$cluster = data_bf_max_df$cluster
      
      }
      
    }
    # hypothesis : cluster 1 = better prognosis
    if (mean(fit$surv[1:fit[['strata']][['cluster=1']]]) > mean(fit$surv[fit[['strata']][['cluster=1']] + 1: fit[['strata']][['cluster=2']]])) {
      short_group = data_cut[which(data_cut$cluster == 2),]
      long_group = data_cut[which(data_cut$cluster == 1),]
      
      short_group_for_fig = data_bf_max_df[which(data_bf_max_df$cluster == 2),]
      long_group_for_fig = data_bf_max_df[which(data_bf_max_df$cluster == 1),]
      
    } else if (mean(fit$surv[1:fit[['strata']][['cluster=1']]]) < mean(fit$surv[fit[['strata']][['cluster=1']] + 1: fit[['strata']][['cluster=2']]])) {
      short_group = data_cut[which(data_cut$cluster == 1),]
      long_group = data_cut[which(data_cut$cluster == 2),]
      
      short_group_for_fig = data_bf_max_df[which(data_bf_max_df$cluster == 2),]
      long_group_for_fig = data_bf_max_df[which(data_bf_max_df$cluster == 1),]
    } else {
      print("I don't know")
    }
    
    short_group = na.omit(short_group)
    long_group = na.omit(long_group)
    
    short_group_for_fig = na.omit(short_group_for_fig)
    long_group_for_fig = na.omit(long_group_for_fig)
    
    # for ttest method (short vs long)
    short_group$cluster = "short" 
    long_group$cluster = "long" 
    
    total_group = rbind(long_group,short_group)
    
    pvals <- as.data.frame(matrix(nrow = c(ncol(total_group)-1))) # -1 means except cluster column
    rownames(pvals) = colnames(total_group)[1:(ncol(total_group)-1)] # -1 means except cluster column
    colnames(pvals) = "pval"
    
    # t_test for divide genes by good or bad 
    for(i in 1:c(ncol(total_group)-1)) { # -1 means except cluster column
      path_short <- short_group[,i]
      path_long <- long_group[,i]
      
      ttest_result <- tryCatch({
        t.test(path_short, path_long)
        
      }, error = function(e) {
        NA
      })
      if (sum(is.na(ttest_result)) != 0) {
        pvals[i,"pval"] <- NA
      }else {
        pvals[i,"pval"] <- ttest_result$p.value
      }
      
    }
    
    deg_short_long = rownames(pvals)[which(pvals < 0.05)]
    deg_group = total_group[,c(deg_short_long,"cluster")]
    
    short_cluster_path = c()
    long_cluster_path = c()
    
    for (deg_path in deg_short_long) {
      if (mean(deg_group[which(deg_group$cluster == "short"), deg_path]) > mean(deg_group[which(deg_group$cluster == "long"), deg_path])) {
        short_cluster_path <- c(short_cluster_path, deg_path)
      } else {
        long_cluster_path <- c(long_cluster_path, deg_path)
      }
    }
    
    hreason = c(0.1,0.2,0.3,0.4,0.5)
    
    if (length(short_cluster_path) != 0 && length(long_cluster_path) != 0) {
      print("There are well divided as short or long pathway")
      
    } else {
      
      for (not_spe in hreason){
        deg_short_long = rownames(pvals)[which(pvals < not_spe)]
        deg_group = total_group[,c(deg_short_long,"cluster")]
        
        short_cluster_path = c()
        long_cluster_path = c()
        
        for (deg_path in deg_short_long) {
          if (mean(deg_group[which(deg_group$cluster == "short"), deg_path]) > mean(deg_group[which(deg_group$cluster == "long"), deg_path])) {
            short_cluster_path <- c(short_cluster_path, deg_path)
          } else {
            long_cluster_path <- c(long_cluster_path, deg_path)
          }
        }
        
        if (length(short_cluster_path) != 0 && length(long_cluster_path) != 0) {
          print(paste0("The least pvalue that are divided by short and long is ",not_spe))
          break
        }
      }
    }
    merge_short_long = c(short_cluster_path,long_cluster_path)

    cancer_short_long = cancer_bf[which(cancer_bf$variable %in% merge_short_long),]
    rownames(cancer_short_long) = NULL
    cancer_short_long$relative_importance =NULL
    cancer_short_long$scaled_importance = NULL
    cancer_short_long$percentage = NULL
    cancer_short_long$classification = NA
    
    cancer_short_long[which(cancer_short_long$variable %in% short_cluster_path),]$classification = "short"
    cancer_short_long[which(cancer_short_long$variable %in% long_cluster_path),]$classification = "long"
    
    if (n == 1) {
      write.csv(cancer_short_long , paste0("~/nas/04.Results/short_long/",CancerType,"_cut50_short_long.csv"))
    } else if (n == 2) {
      write.csv(cancer_short_long , paste0("~/nas/04.Results/short_long/",CancerType,"_cut100_short_long.csv"))
    }
    
    
  }
  
  png(filename = paste0(CancerType,"_short_log_pval.png"),
      width = 35, height = 35,  units = "cm", pointsize = 12,
      bg = "white", res = 1200, family = "")
  
  short_out = pheatmap::pheatmap(short_group_for_fig %>% dplyr::select(-cluster,-duration,-status), 
                                 cluster_cols = T,
                                 cluster_rows = T, 
                                 labels_cols = "",
                                 show_rownames = T,
                                 silent = T)
  
  print(short_out)
  dev.off()
  
  png(filename = paste0(CancerType,"_long_log_pval.png"),
      width = 35, height = 35,  units = "cm", pointsize = 12,
      bg = "white", res = 1200, family = "")
  
  long_out = pheatmap::pheatmap(long_group_for_fig %>% dplyr::select(-cluster,-duration,-status), 
                                cluster_cols = T,
                                cluster_rows = T, 
                                labels_cols = "",
                                show_rownames = T,
                                silent = T)
  
  print(long_out)
  dev.off()

}  


