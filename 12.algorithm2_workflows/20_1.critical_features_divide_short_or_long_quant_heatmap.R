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

setwd("~/nas/04.Results/short_long/quantile")

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

  fit = survfit(Surv(duration, status) ~ cluster, data = data_bf)
  # ggsurvplot(fit, data = data_bf, risk.table = TRUE,
  #            palette = "jco", pval = TRUE, surv.median.line = "hv", xlab = "days")

  # hypothesis : cluster 1 = better prognosis
  if (mean(fit$surv[1:fit[['strata']][['cluster=1']]]) > mean(fit$surv[fit[['strata']][['cluster=1']] + 1: fit[['strata']][['cluster=2']]])) {
    short_group = data_bf[which(data_bf$cluster == 2),]
    long_group = data_bf[which(data_bf$cluster == 1),]
    
    short_group_for_fig = data_bf[which(data_bf$cluster == 2),]
    long_group_for_fig = data_bf[which(data_bf$cluster == 1),]
    
  } else if (mean(fit$surv[1:fit[['strata']][['cluster=1']]]) < mean(fit$surv[fit[['strata']][['cluster=1']] + 1: fit[['strata']][['cluster=2']]])) {
    short_group = data_bf[which(data_bf$cluster == 1),]
    long_group = data_bf[which(data_bf$cluster == 2),]
    
    short_group_for_fig = data_bf[which(data_bf$cluster == 1),]
    long_group_for_fig = data_bf[which(data_bf$cluster == 2),]
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

  # saveRDS(total_group, paste0("~/nas/04.Results/short_long/",CancerType,"_best_features_short_long.rds"))

  wo_num = ncol(total_group) - length(grep("*P", colnames(total_group)))

  pvals <- as.data.frame(matrix(nrow = c(ncol(total_group)-wo_num))) # -1 means except cluster column
  rownames(pvals) = colnames(total_group)[1:(ncol(total_group)-wo_num)] # -1 means except cluster column
  colnames(pvals) = "pval"

  # t_test for divide genes by good or bad
  for(i in 1:c(ncol(total_group)-wo_num)) { # -1 means except cluster column
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
    
    if (length(short_cluster_path) == 0 || length(long_cluster_path) == 0) {
      print(CancerType)
    }
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

  # write.xlsx(cancer_short_long , paste0("~/nas/04.Results/short_long/",CancerType,"_best_features_short_long.xlsx"))

  # fig
  short_group_for_fig$cluster = "short"
  long_group_for_fig$cluster = "long"

  total_group_for_fig = rbind(short_group_for_fig,long_group_for_fig)

  total_group_meta = total_group_for_fig[,c("vitalstatus","cluster", "duration", "status")]
  total_group_for_tmp = total_group_for_fig[,which(!colnames(total_group_for_fig) %in% c("vitalstatus","cluster", "duration", "status"))]

  # # Convert the matrix to a numeric matrix
  
  tmp_numeric <- matrix(as.numeric(unlist(total_group_for_tmp)), nrow = nrow(total_group_for_tmp))
  vec <- as.vector(tmp_numeric)
  # hist(vec)
  total_group_for_tmp[total_group_for_tmp > quantile(vec, probs = 0.75) + 1.5*IQR(vec)] <- round(quantile(vec, probs = 0.75) + 1.5*IQR(vec),1)
  
  total_group_for_fig = cbind(total_group_for_tmp,total_group_meta)

  # pic
  annotation_df <- data.frame(cluster = total_group_for_fig$cluster)
  rownames(annotation_df) <- rownames(total_group_for_fig)

  short_long_colors <- c("short" = "red", "long" = "#009E73")
  # names(short_long_colors) <- unique(annotation_df$cluster)

  short_long_colors = list(cluster = short_long_colors )

  png(filename = paste0(CancerType,"_short_long_complex_pval_quant.png"),
      width = 25, height = 25,  units = "cm" ,pointsize = 12,
      bg = "white", res = 1200, family = "")

  Colors = brewer.pal(9, "YlOrRd")
  
  ComplexHeatmap::pheatmap(as.matrix(t(total_group_for_fig %>%
                                         dplyr::select_if(is.numeric) %>%
                                         dplyr::select(-duration,-status) %>%
                                         as.matrix())),
                           column_split = factor(annotation_df$cluster, levels = c("short","long")),
                           annotation_col = annotation_df,
                           annotation_colors = short_long_colors,
                           cluster_cols = T,
                           legend = T,
                           annotation_legend = F,
                           show_colnames = F,
                           cluster_column_slices = FALSE,
                           color = Colors) 

  total_out = ComplexHeatmap::pheatmap(as.matrix(t(total_group_for_fig %>%
                                                     dplyr::select_if(is.numeric) %>%
                                                     dplyr::select(-duration,-status) %>%
                                                     as.matrix())),
                                       column_split = factor(annotation_df$cluster, levels = c("short","long")),
                                       annotation_col = annotation_df,
                                       annotation_colors = short_long_colors,
                                       cluster_cols = T,
                                       legend = T,
                                       annotation_legend = F,
                                       show_colnames = F,
                                       cluster_column_slices = FALSE,color = Colors) 

  print(total_out)

  dev.off()

  png(filename = paste0(CancerType,"_short_long_cluster_pval_quant.png"),
      width = 25, height = 25,  units = "cm" ,pointsize = 12,
      bg = "white", res = 1200, family = "")

  total_out2 = pheatmap::pheatmap(as.matrix(t(total_group_for_fig %>%
                                                dplyr::select_if(is.numeric) %>%
                                                dplyr::select(-duration,-status) %>%
                                                as.matrix())),
                                  annotation_col = annotation_df,
                                  annotation_colors = short_long_colors,
                                  cluster_cols = TRUE,
                                  legend = T)
  print(total_out2)

  dev.off()

  
}  
