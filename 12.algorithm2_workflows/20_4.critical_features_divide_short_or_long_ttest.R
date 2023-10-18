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

setwd("~/nas/04.Results/short_long/ttest")
# num_CancerType = "19.TCGA-LIHC"
# for all
for (num_CancerType in Cancerlist) {
  
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  
  # call input

  cancer_short_long_features = readRDS(paste0("~/nas/04.Results/short_long/",CancerType,"_critical_features_short_long.rds"))
  short_group = cancer_short_long_features[which(cancer_short_long_features$cluster == "short"),]
  long_group = cancer_short_long_features[which(cancer_short_long_features$cluster == "long"),]
  
  wo_num = ncol(cancer_short_long_features) - length(grep("*P", colnames(cancer_short_long_features)))

  pvals <- as.data.frame(matrix(nrow = c(ncol(cancer_short_long_features)-wo_num))) # -1 means except cluster column
  rownames(pvals) = colnames(cancer_short_long_features)[1:(ncol(cancer_short_long_features)-wo_num)] # -1 means except cluster column
  colnames(pvals) = "pval"
  
  # t_test for divide genes by good or bad
  for(features in colnames(cancer_short_long_features)[which(!colnames(cancer_short_long_features) %in% c("vitalstatus","duration", "status", "cluster"))]) { # -1 means except cluster column
    path_short <- short_group[,features]
    path_long <- long_group[,features]

    ttest_result <- tryCatch({
      t.test(path_short, path_long)

    }, error = function(e) {
      NA
    })
    if (is.na(ttest_result$p.value)) {
      pvals[features,"pval"] <- NA
    } else {
      pvals[features,"pval"] <- ttest_result$p.value
    }

  }

  deg_short_long = rownames(pvals)[which(pvals$pval < 0.05)]
  deg_group = cancer_short_long_features[,c(deg_short_long,"cluster")]
  
  # fig

  total_group_meta = as.data.frame(deg_group[,c("cluster")])
  rownames(total_group_meta) = rownames(deg_group)
  colnames(total_group_meta) = "cluster"
  total_group_for_tmp = deg_group[,which(!colnames(deg_group) %in% c("cluster"))]

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
  # names(short_long_colors) <- unique(annotation_df$short_long)

  short_long_colors = list(cluster = short_long_colors )

  png(filename = paste0(CancerType,"_short_long_complex_ttest_quant.png"),
      width = 35, height = 35,  units = "cm" ,pointsize = 12,
      bg = "white", res = 1200, family = "")

  total_out = ComplexHeatmap::pheatmap(as.matrix(t(total_group_for_fig %>%
                                                     dplyr::select_if(is.numeric) %>%
                                                     as.matrix())),
                                       column_split = factor(annotation_df$cluster, levels = c("short","long")),
                                       annotation_col = annotation_df,
                                       annotation_colors = short_long_colors,
                                       cluster_cols = T,
                                       legend = F,
                                       annotation_legend = F,
                                       show_colnames = F,
                                       cluster_column_slices = FALSE)

  print(total_out)

  dev.off()

  png(filename = paste0(CancerType,"_short_long_cluster_ttest_quant.png"),
      width = 25, height = 25,  units = "cm" ,pointsize = 12,
      bg = "white", res = 1200, family = "")

  total_out2 = pheatmap::pheatmap(as.matrix(t(total_group_for_fig %>%
                                                dplyr::select_if(is.numeric) %>%
                                                as.matrix())),
                                  annotation_col = annotation_df,
                                  annotation_colors = short_long_colors,
                                  cluster_cols = TRUE,
                                  legend = T)
  print(total_out2)

  dev.off()
  
  
}  
