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
surv_total_results = read.xlsx("~/nas/04.Results/Total_results_survpval.xlsx")

setwd("~/nas/04.Results/short_long/")

total_best_features = data.frame()

nor_minmax = function(x){
  result = (x - min(x)) / (max(x) - min(x))
  return(result)
}

# for all
for (num_CancerType in Cancerlist) {
  
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))

  # call input
  
  cancer_bf = read.csv(paste0(filepath,"04.Results/bestfeatures/",CancerType, "_best_features.csv"))
  duration_log_df = readRDS(paste0(main.path_tc, "/", CancerType,"_dual_add_duration_log.rds"))
  
  # cut the number of best features 

  best_features_df = duration_log_df[,cancer_bf$variable]
  # cancer_spe_sum = data.frame(features = colnames(best_features_df) , (colSums(best_features_df) / nrow(best_features_df)), row.names = NULL)

  cancer_spe_sum_minmax = as.data.frame(apply(as.data.frame(colSums(best_features_df) / nrow(best_features_df)), MARGIN = 2, FUN = "nor_minmax"))
  
  cancer_spe_sum_minmax$features = rownames(cancer_spe_sum_minmax)
  colnames(cancer_spe_sum_minmax) = c(CancerType , "features")

  cancer_spe_sum_minmax = cancer_spe_sum_minmax %>% select(all_of(c("features",CancerType)))
  
  if (CancerType == "TCGA-CESC") {
    total_best_features = cancer_spe_sum_minmax
  } else {
    total_best_features = full_join(total_best_features,cancer_spe_sum_minmax,by="features")
  }
 
}
  
total_best_features[is.na(total_best_features)] = 0
total_best_features = total_best_features[order(total_best_features$features),]
rownames(total_best_features) = total_best_features$features
total_best_features$features =NULL

# Convert the matrix to a numeric matrix
total_numeric <- matrix(as.numeric(unlist(total_best_features)), nrow = nrow(total_best_features))

# Convert the numeric matrix to a vector
vec <- as.vector(total_numeric)

# Create the histogram
hist(vec)

boxplot(vec)

total_best_features[total_best_features > quantile(vec, probs = 0.75) + 1.5*IQR(vec)] <- round(quantile(vec, probs = 0.75) + 1.5*IQR(vec),1)

quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}
library(viridis)

mat_breaks <- quantile_breaks(total_numeric, n = 11)

# install.packages("dendsort")
library(dendsort)

sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
mat_cluster_cols <- hclust(dist(t(total_best_features)), method = "complete")
# plot(mat_cluster_cols, main = "Unsorted Dendrogram", xlab = "", sub = "")
mat_cluster_cols <- sort_hclust(mat_cluster_cols)
# plot(mat_cluster_cols, main = "Sorted Dendrogram", xlab = "", sub = "")
mat_cluster_rows <- sort_hclust(hclust(dist(total_best_features), method = "complete"))
# plot(hclust(dist(total_best_features)), main = "Sorted Dendrogram", xlab = "", sub = "")

print(pheatmap::pheatmap(total_best_features ,
                   cluster_cols = mat_cluster_cols,
                   cluster_rows = mat_cluster_rows,
                   labels_cols = "", 
                   show_rownames = F,
                   silent = T,
                   border_color = 'white',
                   color = inferno(length(mat_breaks) - 1)
))
