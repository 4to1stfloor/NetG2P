library(dplyr)
library(RColorBrewer)
library(readxl)
library(openxlsx)
library(igraph)
library(ggforce)
library(ggplot2)
library(graphlayouts)
library(ggraph)
library(igraph)
library(stringr)
filepath = "/home/seokwon/nas/"
ref_path = paste0(filepath, "99.reference/")

Cancerlist = dir(paste0(filepath, "/00.data/filtered_TCGA/"))
surv_total_results = read_xlsx("~/nas/04.Results/Total_results_survpval2.xlsx")
# setwd("~/nas/04.Results/short_long/ttest")

total_features = list()

# for all
for (num_CancerType in Cancerlist) {
  
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  Cancernum = gsub('[.-]','',gsub('[a-zA-Z]','', num_CancerType))
  # call input
  
  cancer_bf = read.csv(paste0(filepath,"04.Results/bestfeatures/",CancerType, "_best_features.csv"))
  
  cancer_bf_cut = cancer_bf[1:surv_total_results[which(surv_total_results$CancerType == CancerType),]$num_of_features,]
  total_features[[CancerType]] = cancer_bf_cut$variable
  
}


total_common = list()

second_Cancerlist = Cancerlist

for (num_CancerType in Cancerlist) {
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  second_Cancerlist = second_Cancerlist[-1]
  
  for (second_cancertype in second_Cancerlist) {
    second_CancerType = gsub('[.]','',gsub('\\d','', second_cancertype))
    
    total_common[[paste0(CancerType,'-',second_CancerType)]] = intersect(total_features[[CancerType]] , total_features[[second_CancerType]])
  }
  
}

total_common_ratio = data.frame(interact_cancer = paste0(names(total_common), 1:length(total_common)),
                                critical_featrues = sapply(total_common, function(x) paste(x, collapse = ",")),
                                count = sapply(total_common, length),
                                stringsAsFactors = FALSE)

total_common_ratio$interact_cancer = rownames(total_common_ratio)

# Split the elements in the interact_cancer column
split_elements <- str_split(total_common_ratio$interact_cancer, "-")

# Extract the from and to values
from <- sapply(split_elements, function(x) paste(x[1:2], collapse = "-"))
to <-  sapply(split_elements, function(x) paste(x[3:4], collapse = "-"))

# Create new columns in the dataframe
total_common_ratio$from <- from
total_common_ratio$to <- to

total_common_ratio$from_ratio = total_common_ratio$count  / surv_total_results$num_of_features[match(total_common_ratio$from , surv_total_results$CancerType)]
total_common_ratio$to_ratio = total_common_ratio$count  / surv_total_results$num_of_features[match(total_common_ratio$to , surv_total_results$CancerType)]
total_common_ratio$mean_ratio = rowMeans(total_common_ratio[, c("from_ratio", "to_ratio")])


write.xlsx(total_common_ratio, paste0(filepath, "04.Results/bestfeatures/shared_critical_features_each_cancer_fixed.xlsx"))
