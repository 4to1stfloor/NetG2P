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
library(tidyr)
library(dplyr)


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

# Get the maximum length
max_length <- max(lengths(total_features))

# Pad or truncate vectors to match the maximum length
total_features_df <- total_features %>%
  purrr::map(function(vec) {
    if (length(vec) < max_length) {
      c(vec, rep(NA, max_length - length(vec)))
    } else if (length(vec) > max_length) {
      head(vec, max_length)
    } else {
      vec
    }
  }) %>%
  bind_cols()

# Set the column names
colnames(total_features_df) <- names(total_features)
saveRDS(total_features_df, "~/nas/04.Results/bestfeatures/total_critical_feartures_merge.rds")
write.xlsx(total_features_df, "~/nas/04.Results/bestfeatures/total_critical_feartures_merge.xlsx")

###
library(tidyverse)

# Reshape the dataframe to long format
df_long <- total_features_df %>%
  pivot_longer(everything(), names_to = "which_cancer", values_to = "features") %>%
  filter(!is.na(features))  # Remove rows with NA values

# Count the occurrences of each feature and aggregate cancer types
feature_counts <- df_long %>%
  group_by(features) %>%
  summarise(which_cancer = paste(which_cancer, collapse = ", "),
            count = n()) %>%
  arrange(desc(count))

# Print the top feature and their counts along with the corresponding cancers
top_feature <- feature_counts$features[1]  # Get the top feature (first row)

result <- feature_counts %>%
  filter(features == top_feature) %>%
  select(features, count, which_cancer)

saveRDS(feature_counts, "~/nas/04.Results/bestfeatures/most_common_critical_features.rds")
write.xlsx(feature_counts, "~/nas/04.Results/bestfeatures/most_common_critical_features.xlsx")

#####
feature_counts_with_importance = df_long %>%
  arrange(features) %>%  group_by(features) %>% mutate(count = n())

feature_counts_with_importance$importance = NA

for (num_CancerType in Cancerlist) {
  
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  Cancernum = gsub('[.-]','',gsub('[a-zA-Z]','', num_CancerType))
  # call input
  
  cancer_bf = read.csv(paste0(filepath,"04.Results/bestfeatures/",CancerType, "_best_features.csv"))
  
  cancer_bf_cut = cancer_bf[1:surv_total_results[which(surv_total_results$CancerType == CancerType),]$num_of_features,]
  cancer_bf_cut = cancer_bf_cut %>% mutate( minmax = (relative_importance - min(relative_importance)) / (max(relative_importance) - min(relative_importance)))
   
  for (critical_feature in feature_counts_with_importance[which(feature_counts_with_importance$which_cancer == CancerType),]$features) {
    feature_counts_with_importance[which(feature_counts_with_importance$which_cancer == CancerType & feature_counts_with_importance$features == critical_feature),]$importance = 
      cancer_bf_cut[cancer_bf_cut$variable == critical_feature,]$minmax
  }
  
}


features_importance_wo = feature_counts_with_importance[which(feature_counts_with_importance$which_cancer %in% c("TCGA-CESC","TCGA-OV", "TCGA-BRCA", "TCGA-UCEC")),]
write.xlsx(features_importance_wo , "~/nas/04.Results/bestfeatures/critical_features_imp_in_woman_cancer.xlsx")


