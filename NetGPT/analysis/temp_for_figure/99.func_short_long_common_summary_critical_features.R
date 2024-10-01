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
fig_path = "~/nas/04.Results/short_long/ttest_common/"
# if(!dir.exists(fig_path)){
#   dir.create(fig_path)
#   print(paste0("Created folder: ", fig_path))
# } else {
#   print(paste0("Folder already exists: ", fig_path))
# }
setwd(fig_path)
total_sl_count = data.frame( matrix(nrow = 3, ncol = length(Cancerlist)))
rownames(total_sl_count) = c("long","short","common")
colnames(total_sl_count) = gsub('TCGA-','',gsub('[.]','', gsub('\\d','',Cancerlist)))

for (num_CancerType in Cancerlist) {
  
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  name_CancerType = gsub('TCGA-','',CancerType)
  
  tmp_cancer_cri = read.csv(paste0("./",CancerType, "_critical_features_short_long_common.csv"))
  tmp_count = tmp_cancer_cri %>% count(classification)
  
  if (any(tmp_count$classification == 'long')) {
    
    total_sl_count['long',name_CancerType] = tmp_count[which(tmp_count$classification == 'long'),]$n
    } else {
      total_sl_count['long',name_CancerType] = 0
    }
  
  if (any(tmp_count$classification == 'short')) {
    total_sl_count['short',name_CancerType] = tmp_count[which(tmp_count$classification == 'short'),]$n
    } else {
      total_sl_count['short',name_CancerType] = 0
      }

  if (any(tmp_count$classification == 'common')) {
    total_sl_count['common',name_CancerType] = tmp_count[which(tmp_count$classification == 'common'),]$n
    } else {
      total_sl_count['common',name_CancerType] = 0
      }
  
}
  
total_sl_count$classification = rownames(total_sl_count)
write_csv(total_sl_count,"summary_of_total_critical_features.csv" )
