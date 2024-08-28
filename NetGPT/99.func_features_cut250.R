library(survival) 
library(dplyr)
library(survminer)

filepath = "/home/seokwon/nas/"
ref_path = paste0(filepath, "99.reference/")

Cancerlist = dir(paste0(filepath, "/00.data/filtered_TCGA/"))
surv_total_results = read_xlsx("~/nas/04.Results/Total_results_survpval2.xlsx")

# num_CancerType = "24.TCGA-OV"
# for all
for (num_CancerType in Cancerlist) {
  
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  
  # call input
  
  cancer_bf = read.csv(paste0(filepath,"04.Results/bestfeatures/",CancerType, "_best_features.csv"))
  cancer_bf_250 = cancer_bf[1:250,]
  cancer_bf_250 = na.omit(cancer_bf_250)
  write.csv(cancer_bf_250 , paste0("/mnt/gluster_server/data/NetGPT/top250/",CancerType,"_top250_from_ml.csv"))
  
}


