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

cut_num = c(50,100)

# for all
for (num_fe in cut_num) {
  for (num_CancerType in Cancerlist) {
    
    main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
    CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
    
    # call input
    
    cancer_bf = read.csv(paste0(filepath,"04.Results/bestfeatures/",CancerType, "_best_features.csv"))
    
    if (nrow(cancer_bf) > num_fe) {
      cancer_bf = cancer_bf[1:num_fe,]
    } else {
      cancer_bf = cancer_bf[1:nrow(cancer_bf),]
    }

    write_csv(cancer_bf , paste0("/mnt/gluster_server/data/network/",CancerType,"_cut",num_fe,"_short_long.csv"))
    
    
  }  
  
}
