library(stringr)
library(data.table)
library(readxl)
library(dplyr)
library(readxl)

filepath = "/home/seokwon/nas/"
ref_path = paste0(filepath, "99.reference/")

Cancerlist = dir(paste0(filepath, "/00.data/filtered_TCGA/"))

surv_total_results = read_xlsx("~/nas/04.Results/Total_results_survpval2.xlsx")

fig_path = paste0(filepath,"04.Results/bestfeatures/for_deep")
if(!dir.exists(fig_path)){
  dir.create(fig_path)
  print(paste0("Created folder: ", fig_path))
} else {
  print(paste0("Folder already exists: ", fig_path))
}
setwd(fig_path)

# for all
for (num_CancerType in Cancerlist) {
  
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  
  # call input
  cancer_bf = read.csv(paste0(filepath,"04.Results/bestfeatures/",CancerType, "_best_features.csv"))
  short_long_features = read_xlsx(paste0(filepath, "04.Results/short_long/",CancerType, "_best_features_short_long.xlsx"))
  bf_short_long = readRDS(paste0(filepath, "04.Results/short_long/",CancerType,"_best_features_short_long.rds"))
  duration_log_df = readRDS(paste0(main.path_tc, "/", CancerType,"_dual_add_duration_log.rds"))
  
  bf_features = cancer_bf[1:surv_total_results[which(surv_total_results$CancerType == CancerType),]$num_of_features,]$variable
  bf_features = c(bf_features,"vitalstatus","duration")
  duration_log_df_bf = duration_log_df[, bf_features]
  
  duration_log_df_bf$status = ifelse(duration_log_df_bf$vitalstatus == "Alive", 0,1)
  duration_log_df_bf$vitalstatus = NULL

  write.xlsx(duration_log_df_bf, paste0(fig_path, "/", CancerType ,"_best_features_cut_log.xlsx"),rowNames=TRUE)
}  


