
library(tidyverse)
library(readxl)
library(tidygraph)
library(data.table)

filepath = "/home/seokwon/nas/"
ref_path = paste0(filepath, "99.reference/")
Cancerlist = dir(paste0(filepath, "/00.data/filtered_TCGA"))
num_CancerType = "10.TCGA-BLCA"
results_df = data.frame(matrix(ncol = 5))
colnames(results_df) = c("original", "three", "five", "median", "total_median")
for (num_CancerType in Cancerlist) {
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  df_dual = read.csv(paste0(main.path_tc,"/",CancerType,"_dual_add_duration_status_log.csv"))
  df_dual = df_dual[complete.cases(df_dual[, "duration"]), ]
  cli = fread(paste0(ref_path,"all_clin_indexed.csv"))
  cli_type = cli[which(cli$project == CancerType),]
  cli_type_filt = cli_type[,c("submitter_id","days_to_last_follow_up","days_to_death","vital_status")]
  cli_type_filt = cli_type_filt[which(cli_type_filt$vital_status !="Not Reported",),]
  cli_type_filt = cli_type_filt[which(!is.na(cli_type_filt$vital_status)),]
  cli_type_filt$duration = ifelse(cli_type_filt$vital_status == "Alive", cli_type_filt$days_to_last_follow_up, cli_type_filt$days_to_death)
  cli_type_filt = cli_type_filt[which(!is.na(cli_type_filt$duration)),]
  cli_type_filt = cli_type_filt[which(cli_type_filt$duration > 1),]
  total_sample_median = median(cli_type_filt$duration)
  
  df_dual$three_years = ifelse(df_dual$duration < 1095 , df_dual$status , ifelse(df_dual$duration > 1095 & df_dual$status == 1, 0, 0))
  df_dual$five_years = ifelse(df_dual$duration < 1825 , df_dual$status , ifelse(df_dual$duration > 1825 & df_dual$status == 1, 0, 0))
  df_dual$median_years = ifelse(df_dual$duration < median(df_dual$duration), df_dual$status ,ifelse(df_dual$duration > median(df_dual$duration) & df_dual$status == 1, 0, 0))
  df_dual$total_median_years = ifelse(df_dual$duration < total_sample_median, df_dual$status ,ifelse(df_dual$duration > total_sample_median & df_dual$status == 1, 0, 0))
  
 
  results_df[CancerType,"original"] = sum(df_dual$status == 1) / length(df_dual$status)
  results_df[CancerType,"three"] = sum(df_dual$three_years == 1) / length(df_dual$three_years)
  results_df[CancerType,"five"] = sum(df_dual$five_years == 1) / length(df_dual$five_years)
  results_df[CancerType,"median"] = sum(df_dual$median_years == 1) / length(df_dual$median_years)
  results_df[CancerType,"total_median"] = sum(df_dual$total_median_years == 1) / length(df_dual$total_median_years)
  saveRDS(df_dual, paste0(main.path_tc,"/",CancerType, "_dual_add_duration_all_status_log.rds"))
  write.csv(df_dual, paste0(main.path_tc,"/",CancerType, "_dual_add_duration_all_status_log.csv"))
}

results_df= results_df[-1,]
write.csv(results_df, paste0(filepath, "/00.data/filtered_TCGA/vias_results.csv"))
  
  
  
  
  
  