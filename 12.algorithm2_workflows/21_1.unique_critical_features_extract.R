library(readxl)
library(tidyverse)
filepath = "/home/seokwon/nas/"
ref_path = paste0(filepath, "99.reference/")

Cancerlist = dir(paste0(filepath, "/00.data/filtered_TCGA/"))
surv_total_results = read_xlsx("~/nas/04.Results/Total_results_survpval2.xlsx")
# setwd("~/nas/04.Results/short_long/ttest")
Cancerlist = Cancerlist[c(-7,-11,-12)]
total_features = list()

# for all
for (num_CancerType in Cancerlist) {
  
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  Cancernum = gsub('[.-]','',gsub('[a-zA-Z]','', num_CancerType))
  # call input
  
  cancer_bf = read.csv(paste0(filepath,"04.Results/bestfeatures/",CancerType, "_critical_features.csv"))
  
  # cancer_bf_cut = cancer_bf[1:surv_total_results[which(surv_total_results$CancerType == CancerType),]$num_of_features,]
  total_features[[CancerType]] = cancer_bf$variable
  
}
total_spe = data.frame()

for (num_CancerType in Cancerlist) {
  
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  
  spe_cancer = total_features[[CancerType]]
  # call input
  for (diff_cancer in names(total_features)[names(total_features) != CancerType]) {
    
    spe_cancer = setdiff(spe_cancer, total_features[[diff_cancer]])
  }
  
  cancer_slc = read.csv(paste0("~/nas/04.Results/short_long/ttest_common/",CancerType,"_critical_features_short_long_common.csv"))
  tmp_spe = cancer_slc[which(cancer_slc$variable %in% spe_cancer),]
  
  tmp_sep_filt = tmp_spe %>% dplyr::select(variable, classification) %>% 
    mutate(long_rat = case_when( classification == "long" ~ 1,.default = 0)) %>% 
    mutate(short_rat = case_when( classification == "short" ~ 1,.default = 0)) %>%
    mutate(common_rat = case_when( classification == "common" ~ 1,.default = 0)) %>% 
    dplyr::select(variable, long_rat, short_rat, common_rat) %>%
    mutate(which_cancer = CancerType) 
  
  total_spe = rbind(total_spe , tmp_sep_filt)
  
  remove(tmp_sep_filt)
  
}

for (num_CancerType in Cancerlist) {
  
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  spe_cancer = total_features[[CancerType]]
  # call input
  for (diff_cancer in names(total_features)[names(total_features) != CancerType]) {
    spe_cancer = setdiff(spe_cancer, total_features[[diff_cancer]])
  }
  
  cancer_slc = read.csv(paste0("~/nas/04.Results/short_long/ttest_common/",CancerType,"_critical_features_short_long_common.csv"))
  tmp_spe = cancer_slc[which(cancer_slc$variable %in% spe_cancer),]
  
  # tmp_spe_filt = tmp_spe %>% dplyr::select(classification)
  tmp_spe_filt = tmp_spe$classification
  
  total_spe[,CancerType] = NA
  total_spe[which(total_spe$which_cancer == CancerType),][,CancerType] = tmp_spe_filt
  
  
}

write.csv(total_spe, "~/nas/04.Results/critical_features/unique_short_long_critical_features_for_SW_ver2.csv" )
