
library(tidyverse)
library(readxl)
library(h2o)

filepath = "/home/seokwon/"
ref_path = paste0(filepath, "99.reference/")
Cancerlist = dir(paste0(filepath, "/00.data/"))
# Cancerlist = Cancerlist[7:9]
Cancerlist = Cancerlist[c(-2,-5,-15)]
folder_name = "h2o_bias_pval_dual_cut_50"
save_folder = "13.analysis"
localH2O = h2o.init(ip = "localhost", port = sample(x=00000:65536,size=1,replace=F), 
                    startH2O = TRUE,min_mem_size = "400G",
                    nthreads = 96,enable_assertions = FALSE)

find_max_model <- function(vec) {
  vec <- vec[!grepl(".csv", vec)]  # remove elements with ".csv" extension
  max_model <- vec[which.max(as.numeric(str_sub(vec, -5,-1) ))]
  return(max_model)
}


for (num_CancerType in Cancerlist) {
  main.path_tc = paste0(filepath, "00.data/", num_CancerType, "/", folder_name)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  dual_cut_50 = list.files(paste0(main.path_tc), pattern = "model")
  best_model_folder = find_max_model(dual_cut_50)
  best_model = h2o.loadModel(paste0(main.path_tc,"/",best_model_folder,"/", gsub('.{6}$', '', best_model_folder)))
  best_model_vari_meta = as.data.frame(h2o.varimp(best_model))

  write_csv(best_model_vari_meta, paste0(filepath,save_folder,"/",CancerType,"_best_features.csv"))
}
  