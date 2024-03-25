
library(tidyverse)
library(readxl)
library(h2o)

filepath = "/home/seokwon/nas/"
ref_path = paste0(filepath, "99.reference/")
Cancerlist = dir(paste0(filepath, "/00.data/filtered_TCGA/"))
# Cancerlist = Cancerlist[7:9]
Cancerlist = Cancerlist[c(-11,-12)]

localH2O = h2o.init(ip = "localhost", port = sample(x=00000:65536,size=1,replace=F), 
                    startH2O = TRUE,min_mem_size = "400G",
                    nthreads = 96,enable_assertions = FALSE)

find_max_model <- function(vec) {
  vec <- vec[!grepl(".csv", vec)]  # remove elements with ".csv" extension
  max_model <- vec[which.max(as.numeric(str_sub(vec, -5,-1) ))]
  return(max_model)
}

# num_CancerType =  "04.TCGA-CESC"

fig_path = paste0(filepath,"04.Results/cell/")
if(!dir.exists(fig_path)){
  dir.create(fig_path)
  print(paste0("Created folder: ", fig_path))
} else {
  print(paste0("Folder already exists: ", fig_path))
}
setwd(fig_path)

for (num_CancerType in Cancerlist) {
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType, "/", folder_name)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  Cancername = gsub('TCGA-', '' ,CancerType)
  dual_cut_50 = list.files(paste0(main.path_tc), pattern = "model")
  best_model_folder = find_max_model(dual_cut_50)
  best_model = h2o.loadModel(paste0(main.path_tc,"/",best_model_folder,"/", gsub('.{6}$', '', best_model_folder)))
  
  cell = readRDS(paste0(filepath, "00.data/filtered_TCGA/", num_CancerType, "/",Cancername, "_cellline_dual_all_log.rds"))
  cell.hex <- as.h2o(x = cell, destination_frame = "cell.hex")
  
  # rownames(cell)
  # h2o.performance
  tmp_predict = as.data.frame(h2o.predict(best_model , newdata = cell.hex))
  rownames(tmp_predict) = rownames(cell)
  tmp_predict_filt = tmp_predict %>% 
    mutate(adjust_predic = ifelse(Alive > Dead , "Alive", "Dead"),
           cellline_id = rownames(tmp_predict)) %>% 
    select(-predict) %>% 
    select(cellline_id, adjust_predic, Alive, Dead)
  
  write_csv(tmp_predict_filt, paste0(fig_path,Cancername,"_cellline_predict_for_bestmodel.csv"),)
}

h2o.shutdown(prompt = F)
