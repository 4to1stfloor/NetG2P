library(tidyverse)
library(readxl)
library(h2o)

saved_model_folder = "saved folder before step - your path"
save_folder = "your path"
localH2O = h2o.init(ip = "localhost", port = sample(x=00000:65536,size=1,replace=F), 
                    startH2O = TRUE,min_mem_size = "400G",
                    nthreads = 96,enable_assertions = FALSE)

find_max_model <- function(vec) {
  vec <- vec[!grepl(".csv", vec)]  # remove elements with ".csv" extension
  max_model <- vec[which.max(as.numeric(str_sub(vec, -5,-1) ))] # you can edit this part 
  return(max_model)
}

CancerType = "your instrested cancer type"
find_list = list.files(saved_model_folder, pattern = "model")
best_model_folder = find_max_model(find_list)
best_model = h2o.loadModel(paste0(main.path_tc,"/",best_model_folder,"/", gsub('.{6}$', '', best_model_folder)))
best_model_vari_meta = as.data.frame(h2o.varimp(best_model))

write_csv(best_model_vari_meta, paste0(save_folder,"/",CancerType,"_best_features.csv"))

