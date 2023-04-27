library(readxl)
# Define function to calculate metrics for a single cancer type
calculate_metrics <- function(cancer_type, count_matrix, model_name) {
  
  # Extract true labels and predicted labels from count matrix
  tp = count_matrix[1,1]
  tn = count_matrix[2,2]
  fp = count_matrix[2,1]
  fn = count_matrix[1,2]
  
  # Calculate evaluation metrics
  
  # precision
  precision <- tp / (tp + fp)
  
  # recall
  recall <- tp / (tp + fn)
  
  # F1
  f1_score <- 2 * precision * recall / (precision + recall)
  
  # MCC
  numerator <- tp * tn - fp * fn
  denominator <- sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
  mcc = numerator / denominator
  
  # accuracy 
  accuracy = (tp + tn) / (tp + tn + fp + fn)
  
  # Create dataframe of results
  metrics_df <- data.frame(cancer_type = cancer_type,
                           model_name = model_name,
                           F1_Score = f1_score,
                           MCC = mcc,
                           Accuracy = accuracy,
                           Precision = precision,
                           Recall = recall)
  
  return(metrics_df)
}

filepath = "/home/seokwon/nas/"
Cancerlist = dir(paste0(filepath, "/00.data/filtered_TCGA/"))
folder_name = "h2o_bias_pval_exp_dual_cut_50/"

# find best model
# total_result = read_xlsx("~/nas/04.Results/algorithm2_results_for_wo_merge.xlsx", sheet = "Sheet5")
total_result = read_xlsx("~/nas/04.Results/algorithm2_results_for_exp.xlsx")

result_each = as.data.frame(total_result)
colnames(result_each) = result_each[1,]
result_each = result_each[-1,-ncol(result_each)]
result_each$`Test-set` = as.numeric(result_each$`Test-set`)
result_each$`Test-set` = round(result_each$`Test-set`,3)

# 
# result_each = as.data.frame(total_result[2:nrow(total_result), c(7, 10)])
# colnames(result_each) = result_each[1,]
# result_each = result_each[-1,]
# result_each$`Test-set` = as.numeric(result_each$`Test-set`)
# result_each$`Test-set` = round(result_each$`Test-set`,3)
# 
# result_link = as.data.frame(total_result[2:nrow(total_result), c(1, 4)])
# colnames(result_link) = result_link[1,]
# result_link = result_link[-1,]
# result_link$`Test-set` = as.numeric(result_link$`Test-set`)
# result_link$`Test-set` = round(result_link$`Test-set`,3)

metrics_list <- list()
n = 0
# for all
for (num_CancerType in Cancerlist) {
  n = n+1
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  
  # call input
  dual_cut_50 = list.files(paste0(main.path_tc,"/",folder_name), pattern = "model")
  
  tmp_num = result_each[which(result_each$`Cancer Type` == num_CancerType),]$`Test-set`
  tmp_num = as.character(tmp_num)
  best_model_csv <- dual_cut_50[grep(paste0("*",tmp_num,"_confusionmatrix*"), dual_cut_50)]
 
  if ( length(best_model_csv) != 1) {
    tmp_model = result_each[which(result_each$`Cancer Type` == num_CancerType),]$model
    best_model_csv = best_model_csv[grep(paste0(tmp_model,"_*"), best_model_csv)]
  }
  model_type <- sub("_grid_filtered_model.*", "", best_model_csv)
  model_name <- sub(".*_filtered_", "", best_model_csv)
  model_name <- sub("\\_\\d+\\.\\d+\\_confusionmatrix\\.csv", "", model_name)
  
  model_name = paste0(model_type,"_",model_name)
  
  best_model = read.csv(paste0(main.path_tc,"/",folder_name, best_model_csv), skip = 4, row.names = 1)
  rownames(best_model) = c("Alive", "Dead", "Totals")
  colnames(best_model) = c("pre_Alive", "pre_Dead", "Error", "Ratio")
  
  metrics <- calculate_metrics(CancerType , best_model, model_name)
  # Add the metrics to the metrics list
  metrics_list[[n]] <- metrics
  
}

metrics_df <- do.call(rbind, metrics_list)
setwd("~/nas/04.Results/")
write.csv(metrics_df,"ml_indicators_for_exp.csv")
