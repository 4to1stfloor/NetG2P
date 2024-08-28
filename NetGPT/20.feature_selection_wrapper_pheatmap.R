library(pheatmap)
library(ComplexHeatmap)
library(ggplot2)
library(ggrepel)

pathwaylink_01 = readRDS("/home/seokwon/00.data/11.TCGA-STAD/TCGA-STAD_pathwaylink_all_0or1.rds")

pathwaylink_01_filt = pathwaylink_01[,-length(pathwaylink_01)]
pathwaylink_01_filt_wo_zero = pathwaylink_01_filt[,-which(colSums(pathwaylink_01_filt) == 0)]

pheatmap(pathwaylink_01_filt_wo_zero)
hist(colSums(pathwaylink_01_filt_wo_zero))

features = rep(colnames(pathwaylink_01[,-length(pathwaylink_01)]),each = 4)
vitalstatus =  rep(c("Alive","Alive","Dead","Dead"), time = length(colnames(pathwaylink_01[,-length(pathwaylink_01)])))
exist = rep(c("yes","no"), time = length(colnames(pathwaylink_01[,-length(pathwaylink_01)])) *2)
feature_1_count_adjust_cluster = data.frame(features,vitalstatus,exist)

feature_1_count_adjust_cluster$counts = NA

cluster_alive = as.data.frame(colSums(pathwaylink_01[which(pathwaylink_01$vitalstatus == "Alive"),][,-length(pathwaylink_01)]))
colnames(cluster_alive) = "counts"
cluster_dead= as.data.frame(colSums(pathwaylink_01[which(pathwaylink_01$vitalstatus == "Dead"),][,-length(pathwaylink_01)]))
colnames(cluster_dead) = "counts"

feature_1_count_adjust_cluster[which(feature_1_count_adjust_cluster$vitalstatus == "Alive" & feature_1_count_adjust_cluster$exist == "yes"),]$counts = cluster_alive$counts
feature_1_count_adjust_cluster[which(feature_1_count_adjust_cluster$vitalstatus == "Dead" & feature_1_count_adjust_cluster$exist == "yes"),]$counts = cluster_dead$counts

feature_1_count_adjust_cluster[which(feature_1_count_adjust_cluster$vitalstatus == "Alive" & feature_1_count_adjust_cluster$exist == "no"),]$counts = 
  nrow(pathwaylink_01[which(pathwaylink_01$vitalstatus == "Alive"),]) - cluster_alive$counts

feature_1_count_adjust_cluster[which(feature_1_count_adjust_cluster$vitalstatus == "Dead" & feature_1_count_adjust_cluster$exist == "no"),]$counts = 
  nrow(pathwaylink_01[which(pathwaylink_01$vitalstatus == "Dead"),]) - cluster_dead$counts

# fisher test

feature_pval = data.frame(matrix(nrow = length(unique(feature_1_count_adjust_cluster$features))))
rownames(feature_pval) = unique(feature_1_count_adjust_cluster$features)
colnames(feature_pval) = "pval"
if (sum(is.na(feature_1_count_adjust_cluster$count)) == 0) {
  for (features_tmp in unique(feature_1_count_adjust_cluster$features)) {
    tmp_tab = feature_1_count_adjust_cluster[which(feature_1_count_adjust_cluster$features == features_tmp),][,2:4]
    
    result_fisher = fisher.test(xtabs(counts~exist+vitalstatus,data = tmp_tab), workspace = 2e8)
    feature_pval[features_tmp,] =result_fisher$p.value
    remove(tmp_tab)
  }
  
} else {
  print("It has NA")
}

feature_pval_specific = subset(feature_pval, pval < 0.05)
cluster_alive_counts = subset(cluster_alive, subset =  rownames(cluster_alive) %in% rownames(feature_pval_specific))
cluster_dead_counts = subset(cluster_dead, subset =  rownames(cluster_dead) %in% rownames(feature_pval_specific))

feature_pval_specific$alive_counts = cluster_alive_counts$counts
feature_pval_specific$dead_counts = cluster_dead_counts$counts
feature_pval_specific$alive_to_dead_ratio = feature_pval_specific$alive_counts /feature_pval_specific$dead_counts
feature_pval_specific[which(is.infinite(feature_pval_specific$alive_to_dead_ratio)),]$alive_to_dead_ratio = max(feature_pval_specific$alive_to_dead_ratio)

ggplot(data=feature_pval_specific, aes(x=pval,y= alive_to_dead_ratio , label = rownames(feature_pval_specific), )) + 
  geom_point(shape=19) + # shape 19: solid circle
  geom_text_repel() +theme_classic() 

speci_features = rownames(feature_pval_specific)

filtered_pathwaylink = cbind(pathwaylink_01_filt_wo_zero[,speci_features], vitalstatus = pathwaylink_01$vitalstatus)
col_ann = as.data.frame(filtered_pathwaylink[order(filtered_pathwaylink$vitalstatus),]$vitalstatus)
colnames(col_ann) = "vital_status"
rownames(col_ann) = rownames(filtered_pathwaylink[order(filtered_pathwaylink$vitalstatus),])
col_colors = list(vital_status = c("black", "yellow"))

names(col_colors$vital_status) = unique(col_ann$vital_status)
col_ann = subset(filtered_pathwaylink[order(filtered_pathwaylink$vitalstatus),],select = vitalstatus)
filtered_pathwaylink = filtered_pathwaylink[rownames(col_ann),]

ComplexHeatmap::pheatmap(as.matrix(t(filtered_pathwaylink[,-length(filtered_pathwaylink)])), 
                         column_split = col_ann$vitalstatus,
                         labels_col = "",
                         show_rownames = T, 
                         show_colnames = F, 
                         annotation_col = col_ann,
                         annotation_colors = col_colors,
                         # clustering_method = "average",
                         cluster_cols = T,
                         cluster_rows = T
)

ComplexHeatmap::pheatmap(as.matrix(t(filtered_pathwaylink[,-length(filtered_pathwaylink)])), 
                         column_split = col_ann$vitalstatus,
                         labels_col = "",
                         show_rownames = T, 
                         show_colnames = F, 
                         annotation_col = col_ann,
                         annotation_colors = col_colors,
                         # clustering_method = "average",
                         cluster_cols = T,
                         cluster_rows = T,
                         row_km = 6
)

# mldl

h2o.shutdown()
# 0) load library

library(h2o)
set.seed(1)

# 1) call input

filepath = "/home/seokwon/"
ref_path = paste0(filepath, "99.reference/")
Cancerlist = dir(paste0(filepath, "/00.data/"))
Cancerlist = Cancerlist[6]
# can choose!
# type = "each"
type = "link"

bi_num_mode ="1"
# bi_num_mode ="1000000"
# num_CancerType = Cancerlist
for (num_CancerType in Cancerlist) {
  
  # 1) prepare for mldl
  # (1) for TCGA data upload
  
  localH2O = h2o.init(ip = "localhost", port = sample(x=00000:65536,size=1,replace=F), startH2O = TRUE,min_mem_size = "50G",nthreads = 96,enable_assertions = FALSE)
  h2o.removeAll()
  
  Sys.sleep(runif(1,min=1,max=10))
  
  main.path_tc = paste0(filepath, "00.data/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  
  if (type == "link" && bi_num_mode == "1000000") {
    
    data_tc = readRDS(file = paste0(main.path_tc,"/",CancerType,"_pathwaylink_all_1or1000000.rds"))
    # data_tc = data_tc[,-which(dim(data_tc)[1] == colSums(data_tc[,-length(colnames(data_tc))]))]
    
  } else if (type == "link" && bi_num_mode == "1") {
    data_tc = readRDS(file = paste0(main.path_tc,"/",CancerType,"_pathwaylink_all_0or1.rds"))
    # data_tc = data_tc[,-which(dim(data_tc)[1] == colSums(data_tc[,-length(colnames(data_tc))]))]
    
  } else if (type == "each" && bi_num_mode == "1000000") {
    
    data_tc = readRDS(file = paste0(main.path_tc,"/",CancerType,"_pathwayeach_all_1or1000000.rds"))
    # data_tc = data_tc[,-which(dim(data_tc)[1] == colSums(data_tc[,-length(colnames(data_tc))]))]
  } else if (type == "each" && bi_num_mode == "1") {
    data_tc = readRDS(file = paste0(main.path_tc,"/",CancerType,"_pathwayeach_all_0or1.rds"))
    # data_tc = data_tc[,-which(dim(data_tc)[1] == colSums(data_tc[,-length(colnames(data_tc))]))
  }
  
  if ("Not Reported" %in% data_tc$vitalstatus ) {
    if (length(which(data_tc$vitalstatus == "Not Reported")) == 1) {
      data_tc = data_tc[-which(data_tc$vitalstatus == "Not Reported"),]
    } else {
      while ("Not Reported" %in% data_tc$vitalstatus) {
        data_tc = data_tc[-which(data_tc$vitalstatus == "Not Reported")[1],]
      }
      
    }
    
  }
  
  if (sum(colSums(data_tc[,-length(data_tc)]) == 0) != 0 ) {
    data_tc = data_tc[,-which(colSums(data_tc[,-length(data_tc)]) == 0)]
  } else {
    data_tc = data_tc
  }
  data_tc = cbind(data_tc[,speci_features], vitalstatus = data_tc$vitalstatus)
  
  # 2) transfer data for model build
  
  train_val = sample(1:nrow(data_tc), nrow(data_tc)*0.6) 
  testtmp_val = -train_val
  
  # A - data split
  
  train_tc <- data_tc[train_val,]
  testtmp_tc <- data_tc[testtmp_val,]
  valid_val = sample(1:nrow(testtmp_tc), nrow(testtmp_tc)*0.5) 
  
  test_val = -valid_val
  valid_tc <- testtmp_tc[valid_val,]
  test_tc <-  testtmp_tc[test_val,]
  
  # B - Define input (features) and output (response) variables"
  
  features=setdiff(colnames(data_tc),"vitalstatus")
  response = colnames(data_tc)[length(colnames(data_tc))]
  
  # C - Convert the outcome variable into factor
  
  train_tc$vitalstatus <- as.factor(train_tc$vitalstatus)
  valid_tc$vitalstatus <- as.factor(valid_tc$vitalstatus)
  test_tc$vitalstatus <- as.factor(test_tc$vitalstatus)
  
  # D - Convert the data into h2o form (TCGA)
  
  train_tc.hex <- as.h2o(x = train_tc, destination_frame = "train_tc.hex")
  Sys.sleep(runif(1,min=1,max=10))
  test_tc.hex <- as.h2o(x = test_tc, destination_frame = "test_tc.hex")
  Sys.sleep(runif(1,min=1,max=10))
  valid_tc.hex <- as.h2o(x = valid_tc, destination_frame = "valid_tc.hex")
  Sys.sleep(runif(1,min=1,max=10))
  
  # 3) model building 
  # A - Deep learning 
  # a - set hyper-parameter
  
  dl_params <- list(
    activation = c("Rectifier", "Maxout", "Tanh", "RectifierWithDropout", "MaxoutWithDropout", "TanhWithDropout"), 
    hidden = list(c(5, 5, 5), c(10, 10, 10), c(50, 50, 50)),
    epochs = c(50, 100, 200, 500),
    l1 = c(0, 0.00001, 0.0001, 0.001, 0.01), 
    l2 = c(0, 0.00001, 0.0001, 0.001,0.01),
    rate = c(0, 01, 0.005, 0.001),
    rate_annealing = c(1e-8, 1e-7, 1e-6),
    rho = c(0.9, 0.95, 0.99, 0.999),
    epsilon = c(1e-10, 1e-8, 1e-6, 1e-4),
    momentum_start = c(0, 0.5),
    momentum_stable = c(0.99, 0.5, 0),
    input_dropout_ratio = c(0, 0.1, 0.2),
    max_w2 = c(10, 100, 1000, 3.4028235e+38)
  )
  
  search_criteria <- list(strategy = "RandomDiscrete", 
                          max_models = 100,
                          max_runtime_secs = 900,
                          stopping_tolerance = 0.001,
                          stopping_rounds = 15)
  
  # a - train
  
  dl_grid <- h2o.grid(algorithm = "deeplearning", 
                      x = features,
                      y = response,
                      grid_id = "dl_grid",
                      training_frame = train_tc.hex,
                      validation_frame = valid_tc.hex,
                      nfolds = 5,                           
                      fold_assignment = "Stratified",
                      hyper_params = dl_params,
                      search_criteria = search_criteria,
                      seed = 1
  )
  
  Sys.sleep(runif(1,min=1,max=10))
  
  tmp.dl_test_df = as.data.frame(h2o.make_leaderboard(dl_grid, test_tc.hex))
  
  best_dl = h2o.getModel(tmp.dl_test_df[which(tmp.dl_test_df$aucpr ==max(tmp.dl_test_df$aucpr))[1],]$model_id)
  
  best_dl_perf = h2o.performance(best_dl, newdata = test_tc.hex)
  best_dl_conma = as.data.frame(h2o.confusionMatrix(best_dl, test_tc.hex))
  df_dl_merge_confusion = rbind(as.data.frame(h2o.confusionMatrix(best_dl)),"\n",
                                as.data.frame(h2o.confusionMatrix(object = best_dl,test_tc.hex)))
  dir.create(paste0(main.path_tc,"/h2o_filt"))
  write.csv(df_dl_merge_confusion,paste0(main.path_tc,"/h2o_filt/", best_dl@model_id,"_", 1-best_dl_conma$Error[3] ,"_confusionmatrix.csv"))
  h2o.saveModel(best_dl, path = paste0(main.path_tc, "/h2o_filt/", best_dl@model_id,"_", 1-best_dl_conma$Error[3]))
  
  # figure
  
  if (best_dl@model$training_metrics@metrics$AUC == "NaN") {
    print(paste0(num_CancerType," " ,best_dl@algorithm," has no specific model"))
    
  } else {
    
    dir.create(paste0(main.path_tc,"/h2o_filt/dl_figure"))
    fig_path = paste0(main.path_tc,"/h2o_filt/dl_figure")
    setwd(fig_path)
    png(filename = paste0(best_dl@model_id, "roc_test_set.png"),
        width = 2000, height = 2000, units = "px", pointsize = 12,
        bg = "white", res = NA, family = "")
    
    plot(h2o.performance(model = best_dl, newdata = test_tc.hex) , type = "roc")
    dev.off()
    
    png(filename = paste0(best_dl@model_id, "varimp_plot.png"),
        width = 2000, height = 4000, units = "px", pointsize = 12,
        bg = "white", res = NA, family = "")
    
    h2o.varimp_plot(best_dl, num_of_features = 150)
    dev.off()
    
    png(filename = paste0(best_dl@model_id, "gains_lift_plot.png"),
        width = 2000, height = 2000, units = "px", pointsize = 12,
        bg = "white", res = NA, family = "")
    
    h2o.gains_lift_plot(best_dl)
    dev.off()
  } 
  
  remove(tmp.dl_test_df,best_dl,best_dl_perf, best_dl_conma)
  
  # B - RandomForest 
  
  # b - set hyper-parameter
  rf_params <- list(ntrees = seq(50, 500, by = 50),
                    mtries = seq(3, 5, by = 1),
                    max_depth = seq(10, 30, by = 10),
                    min_rows = seq(1, 3, by = 1),
                    nbins = seq(20, 30, by = 10),
                    sample_rate = c(0.55, 0.632, 0.75))
  
  # b - train
  rf_grid <- h2o.grid(algorithm = "randomForest", 
                      x = features,
                      y = response,
                      grid_id = "rf_grid",
                      training_frame = train_tc.hex,
                      validation_frame = valid_tc.hex,
                      nfolds = 5,                           
                      fold_assignment = "Stratified",
                      hyper_params = rf_params,
                      search_criteria = search_criteria,
                      seed = 1
  )
  Sys.sleep(runif(1,min=1,max=10))
  # Now let's evaluate the model performance on a test set
  # so we get an honest estimate of top model performance
  # rf_gridr2 <- h2o.getGrid(grid_id = "rf_grid",
  #                          sort_by = "r2",
  #                          decreasing = TRUE)
  
  tmp.rf_test_df = as.data.frame(h2o.make_leaderboard(rf_grid, test_tc.hex))
  best_rf = h2o.getModel(tmp.rf_test_df[which(tmp.rf_test_df$aucpr ==max(tmp.rf_test_df$aucpr))[1],]$model_id)
  best_rf_perf = h2o.performance(best_rf, newdata = test_tc.hex)
  best_rf_conma = as.data.frame(h2o.confusionMatrix(best_rf, test_tc.hex))
  
  df_rf_merge_confusion = rbind(as.data.frame(h2o.confusionMatrix(best_rf)),"\n",
                                as.data.frame(h2o.confusionMatrix(object = best_rf,test_tc.hex)))
  write.csv(df_rf_merge_confusion,paste0(main.path_tc, "/h2o_filt/", best_rf@model_id,"_", 1-best_rf_conma$Error[3] ,"_confusionmatrix.csv"))
  h2o.saveModel(best_rf, path = paste0(main.path_tc, "/h2o_filt/", best_rf@model_id,"_", 1-best_rf_conma$Error[3]))
  
  # figure
  
  if (best_rf@model$training_metrics@metrics$AUC == "NaN") {
    print(paste0(num_CancerType," " ,best_rf@algorithm," has no specific model"))
    
  } else {
    
    dir.create(paste0(main.path_tc,"/h2o_filt/rf_figure"))
    fig_path = paste0(main.path_tc,"/h2o_filt/rf_figure")
    setwd(fig_path)
    png(filename = paste0(best_rf@model_id, "roc_test_set.png"),
        width = 2000, height = 2000, units = "px", pointsize = 12,
        bg = "white", res = NA, family = "")
    
    plot(h2o.performance(model = best_rf, newdata = test_tc.hex) , type = "roc")
    dev.off()
    
    png(filename = paste0(best_rf@model_id, "varimp_plot.png"),
        width = 2000, height = 4000, units = "px", pointsize = 12,
        bg = "white", res = NA, family = "")
    
    h2o.varimp_plot(best_rf, num_of_features = 150)
    dev.off()
    
    png(filename = paste0(best_rf@model_id, "gains_lift_plot.png"),
        width = 2000, height = 2000, units = "px", pointsize = 12,
        bg = "white", res = NA, family = "")
    
    h2o.gains_lift_plot(best_rf)
    dev.off()
    
  }
  
  remove(tmp.rf_test_df,best_rf,best_rf_perf,best_rf_conma,df_rf_merge_confusion)
  
  # C - Gradient Boosting Machine 
  
  # c - set hyper-parameter
  
  gbm_params <- list(learn_rate = seq(0.01, 1, 0.01),
                     max_depth = seq(2, 30, 1),
                     sample_rate = seq(0, 2, 0.1),
                     col_sample_rate = seq(0.1, 1.0, 0.1))
  
  search_criteria2 <- list(strategy = "RandomDiscrete", max_models = 100, seed = 1)
  
  # c - train
  
  gbm_grid <- h2o.grid("gbm", x = features, y = response,
                       grid_id = "gbm_grid",
                       training_frame = train_tc.hex,
                       validation_frame = valid_tc.hex,
                       ntrees = 100,
                       seed = 1,
                       hyper_params = gbm_params,
                       search_criteria = search_criteria)
  Sys.sleep(runif(1,min=1,max=10))
  # Now let's evaluate the model performance on a test set
  # so we get an honest estimate of top model performance
  # gbm_gridr2 <- h2o.getGrid(grid_id = "gbm_grid",
  #                           sort_by = "r2",
  #                           decreasing = TRUE)
  
  tmp.gbm_test_df = as.data.frame(h2o.make_leaderboard(gbm_grid, test_tc.hex))
  best_gbm = h2o.getModel(tmp.gbm_test_df[which(tmp.gbm_test_df$aucpr ==max(tmp.gbm_test_df$aucpr))[1],]$model_id)
  
  best_gbm_perf = h2o.performance(best_gbm, newdata = test_tc.hex)
  best_gbm_conma = as.data.frame(h2o.confusionMatrix(best_gbm, test_tc.hex))
  plot(best_gbm_perf, type = "roc")
  df_gbm_merge_confusion = rbind(as.data.frame(h2o.confusionMatrix(best_gbm)),"\n",
                                 as.data.frame(h2o.confusionMatrix(object = best_gbm,test_tc.hex)))
  write.csv(df_gbm_merge_confusion,paste0(main.path_tc, "/h2o_filt/", best_gbm@model_id,"_", 1-best_gbm_conma$Error[3] ,"_confusionmatrix.csv"))
  h2o.saveModel(best_gbm, path = paste0(main.path_tc, "/h2o_filt/", best_gbm@model_id,"_", 1-best_gbm_conma$Error[3]))
  
  # figure
  
  if (best_gbm@model$training_metrics@metrics$AUC == "NaN") {
    print(paste0(num_CancerType," " ,best_gbm@algorithm," has no specific model"))
  } else {
    
    dir.create(paste0(main.path_tc,"/h2o_filt/gbm_figure"))
    fig_path = paste0(main.path_tc,"/h2o_filt/gbm_figure")
    setwd(fig_path)
    png(filename = paste0(best_gbm@model_id, "roc_test_set.png"),
        width = 2000, height = 2000, units = "px", pointsize = 12,
        bg = "white", res = NA, family = "")
    
    plot(h2o.performance(model = best_gbm, newdata = test_tc.hex) , type = "roc")
    dev.off()
    
    png(filename = paste0(best_gbm@model_id, "varimp_plot.png"),
        width = 2000, height = 4000, units = "px", pointsize = 12,
        bg = "white", res = NA, family = "")
    
    h2o.varimp_plot(best_gbm, num_of_features = 150)
    dev.off()
    
    png(filename = paste0(best_gbm@model_id, "gains_lift_plot.png"),
        width = 2000, height = 2000, units = "px", pointsize = 12,
        bg = "white", res = NA, family = "")
    
    h2o.gains_lift_plot(best_gbm)
    dev.off()
    
  }
  remove(tmp.bgm_test_df,best_gbm,best_gbm_perf,best_gbm_conma, df_gbm_merge_confusion)
  
  
  # D - eXtreme Gradient Boosting 
  
  # d - set hyper-parameter
  
  xgb_params <- list(ntrees = seq(10, 1000, 1),
                     learn_rate = seq(0.01, 0.5, 0.01),
                     max_depth = seq(1, 10, 1),
                     sample_rate = seq(0.5, 0.8, 0.0001),
                     col_sample_rate = seq(0.2, 1.0, 0.0001),
                     tree_method="hist", 
                     grow_policy="lossguide"
  )
  
  search_criteria2 <- list(strategy = "RandomDiscrete", max_models = 300,
                           stopping_tolerance = 0.001,
                           seed = 1)
  
  # d - train
  
  xgb_grid <- h2o.grid("xgboost", x = features, y = response,
                       grid_id = "xgb_grid",
                       training_frame = train_tc.hex,
                       validation_frame = valid_tc.hex,
                       min_split_improvement = 1e-6,
                       seed = 1,
                       nfolds=5,

                       keep_cross_validation_predictions=TRUE,
                       hyper_params = xgb_params,
                       search_criteria = search_criteria2)
  Sys.sleep(runif(1,min=1,max=10))
  # xgb_gridr2 <- h2o.getGrid(grid_id = "xgb_grid",
  #                           sort_by = "r2",
  #                           decreasing = TRUE)
  
  tmp.xgb_test_df = as.data.frame(h2o.make_leaderboard(xgb_grid, test_tc.hex))
  best_xgb = h2o.getModel(tmp.xgb_test_df[which(tmp.xgb_test_df$aucpr ==max(tmp.xgb_test_df$aucpr))[1],]$model_id)
  best_xgb = h2o.getModel("xgb_grid_model_171")
  best_xgb_perf = h2o.performance(best_xgb, newdata = test_tc.hex)
  best_xgb_conma = as.data.frame(h2o.confusionMatrix(best_xgb, test_tc.hex))
  
  df_xgb_merge_confusion = rbind(as.data.frame(h2o.confusionMatrix(best_xgb)),"\n",
                                 as.data.frame(h2o.confusionMatrix(object = best_xgb,test_tc.hex)))
  write.csv(df_xgb_merge_confusion,paste0(main.path_tc, "/h2o_filt/", best_xgb@model_id,"_", 1-best_xgb_conma$Error[3] ,"_confusionmatrix.csv"))
  h2o.saveModel(best_xgb, path = paste0(main.path_tc, "/h2o_filt/", best_xgb@model_id,"_", 1-best_xgb_conma$Error[3]))
  
  # figure
  
  if (best_xgb@model$training_metrics@metrics$AUC == "NaN") {
    print(paste0(num_CancerType," " ,best_xgb@algorithm," has no specific model"))
    
  } else {
    dir.create(paste0(main.path_tc,"/h2o_filt/xgb_figure"))
    fig_path = paste0(main.path_tc,"/h2o_filt/xgb_figure")
    setwd(fig_path)
    png(filename = paste0(best_xgb@model_id, "roc_test_set.png"),
        width = 2000, height = 2000, units = "px", pointsize = 12,
        bg = "white", res = NA, family = "")
    
    plot(h2o.performance(model = best_xgb, newdata = test_tc.hex) , type = "roc")
    dev.off()
    
    png(filename = paste0(best_xgb@model_id, "varimp_plot.png"),
        width = 2000, height = 4000, units = "px", pointsize = 12,
        bg = "white", res = NA, family = "")
    
    h2o.varimp_plot(best_xgb, num_of_features = 150)
    dev.off()
    
    png(filename = paste0(best_xgb@model_id, "gains_lift_plot.png"),
        width = 2000, height = 2000, units = "px", pointsize = 12,
        bg = "white", res = NA, family = "")
    
    h2o.gains_lift_plot(best_xgb)
    dev.off()
    
  }
  remove(tmp.xgb_test_df,best_xgb,best_xgb_perf,best_xgb_conma, df_xgb_merge_confusion)
  
  # lasso & ridge & elasticnet (glm)
  
  hyper_glm_params <- list(alpha = c(0,0.1,0.25,0.3,0.5,0.7,0.75,0.9,1),
                           lambda = c(seq(0,0.5,.1), 1e-3,1e-5,1e-7,1e-9))
  
  search_criteria2 <- list(strategy = "RandomDiscrete", max_models = 100 ,seed = 1)
  
  glm_grid <- h2o.grid("glm", x = features, y = response,
                       grid_id = "glm_grid",
                       training_frame = train_tc.hex,
                       validation_frame = valid_tc.hex,
                       seed = 1,
                       nfolds=5,
                       family = c("binomial"),
                       keep_cross_validation_predictions=TRUE,
                       hyper_params = hyper_glm_params,
                       search_criteria = search_criteria2)
  Sys.sleep(runif(1,min=1,max=10))
  # glm_gridr2 <- h2o.getGrid(grid_id = "glm_grid",
  #                           sort_by = "r2",
  #                           decreasing = TRUE)
  
  tmp.glm_test_df = as.data.frame(h2o.make_leaderboard(glm_grid, test_tc.hex))
  best_glm = h2o.getModel(tmp.glm_test_df[which(tmp.glm_test_df$mse ==min(tmp.glm_test_df$mse))[1],]$model_id)
  
  best_glm_perf = h2o.performance(best_glm, newdata = test_tc.hex)
  best_glm_conma = as.data.frame(h2o.confusionMatrix(best_glm, test_tc.hex))
  
  df_glm_merge_confusion = rbind(as.data.frame(h2o.confusionMatrix(best_glm)),"\n",
                                 as.data.frame(h2o.confusionMatrix(object = best_glm,test_tc.hex)))
  write.csv(df_glm_merge_confusion,paste0(main.path_tc, "/h2o_filt/", best_glm@model_id,"_", 1-best_glm_conma$Error[3] ,"_confusionmatrix.csv"))
  h2o.saveModel(best_glm, path = paste0(main.path_tc, "/h2o_filt/", best_glm@model_id,"_", 1-best_glm_conma$Error[3]))
  
  
  # figure
  
  if (best_glm@model$training_metrics@metrics$AUC == "NaN") {
    
    print(paste0(num_CancerType," " ,best_glm@algorithm," has no specific model"))
    
  } else {
    
    dir.create(paste0(main.path_tc,"/h2o_filt/glm_figure"))
    fig_path = paste0(main.path_tc,"/h2o_filt/glm_figure")
    setwd(fig_path)
    png(filename = paste0(best_glm@model_id, "roc_test_set.png"),
        width = 2000, height = 2000, units = "px", pointsize = 12,
        bg = "white", res = NA, family = "")
    
    plot(h2o.performance(model = best_glm, newdata = test_tc.hex) , type = "roc")
    dev.off()
    
    png(filename = paste0(best_glm@model_id, "varimp_plot.png"),
        width = 2000, height = 4000, units = "px", pointsize = 12,
        bg = "white", res = NA, family = "")
    
    try(h2o.varimp_plot(best_glm, num_of_features = 150))
    
    dev.off()
    
    png(filename = paste0(best_glm@model_id, "gains_lift_plot.png"),
        width = 2000, height = 2000, units = "px", pointsize = 12,
        bg = "white", res = NA, family = "")
    
    try(h2o.gains_lift_plot(best_glm)) 
    
    dev.off()
    
  }
  remove(tmp.glm_test_df,best_glm,best_glm_perf,best_glm_conma, df_glm_merge_confusion)
  
  h2o.shutdown(prompt = F)
  
}



### fin

# 
# cancer_complexheat = ComplexHeatmap::pheatmap(as.matrix(t(filtered_pathwaylink[,-length(filtered_pathwaylink)])), 
#                                               column_split = col_ann$vitalstatus,
#                                               labels_col = "",
#                                               show_rownames = T, 
#                                               show_colnames = F, 
#                                               annotation_col = col_ann,
#                                               annotation_colors = col_colors,
#                                               # clustering_method = "average",
#                                               cluster_cols = T,
#                                               cluster_rows = T,
#                                               row_km = 10)
# 
# draw_ht = draw(cancer_complexheat)
# n = 0
# for (cluster in row_order(draw_ht)) {
#   # print(cluster)
#   n = n +1
#   filtered_features = colnames(as.data.frame(pathwaylink_01_filt_wo_zero[rownames(col_ann),])[,cluster])
#   tmp_cutrow = cbind(pathwaylink_01[,filtered_features],vitalstatus = pathwaylink_01$vitalstatus)
#   tmp_cutrow_order = tmp_cutrow[order(tmp_cutrow$vitalstatus),]
#   tmp_cutrow_order_alive = tmp_cutrow_order[which(tmp_cutrow_order$vitalstatus == "Alive"),]
#   tmp_cutrow_order_dead = tmp_cutrow_order[which(tmp_cutrow_order$vitalstatus == "Dead"),]
#   tmp_cutrow_order_alive_filt = tmp_cutrow_order_alive[,-length(tmp_cutrow_order_alive)]
#   tmp_cutrow_order_dead_filt = tmp_cutrow_order_dead[,-length(tmp_cutrow_order_dead)]
#   if (t.test(colMeans(tmp_cutrow_order_alive_filt), colMeans(tmp_cutrow_order_dead_filt))$p.value < 0.05) {
#     print(names(row_order(draw_ht))[n])
#     print(t.test(colMeans(tmp_cutrow_order_alive_filt), colMeans(tmp_cutrow_order_dead_filt))$p.value)
#   }
#   remove(tmp_cutrow_order_alive_filt,tmp_cutrow_order_dead_filt)
# }
# colnames(pathwaylink_01_filt)[row_order(draw_ht)$'2']
# pathwaylink_01
