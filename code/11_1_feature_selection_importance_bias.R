library(pheatmap)
library(ComplexHeatmap)
library(ggplot2)
library(ggrepel)


# mldl

h2o.shutdown()
# 0) load library

library(h2o)
set.seed(1)

# 1) call input

filepath = "/home/seokwon/"
ref_path = paste0(filepath, "99.reference/")
Cancerlist = dir(paste0(filepath, "/00.data/"))
Cancerlist = Cancerlist[c(-2,-3)]

# can choose!
# type = "each"
type = "link"

bi_num_mode ="1"
# bi_num_mode ="1000000"
# num_CancerType = Cancerlist
folder_name = "h2o_filtered_bias"
for (num_CancerType in Cancerlist) {
  
  
  # 1) prepare for mldl
  # (1) for TCGA data upload
  
  localH2O = h2o.init(ip = "localhost", port = sample(x=00000:65536,size=1,replace=F), startH2O = TRUE,min_mem_size = "200G",nthreads = 96, enable_assertions = FALSE)
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
  
  if (sum((colSums(data_tc[,-ncol(data_tc)]) == 0)) != 0 ) {
    data_tc_ori = data_tc[,-which(colSums(data_tc[,-ncol(data_tc)]) == 0)]
  } else {
    data_tc_ori = data_tc
  }
  
  # A. data_split_with_bias
  
  data_tc_ori_alive = data_tc_ori[data_tc_ori$vitalstatus == "Alive",]
  data_tc_ori_dead = data_tc_ori[data_tc_ori$vitalstatus == "Dead",]
  
  train_alive = sample(1:nrow(data_tc_ori_alive), nrow(data_tc_ori_alive)*0.6) 
  train_dead = sample(1:nrow(data_tc_ori_dead), nrow(data_tc_ori_dead)*0.6) 
  testtmp_val_alive = -train_alive
  testtmp_val_dead = -train_dead
  
  train_tc_alive = data_tc_ori_alive[train_alive,]
  train_tc_dead = data_tc_ori_dead[train_dead,]
  
  train_tc = rbind(train_tc_alive,train_tc_dead)
  
  testtmp_tc_alive = data_tc_ori_alive[testtmp_val_alive,]
  testtmp_tc_dead = data_tc_ori_dead[testtmp_val_dead,]
  
  valid_val_alive = sample(1:nrow(testtmp_tc_alive), nrow(testtmp_tc_alive)*0.5) 
  valid_val_dead = sample(1:nrow(testtmp_tc_dead), nrow(testtmp_tc_dead)*0.5) 
  
  test_val_alive = -valid_val_alive
  test_val_dead = -valid_val_dead
  
  valid_tc_alive = testtmp_tc_alive[valid_val_alive,]
  valid_tc_dead = testtmp_tc_dead[valid_val_dead,]
  valid_tc = rbind(valid_tc_alive,valid_tc_dead)
  
  test_tc_alive = testtmp_tc_alive[test_val_alive,]
  test_tc_dead = testtmp_tc_dead[test_val_dead,]
  test_tc = rbind(test_tc_alive,test_tc_dead)
  
  # B - Define input (features) and output (response) variables"
  
  features=setdiff(colnames(data_tc_ori),"vitalstatus")
  response = colnames(data_tc_ori)[length(colnames(data_tc_ori))]
  # response = "status"
  # C - Convert the outcome variable into factor
  
  train_tc$vitalstatus <- as.factor(train_tc$vitalstatus)
  valid_tc$vitalstatus <- as.factor(valid_tc$vitalstatus)
  test_tc$vitalstatus <- as.factor(test_tc$vitalstatus)

  # train_tc$status <- as.factor(train_tc$status)
  # valid_tc$status <- as.factor(valid_tc$status)
  # test_tc$status <- as.factor(test_tc$status)
  # 
  # D - Convert the data into h2o form (TCGA)
  
  train_tc.hex <- as.h2o(x = train_tc, destination_frame = "train_tc.hex")
  Sys.sleep(runif(1,min=1,max=10))
  test_tc.hex <- as.h2o(x = test_tc, destination_frame = "test_tc.hex")
  Sys.sleep(runif(1,min=1,max=10))
  valid_tc.hex <- as.h2o(x = valid_tc, destination_frame = "valid_tc.hex")
  Sys.sleep(runif(1,min=1,max=10))
  
  # 3) model building 
  # test 
  
  # aml <- h2o.automl(x = features, y = response,
  #                   training_frame = train_tc.hex,
  #                   validation_frame = valid_tc.hex,
  #                   leaderboard_frame = test_tc.hex,
  #                   nfolds = 5,
  #                   distribution = "AUTO",
  #                   sort_metric = "AUTO",
  #                   keep_cross_validation_fold_assignment = TRUE,
  #                   keep_cross_validation_models = TRUE,
  #                   keep_cross_validation_predictions = TRUE,
  #                   max_runtime_secs = 3000)
  # lb <- aml@leaderboard
  # lb
  # A - Deep learning 
  # a - set hyper-parameter
  
  dl_params <- list(
    activation = c("Rectifier", "Maxout", "Tanh", "RectifierWithDropout", "MaxoutWithDropout", "TanhWithDropout"), 
    hidden = list(c(5, 5, 5), c(10, 10, 10), c(50, 50, 50), c(100,100,100)),
    epochs = c(50, 100, 200, 500, 1000, 2000,5000),
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
                          max_models = 300,
                          max_runtime_secs = 1000,
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
                      fold_assignment = "Modulo",
                      keep_cross_validation_predictions=TRUE,
                      export_weights_and_biases=T,
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
  dir.create(paste0(main.path_tc,"/",folder_name))
  write.csv(df_dl_merge_confusion,paste0(main.path_tc,"/",folder_name,"/", best_dl@model_id,"_", round(1-best_dl_conma$Error[3],3) ,"_confusionmatrix.csv"))
  h2o.saveModel(best_dl, path = paste0(main.path_tc,"/",folder_name,"/", best_dl@model_id,"_", round(1-best_dl_conma$Error[3],3) ))
  
  # figure
  
  if (best_dl@model$training_metrics@metrics$AUC == "NaN") {
    print(paste0(num_CancerType," " ,best_dl@algorithm," has no specific model"))
    
  } else {
    
    dir.create(paste0(main.path_tc,"/",folder_name,"/dl_figure"))
    fig_path = paste0(main.path_tc,"/",folder_name,"/dl_figure")
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
  
  remove(tmp.dl_test_df,best_dl_perf, best_dl_conma)
  
  # B - Gradient Boosting Machine 
  
  # b - set hyper-parameter
  
  gbm_params <- list(learn_rate = seq(0.01, 1, 0.01),
                     max_depth = seq(2, 30, 1),
                     ntrees = seq(10, 5000, 10),
                     sample_rate = seq(0, 2, 0.01),
                     col_sample_rate = seq(0.2, 2.0 ,0.01),
                     col_sample_rate_change_per_level = seq(0.9,1.1,0.01),  
                     nbins = round(2 ^ seq(2, 6, length = 15)),
                     min_split_improvement = c(0,1e-8,1e-6,1e-4),
                     histogram_type = c("UniformAdaptive","Random","QuantilesGlobal","RoundRobin")  )
  
  search_criteria <- list(strategy = "RandomDiscrete",
                           max_runtime_secs = 60* 60, 
                           max_models = 100,  
                           stopping_rounds = 5,                
                           stopping_metric = "logloss",
                           stopping_tolerance = 1e-5,seed = 1)
  
  # b - train
  
  gbm_grid <- h2o.grid("gbm", x = features, y = response,
                       grid_id = "gbm_grid",
                       training_frame = train_tc.hex,
                       validation_frame = valid_tc.hex,
                       fold_assignment = "Modulo",
                       nfolds = 5,
                       seed = 1,
                       keep_cross_validation_predictions=TRUE,
                       hyper_params = gbm_params,
                       search_criteria = search_criteria)
  
  Sys.sleep(runif(1,min=1,max=10))
  
  # Now let's evaluate the model performance on a test set

  ##
  tmp.gbm_test_df = as.data.frame(h2o.make_leaderboard(gbm_grid, test_tc.hex))
  best_gbm = h2o.getModel(tmp.gbm_test_df[which(tmp.gbm_test_df$aucpr ==max(tmp.gbm_test_df$aucpr))[1],]$model_id)
  
  best_gbm_perf = h2o.performance(best_gbm, newdata = test_tc.hex)
  best_gbm_conma = as.data.frame(h2o.confusionMatrix(best_gbm, test_tc.hex))

  df_gbm_merge_confusion = rbind(as.data.frame(h2o.confusionMatrix(best_gbm)),"\n",
                                 as.data.frame(h2o.confusionMatrix(object = best_gbm,test_tc.hex)))
  write.csv(df_gbm_merge_confusion,paste0(main.path_tc, "/",folder_name,"/", best_gbm@model_id,"_", round(1-best_gbm_conma$Error[3],3) ,"_confusionmatrix.csv"))
  h2o.saveModel(best_gbm, path = paste0(main.path_tc, "/",folder_name,"/", best_gbm@model_id,"_", round(1-best_gbm_conma$Error[3],3)))
  
  # figure
  
  if (best_gbm@model$training_metrics@metrics$AUC == "NaN") {
    print(paste0(num_CancerType," " ,best_gbm@algorithm," has no specific model"))
  } else {
    
    dir.create(paste0(main.path_tc,"/",folder_name,"/gbm_figure"))
    fig_path = paste0(main.path_tc,"/",folder_name,"/gbm_figure")
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
  remove(tmp.gbm_test_df,best_gbm_perf,best_gbm_conma, df_gbm_merge_confusion)
  
  
  # D - eXtreme Gradient Boosting 
  
  # d - set hyper-parameter
  
  xgb_params <- list(ntrees = seq(10, 5000, 10),
                     learn_rate = seq(0.01, 1, 0.01),
                     max_depth = seq(2, 30, 1),
                     sample_rate = seq(0, 2, 0.01),
                     tree_method="hist", 
                     grow_policy="lossguide"
  )
  
  search_criteria <- list(strategy = "RandomDiscrete",
                          max_runtime_secs = 60* 60, 
                          max_models = 100,  
                          stopping_rounds = 5,                
                          stopping_metric = "logloss",
                          stopping_tolerance = 1e-5,seed = 1)
  
  # d - train
  
  xgb_grid <- h2o.grid("xgboost", x = features, y = response,
                       grid_id = "xgb_grid",
                       training_frame = train_tc.hex,
                       validation_frame = valid_tc.hex,
                       min_split_improvement = 1e-6,
                       fold_assignment = "Modulo",
                       seed = 1,
                       nfolds=5,
                       keep_cross_validation_predictions=TRUE,
                       hyper_params = xgb_params,
                       search_criteria = search_criteria)
  
  Sys.sleep(runif(1,min=1,max=10))
  # xgb_gridr2 <- h2o.getGrid(grid_id = "xgb_grid",
  #                           sort_by = "r2",
  #                           decreasing = TRUE)
  
  tmp.xgb_test_df = as.data.frame(h2o.make_leaderboard(xgb_grid, test_tc.hex))
  best_xgb = h2o.getModel(tmp.xgb_test_df[which(tmp.xgb_test_df$aucpr ==max(tmp.xgb_test_df$aucpr))[1],]$model_id)

  best_xgb_perf = h2o.performance(best_xgb, newdata = test_tc.hex)
  best_xgb_conma = as.data.frame(h2o.confusionMatrix(best_xgb, test_tc.hex))
  
  df_xgb_merge_confusion = rbind(as.data.frame(h2o.confusionMatrix(best_xgb)),"\n",
                                 as.data.frame(h2o.confusionMatrix(object = best_xgb,test_tc.hex)))
  write.csv(df_xgb_merge_confusion,paste0(main.path_tc, "/",folder_name,"/", best_xgb@model_id,"_", round(1-best_xgb_conma$Error[3],3) ,"_confusionmatrix.csv"))
  h2o.saveModel(best_xgb, path = paste0(main.path_tc, "/",folder_name,"/", best_xgb@model_id,"_", round(1-best_xgb_conma$Error[3],3)))
  
  # figure
  
  if (best_xgb@model$training_metrics@metrics$AUC == "NaN") {
    print(paste0(num_CancerType," " ,best_xgb@algorithm," has no specific model"))
    
  } else {
    dir.create(paste0(main.path_tc,"/",folder_name,"/xgb_figure"))
    fig_path = paste0(main.path_tc,"/",folder_name,"/xgb_figure")
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
  remove(tmp.xgb_test_df,best_xgb_perf,best_xgb_conma, df_xgb_merge_confusion)
  
  # lasso & ridge & elasticnet (glm)
  
  hyper_glm_params <- list(alpha = c(0,0.1,0.25,0.3,0.5,0.7,0.75,0.9,1),
                           lambda = c(seq(0,0.5,.1), 1e-3,1e-5,1e-7,1e-9))
  
  search_criteria2 <- list(strategy = "RandomDiscrete", max_models = 100 ,seed = 1)
  
  glm_grid <- h2o.grid("glm", x = features, y = response,
                       grid_id = "glm_grid",
                       training_frame = train_tc.hex,
                       validation_frame = valid_tc.hex,
                       fold_assignment = "Modulo",
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
  best_glm = h2o.getModel(tmp.glm_test_df[which(tmp.glm_test_df$aucpr ==max(tmp.glm_test_df$aucpr))[1],]$model_id)
  
  best_glm_perf = h2o.performance(best_glm, newdata = test_tc.hex)
  best_glm_conma = as.data.frame(h2o.confusionMatrix(best_glm, test_tc.hex))
  
  df_glm_merge_confusion = rbind(as.data.frame(h2o.confusionMatrix(best_glm)),"\n",
                                 as.data.frame(h2o.confusionMatrix(object = best_glm,test_tc.hex)))
  write.csv(df_glm_merge_confusion,paste0(main.path_tc, "/",folder_name,"/", best_glm@model_id,"_", round(1-best_glm_conma$Error[3],3) ,"_confusionmatrix.csv"))
  h2o.saveModel(best_glm, path = paste0(main.path_tc, "/",folder_name,"/", best_glm@model_id,"_", round(1-best_glm_conma$Error[3], 3)))
  
  
  # figure
  
  if (best_glm@model$training_metrics@metrics$AUC == "NaN") {
    
    print(paste0(num_CancerType," " ,best_glm@algorithm," has no specific model"))
    
  } else {
    
    dir.create(paste0(main.path_tc,"/",folder_name,"/glm_figure"))
    fig_path = paste0(main.path_tc,"/",folder_name,"/glm_figure")
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
  remove(tmp.glm_test_df,best_glm_perf,best_glm_conma, df_glm_merge_confusion)
  # best_dl = h2o.loadModel(paste0(main.path_tc, "/",folder_name,"/","dl_grid_model_105_0.764/dl_grid_model_105"))
  # best_gbm = h2o.loadModel(paste0(main.path_tc, "/",folder_name,"/","gbm_grid_model_251_1/gbm_grid_model_251"))
  # best_xgb = h2o.loadModel(paste0(main.path_tc, "/",folder_name,"/","xgb_grid_model_2_1/xgb_grid_model_2"))
  # best_glm = h2o.loadModel(paste0(main.path_tc, "/",folder_name,"/","glm_grid_model_41_0.764/glm_grid_model_41"))
  # 
  # stack ensemble
  ensemble <- h2o.stackedEnsemble(x = features,
                                  y = response,
                                  training_frame = train_tc.hex,
                                  validation_frame = valid_tc.hex,
                                  metalearner_algorithm = "AUTO",
                                  base_models = list(best_dl,best_gbm,best_xgb,best_glm))
                                  

  
  # tmp.ensemble_test_df = as.data.frame(h2o.make_leaderboard(ensemble, test_tc.hex))
  # best_ensemble = h2o.getModel(tmp.ensemble_test_df[which(tmp.ensemble_test_df$aucpr == max(tmp.ensemble_test_df$aucpr))[1],]$model_id)
  
  best_ensemble_perf = h2o.performance(ensemble, newdata = test_tc.hex)
  ensemble_auc_test <- h2o.auc(best_ensemble_perf)
  print(sprintf("Ensemble Test AUC:  %s", ensemble_auc_test))
  
  best_ensemble_conma = as.data.frame(h2o.confusionMatrix(ensemble, test_tc.hex))
  df_ensemble_merge_confusion = rbind(as.data.frame(h2o.confusionMatrix(ensemble)),"\n",
                                      as.data.frame(h2o.confusionMatrix(object = ensemble,test_tc.hex)))
  write.csv(df_ensemble_merge_confusion,paste0(main.path_tc, "/",folder_name,"/ensemble_", round(ensemble_auc_test,3) ,"_confusionmatrix.csv"))
  h2o.saveModel(ensemble, path = paste0(main.path_tc, "/",folder_name,"/ensemble_", round(ensemble_auc_test,3) ))
  
  h2o.shutdown(prompt = F)
  
}



### fin