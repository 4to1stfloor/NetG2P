h2o.shutdown()
# 0) load library

library(h2o)
set.seed(1)

# 1) call input

filepath = "/home/seokwon/"
ref_path = paste0(filepath, "99.reference/")

Cancerlist = dir(paste0(filepath, "/00.data/"))
Cancerlist = Cancerlist[c(-1,-3,-7,-20,-21,-34,-35)]
# Cancerlist= Cancerlist[seq(-1,-19)]
Cancerlist= Cancerlist[-c(1:10)]
Cancerlist = Cancerlist[c(1:8)]
Cancerlist = Cancerlist[13]
# can choose!
type = "each"
# type = "link"

bi_num_mode ="1"
# bi_num_mode ="1000000"
num_CancerType = Cancerlist
for (num_CancerType in Cancerlist) {
  
  # 1) prepare for mldl
  # (1) for TCGA data upload
  
  localH2O = h2o.init(ip = "localhost", port = sample(x=00000:65536,size=1,replace=F), startH2O = TRUE,min_mem_size = "50G",nthreads = 96,enable_assertions = FALSE)
  
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
    epochs = c(50, 100),
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
  dir.create(paste0(main.path_tc,"/h2o_each"))
  write.csv(df_dl_merge_confusion,paste0(main.path_tc,"/h2o_each/", best_dl@model_id,"_", 1-best_dl_conma$Error[3] ,"_confusionmatrix.csv"))
  h2o.saveModel(best_dl, path = paste0(main.path_tc, "/h2o_each/", best_dl@model_id,"_", 1-best_dl_conma$Error[3]))
  
  # figure
  
  if (best_dl@model$training_metrics@metrics$AUC == "NaN") {
    print(paste0(num_CancerType," " ,best_dl@algorithm," has no specific model"))
    
    } else {
      
      dir.create(paste0(main.path_tc,"/h2o_each/dl_figure"))
      fig_path = paste0(main.path_tc,"/h2o_each/dl_figure")
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
  write.csv(df_rf_merge_confusion,paste0(main.path_tc, "/h2o_each/", best_rf@model_id,"_", 1-best_rf_conma$Error[3] ,"_confusionmatrix.csv"))
  h2o.saveModel(best_rf, path = paste0(main.path_tc, "/h2o_each/", best_rf@model_id,"_", 1-best_rf_conma$Error[3]))
  
  # figure
  
  if (best_rf@model$training_metrics@metrics$AUC == "NaN") {
    print(paste0(num_CancerType," " ,best_rf@algorithm," has no specific model"))
    
    } else {
    
      dir.create(paste0(main.path_tc,"/h2o_each/rf_figure"))
      fig_path = paste0(main.path_tc,"/h2o_each/rf_figure")
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
  write.csv(df_gbm_merge_confusion,paste0(main.path_tc, "/h2o_each/", best_gbm@model_id,"_", 1-best_gbm_conma$Error[3] ,"_confusionmatrix.csv"))
  h2o.saveModel(best_gbm, path = paste0(main.path_tc, "/h2o_each/", best_gbm@model_id,"_", 1-best_gbm_conma$Error[3]))
  
  # figure
  
  if (best_gbm@model$training_metrics@metrics$AUC == "NaN") {
    print(paste0(num_CancerType," " ,best_gbm@algorithm," has no specific model"))
  } else {
    
    dir.create(paste0(main.path_tc,"/h2o_each/gbm_figure"))
    fig_path = paste0(main.path_tc,"/h2o_each/gbm_figure")
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
                     dmatrix_type = c("dense")
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
  
  best_xgb_perf = h2o.performance(best_xgb, newdata = test_tc.hex)
  best_xgb_conma = as.data.frame(h2o.confusionMatrix(best_xgb, test_tc.hex))
  
  df_xgb_merge_confusion = rbind(as.data.frame(h2o.confusionMatrix(best_xgb)),"\n",
                                 as.data.frame(h2o.confusionMatrix(object = best_xgb,test_tc.hex)))
  write.csv(df_xgb_merge_confusion,paste0(main.path_tc, "/h2o_each/", best_xgb@model_id,"_", 1-best_xgb_conma$Error[3] ,"_confusionmatrix.csv"))
  h2o.saveModel(best_xgb, path = paste0(main.path_tc, "/h2o_each/", best_xgb@model_id,"_", 1-best_xgb_conma$Error[3]))
  
  # figure
  
  if (best_xgb@model$training_metrics@metrics$AUC == "NaN") {
    print(paste0(num_CancerType," " ,best_xgb@algorithm," has no specific model"))
    
    } else {
      dir.create(paste0(main.path_tc,"/h2o_each/xgb_figure"))
      fig_path = paste0(main.path_tc,"/h2o_each/xgb_figure")
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
  write.csv(df_glm_merge_confusion,paste0(main.path_tc, "/h2o_each/", best_glm@model_id,"_", 1-best_glm_conma$Error[3] ,"_confusionmatrix.csv"))
  h2o.saveModel(best_glm, path = paste0(main.path_tc, "/h2o_each/", best_glm@model_id,"_", 1-best_glm_conma$Error[3]))
  
  
  # figure
  
  if (best_glm@model$training_metrics@metrics$AUC == "NaN") {
    
    print(paste0(num_CancerType," " ,best_glm@algorithm," has no specific model"))
    
    } else {
    
      dir.create(paste0(main.path_tc,"/h2o_each/glm_figure"))
      fig_path = paste0(main.path_tc,"/h2o_each/glm_figure")
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


# library(mlr)
# 
# task=makeClassifTask(data=test_tc,target="vitalstatus")
# 
# dummylrner=makeLearner("classif.rpart",predict.type ="prob")
# dummymod=train(dummylrner,task)
# dummypred=predict(dummymod,task)
# 
# LASSOPRED=predict(lassomod,test_tc.hex)%>%as.data.frame()
# LASSOPRED$cal_predic = ifelse(LASSOPRED$Alive > LASSOPRED$Dead, "Alive", "Dead")
# LASSOPRED=LASSOPRED%>%mutate(.,id=c(1:length(rownames(test_tc))),truth=test_tc$vitalstatus,prob.Alive=.[,2],prob.Dead=.[,3],response=.$cal_predic)%>%.[,-c(1:3)]
# 
# 
# RIDGEPRED=predict(ridgemod,test_tc.hex)%>%as.data.frame()
# RIDGEPRED=RIDGEPRED%>%mutate(.,id=c(1:length(rownames(test_tc))),truth=test_tc$vitalstatus,prob.Alive=.[,2],prob.Dead=.[,3],response=.$predict)%>%.[,-c(1:3)]
# remove(ENETPRED)
# remove(enetclass)
# enetclass=predict(enet,as.matrix(test_tc[,-length(colnames(test_tc))]),type="class",s = "lambda.min")%>%as.data.frame()
# ENETPRED=predict(enet,as.matrix(test_tc[,-length(colnames(test_tc))]) ,s = "lambda.min")%>%as.data.frame()%>%mutate(.,id=c(1:length(rownames(test_tc))),truth=test_tc$vitalstatus,prob.Alive=1-.[,1],prob.Dead=.[,1],response=enetclass[,1])%>%.[,-1]
# 
# RFPRED=predict(rfmod,test_tc.hex)%>%as.data.frame()
# RFPRED$cal_predic = ifelse(RFPRED$Alive > RFPRED$Dead, "Alive", "Dead")
# RFPRED=RFPRED%>%mutate(.,id=c(1:length(rownames(test_tc))),truth=test_tc$vitalstatus,prob.Alive=.[,2],prob.Dead=.[,3],response=.$cal_predic)%>%.[,-c(1:3)]
# 
# pred1=dummypred
# pred1$data=LASSOPRED
# pred2=dummypred
# pred2$data=ENETPRED
# pred3=dummypred
# pred3$data=RFPRED
# pred4=dummypred
# pred4$data=RIDGEPRED%>%na.omit()
# 
# measure=list(mlr::bac,mlr::auc,mlr::acc,mlr::f1,mlr::ppv,mlr::npv,mlr::tnr,mlr::tpr)
# p1 = NULL
# p2 = NULL
# p3 = NULL
# p4 = NULL
# 
# p1=performance(pred1,measure)
# p2=performance(pred2,measure)
# p3=performance(pred3,measure)
# p4=performance(pred4,measure)
# 
# lasso=p1%>%as.vector()
# enet=p2%>%as.vector()
# rf=p3%>%as.vector()
# rdg=p4%>%as.vector()
# 
# labs=c("BAC","AUC","ACC","F1","PPV","NPV","TPR","TNR")
# 
# scores=list("Lasso"=lasso*10,"Elasticnet"=enet*10,"RandomForest"=rf*10,"Ridge"=rdg*10)
# scores=list("Elasticnet"=enet*10,"RandomForest"=rf*10,"Ridge"=rdg*10)
# scores=list("Lasso"=lasso*10,"RandomForest"=rf*10,"Ridge"=rdg*10)
# library(radarchart)
# 
# chartJSRadar(scores = scores, labs = labs, maxScale = 10)
# 
# measure2=list(ber,brier,fnr,fpr,mmce,fdr)
# 
# e1 = NULL
# e2 = NULL
# e3 = NULL
# e4 = NULL
# 
# e1=performance(pred1,measure2)
# e2=performance(pred2,measure2)
# e3=performance(pred3,measure2)
# e4=performance(pred4,measure2)
# 
# lasso=e1%>%as.vector()
# enet=e2%>%as.vector()
# rf=e3%>%as.vector()
# rdg=e4%>%as.vector()
# 
# labs=c("BER","BRIER","FNR","FPR","MMCE","FDR")
# 
# scores=list("Lasso"=lasso,"Elasticnet"=enet,"RandomForest"=rf,"Ridge"=rdg)
# scores=list("Elasticnet"=enet,"RandomForest"=rf,"Ridge"=rdg)
# scores=list("Lasso"=lasso,"RandomForest"=rf,"Ridge"=rdg)
# chartJSRadar(scores = scores, labs = labs, maxScale =0.8)

# plot setting 
# library(ggplot2)
# library(ggrepel)
# 
# my_theme <- function(base_size = 12, base_family = "sans"){
#   theme_minimal(base_size = base_size, base_family = base_family) +
#     theme(
#       axis.text = element_text(size = 12),
#       axis.title = element_text(size = 14),
#       panel.grid.major = element_line(color = "grey"),
#       panel.grid.minor = element_blank(),
#       panel.background = element_rect(fill = "aliceblue"),
#       strip.background = element_rect(fill = "darkgrey", color = "grey", size = 1),
#       strip.text = element_text(face = "bold", size = 12, color = "white"),
#       legend.position = "right",
#       legend.justification = "top", 
#       panel.border = element_rect(color = "grey", fill = NA, size = 0.5)
#     )
# }
# 
# # performance metrics where smaller is better -> order with decreasing = FALSE
# sort_options_1 <- c("mean_per_class_error", "mse", "err")
# 
# for (sort_by_1 in sort_options_1) {
#   
#   grid <- h2o.getGrid("dl_grid", sort_by = sort_by_1, decreasing = FALSE)
#   
#   model_ids <- grid@model_ids
#   best_model <- h2o.getModel(model_ids[[1]])
#   
#   assign(paste0("best_model_", sort_by_1), best_model)
#   
# }
# 
# # performance metrics where bigger is better -> order with decreasing = TRUE
# sort_options_2 <- c("auc", "precision", "accuracy", "recall", "specificity")
# 
# for (sort_by_2 in sort_options_2) {
#   
#   grid <- h2o.getGrid("dl_grid", sort_by = sort_by_2, decreasing = TRUE)
#   
#   model_ids <- grid@model_ids
#   best_model <- h2o.getModel(model_ids[[1]])
#   
#   assign(paste0("best_model_", sort_by_2), best_model)
#   
# }
# 
# library(tibble)
# 
# sort_options <- c("mean_per_class_error", "mse", "err", "auc", "precision", "accuracy", "recall", "specificity")
# best_model = NULL
# errors_df_comb = NULL
# for (sort_by in sort_options) {
#   
#   best_model <- get(paste0("best_model_", sort_by))
#   errors <- h2o.mean_per_class_error(best_model, train = TRUE, valid = TRUE, xval = TRUE)
#   
#   errors_df <- data.frame(model_id = best_model@model_id, sort = sort_by, errors) %>%
#     rownames_to_column(var = "rowname")
#   
#   if (sort_by == "mean_per_class_error") {
#     
#     errors_df_comb <- errors_df
#     
#   } else {
#     
#     errors_df_comb <- rbind(errors_df_comb, errors_df)
#     
#   }
# }
# 
# order <- subset(errors_df_comb, rowname == "xval") %>%
#   arrange(errors)
# 
# errors_df_comb %>%
#   mutate(sort = factor(sort, levels = order$sort)) %>%
#   ggplot(aes(x = sort, y = errors, fill = model_id)) +
#   facet_grid(rowname ~ ., scales = "free") +
#   geom_bar(stat = "identity", alpha = 0.8) +
#   scale_fill_brewer(palette = "Set1") +
#   my_theme() +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
#         plot.margin = unit(c(0.5, 0, 0, 1), "cm")) +
#   labs(x = "")
# 
# 
# # model performance
# 
# for (sort_by in sort_options) {
#   
#   best_model <- get(paste0("best_model_", sort_by))
#   mse_auc_test <- data.frame(model_id = best_model@model_id,
#                              sort = sort_by, 
#                              mse = h2o.mse(h2o.performance(best_model, test_tc.hex)),
#                              auc = h2o.auc(h2o.performance(best_model, test_tc.hex)))
#   
#   if (sort_by == "mean_per_class_error") {
#     
#     mse_auc_test_comb <- mse_auc_test
#     
#   } else {
#     
#     mse_auc_test_comb <- rbind(mse_auc_test_comb, mse_auc_test)
#     
#   }
# }
# 
# library(tidyr)
# 
# mse_auc_test_comb %>%
#   gather(x, y, mse:auc) %>%
#   mutate(sort = factor(sort, levels = order$sort))%>%
#   ggplot(aes(x = sort, y = y, fill = model_id)) +
#   facet_grid(x ~ ., scales = "free") +
#   geom_bar(stat = "identity", alpha = 0.8, position = "dodge") +
#   scale_fill_brewer(palette = "Set1") +
#   my_theme() +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
#         plot.margin = unit(c(0.5, 0, 0, 1.5), "cm")) +
#   labs(x = "", y = "value", fill = "")
# 
# # 
# best_model = NULL
# finalRf_predictions = NULL
# finalRf_predictions_comb = NULL
# for (sort_by in sort_options) {
#   
#   best_model <- get(paste0("best_model_", sort_by))
#   
#   finalRf_predictions <- data.frame(model_id = rep(best_model@model_id, nrow(test_tc.hex)),
#                                     sort = rep(sort_by, nrow(test_tc.hex)),
#                                     
#                                     actual = as.vector(test_tc$vitalstatus), 
#                                     as.data.frame(h2o.predict(object = best_model, newdata = test_tc.hex)))
#   
#   
#   
#   finalRf_predictions$predict_stringent <- ifelse(finalRf_predictions$Dead > 0.8, "Dead", 
#                                                   ifelse(finalRf_predictions$Alive > 0.8, "Alive", "uncertain"))
#   finalRf_predictions$accurate_stringent <- ifelse(finalRf_predictions$actual == finalRf_predictions$predict_stringent, "yes", 
#                                                    ifelse(finalRf_predictions$predict_stringent == "uncertain", "na", "no"))
#   finalRf_predictions$cal_acc <- ifelse(finalRf_predictions$Alive > finalRf_predictions$Dead , "Alive", "Dead")
#   finalRf_predictions$accurate <- ifelse(finalRf_predictions$actual == finalRf_predictions$cal_acc, "yes", "no")
#   
#   
#   if (sort_by == "mean_per_class_error") {
#     
#     finalRf_predictions_comb <- finalRf_predictions
#     
#   } else {
#     
#     finalRf_predictions_comb <- rbind(finalRf_predictions_comb, finalRf_predictions)
#     
#   }
#   
#   tmp_dl_best = h2o.performance(best_model, test_tc.hex)
#   dl_best_data = h2o.confusionMatrix(tmp_dl_best)
#   print(dl_best_data)
#   print(best_model@model_id)
# }
# 
# 
# 
# finalRf_predictions_comb %>%
#   ggplot(aes(x = actual, fill = accurate)) +
#   geom_bar(position = "dodge") +
#   scale_fill_brewer(palette = "Set1") +
#   my_theme() +
#   facet_wrap(~ sort, ncol = 4) +
#   labs(fill = "Were\npredictions\naccurate?",
#        title = "Default predictions")
# 
# finalRf_predictions_comb %>%
#   subset(accurate_stringent != "na") %>%
#   ggplot(aes(x = actual, fill = accurate_stringent)) +
#   geom_bar(position = "dodge") +
#   scale_fill_brewer(palette = "Set1") +
#   my_theme() +
#   facet_wrap(~ sort, ncol = 4) +
#   labs(fill = "Were\npredictions\naccurate?",
#        title = "Stringent predictions")
# 
# 
# result_sum = data.frame()
# 
# for (model_id in unique(finalRf_predictions_comb$model_id)) {
#   tmp.finaldl = finalRf_predictions_comb[which(finalRf_predictions_comb$model_id == model_id),]
#   for (sort_id in unique(tmp.finaldl$sort)) {
#     tmp2.finaldl = tmp.finaldl[which(tmp.finaldl$sort == sort_id),]
#     acc_tmp = sum(tmp2.finaldl$accurate == "yes") / nrow(test_tc) * 100
#     print(acc_tmp)
#     tmp_df = data.frame(model_id = model_id, sort = sort_id, acc = acc_tmp)
#     
#     result_sum = rbind(result_sum,tmp_df)
#     tmp_df = NULL
#     
#   } 
#   
#   acc_tmp = NULL
#   tmp2.finaldl =NULL
#   tmp.finaldl = NULL
#   
# } 
# 
# result_sum

# set.seed(12345)
# rfseed=sample(10000,100,replace=F)
# 
# res=NULL
# pdf=NULL
# train_tc
# for (i in 1:100){
#   
#   rfmod=h2o.randomForest(x = features,
#                          y = "vitalstatus",
#                          training_frame = train_tc.hex,
#                          validation_frame = valid_tc.hex,
#                          nfolds=10,
#                          fold_assignment = "Stratified",
#                          balance_classes = TRUE,
#                          ntrees = 100, max_depth = 50,mtries = -1,sample_rate = 0.7,
#                          stopping_metric = "misclassification",
#                          stopping_tolerance = 0.01,
#                          stopping_rounds = 3,
#                          keep_cross_validation_fold_assignment = F,
#                          keep_cross_validation_predictions=F,
#                          score_each_iteration = TRUE,
#                          seed=rfseed[i])
#   
#   pdftemp=predict(rfmod,test_tc.hex)%>%as_tibble()%>%mutate(.,Id=row.names(test_tc),Truth=test_tc$vitalstatus,Accuracy=ifelse(test_tc$vitalstatus ==.$predict, "Correct", "Wrong"),Model=i)
#   
#   vimp=rfmod@model$variable_importances%>%as_tibble()%>%mutate(.,Seed=rfseed[i],Model=i)
#   
#   res=rbind(res,vimp)
#   
#   pdf=rbind(pdf,pdftemp)
#   
# }
# 
# res=as_tibble(res)
# res$Model=as.factor(res$Model)
# pdf$Id=as.factor(pdf$Id)
# 
# pdf%>%ggplot(aes(x=Model,y=reorder(Id,Dead),fill=Dead))+geom_tile()+theme_bw()+scale_fill_gradient2(low="green",mid="gold",high="red",midpoint = 0.5)+facet_wrap(~Truth,shrink=TRUE,scale="free",ncol=1)+theme_bw(8)+theme(axis.text.x = element_text(angle =45, hjust = 1))
# 
# pdf%>%ggplot(aes(x=Model,y=reorder(Id,1-Dead),fill=Accuracy))+geom_tile()+theme_bw()+scale_fill_manual(values=c("skyblue","red"))+theme_bw(8)+theme(axis.text.x = element_text(angle =45, hjust = 1))+facet_wrap(~Truth,shrink=TRUE,scale="free",ncol=1)+scale_y_discrete("Idx")
# 
# # average predicted probability and averaged accuracy of 100 models
# pdf2=NULL
# getwd()
# library(Hmisc)
# 
# 
# write.csv(pdf, paste0(main.path_tc,"randomforest_1st_results3.csv"))
# for (i in 1:length(rownames(test_tc))){
#   Idt=pdf$Id[i]
#   dtemp=subset(pdf,Id==Idt)
#   sum=Hmisc::describe(dtemp)
#   
#   predprob=sum$Dead$counts[[5]]
#   acc=sum$Accuracy$values[[2]][1]/100
#   wrong=sum$Accuracy$values[[2]][2]/100
#   
#   pdftemp=cbind.data.frame(Id=Idt,predprob,acc,wrong)
#   pdf2=rbind(pdf2,pdftemp)
# }
# 
# pdf2$class=test_tc$vitalstatus
# pdf2$predprob=as.numeric(as.character(pdf2$predprob))
# 
# pdf2$Id=as.factor(pdf2$Id)
# 
# pdf2$wrong[is.na(pdf2$wrong)] <- 0
# 
# pdf2%>%ggplot(aes(x=reorder(Id,predprob),y=predprob,fill=predprob))+geom_bar(alpha=0.8,stat="identity")+theme_bw(8)+scale_fill_gradient2(low="gold",mid="red",high="purple",midpoint = 0.4)+facet_wrap(~class,shrink = TRUE,ncol=1,scales = "free")+geom_hline(yintercept = 0.5,linetype=2,size=1)
# 
# pdf2%>%gather(wrong,acc,key="Accuracy",value="Prediction")%>%ggplot(aes(x=reorder(Id,-predprob),y=Prediction,fill=Accuracy))+geom_bar(alpha=0.7,position="fill",stat="Identity")+theme_bw(8)+facet_wrap(~class,shrink = TRUE,ncol=1,scales = "free")+scale_fill_manual(values=c("skyblue","red"))+scale_x_discrete("Id")
# 
# # relative importance of variables (more important a feature is, more saturated in color that feature would be)
# 
# res%>%ggplot(aes(x=reorder(Model,relative_importance),y=reorder(variable,-relative_importance),fill=relative_importance,color=relative_importance))+geom_tile(show.legend = F)+theme_bw()+scale_fill_gradient2(low="white",mid="red",high="purple",midpoint = 50)+scale_color_gradient2(low="white",mid="red",high="purple",midpoint = 50)+theme_bw(5)+theme(axis.text.x = element_text(angle =45, hjust = 1))+scale_x_discrete("Models")+scale_y_discrete("Features")
# 
# 
# res%>%ggplot(aes(x=reorder(Model,scaled_importance),y=reorder(variable,-scaled_importance),fill=scaled_importance,color=scaled_importance))+geom_tile(show.legend = F)+theme_bw()+scale_fill_gradient2(low="white",mid="gold",high="red",midpoint = 0.5)+scale_color_gradient2(low="white",mid="gold",high="red",midpoint =0.5)+theme_bw(5)+theme(axis.text.x = element_text(angle =45, hjust = 1))+scale_x_discrete("Models")+scale_y_discrete("Features")
# 
# # extract averaged importance score
# 
# res2=NULL
# 
# for (i in 1:length(colnames(test_tc[,-length(test_tc)]))){
#   dtemp=subset(res,variable==features[i])
#   sum=psych::describe(dtemp)
#   rei=sum[[5]][2]
#   sci=sum[[5]][3]
#   pct=sum[[5]][4]
#   
#   rtemp=cbind.data.frame(rei,sci,pct,var=dtemp[[1]][1])
#   res2=rbind(res2,rtemp)
#   
# }
# 
# res2%>%ggplot(aes(x=reorder(var,rei),y=rei,fill=rei,color=rei))+geom_segment(aes(x=reorder(var,rei),xend=var,y=0,yend=rei),size=1,show.legend = F)+geom_point(size=1,show.legend = F)+scale_x_discrete("Features")+scale_y_continuous("Relative importance")+coord_flip()+scale_fill_gradient(low="blue",high="red")+scale_color_gradient(low="blue",high="red")+ggtitle("Relative importance")+theme_bw(5)
# 
# # extract top150 variables
# res3=res2%>%arrange(desc(rei))%>%.[c(1:150),]
# 
# res3%>%ggplot(aes(x=reorder(var,rei),y=rei,fill=rei,color=rei))+geom_segment(aes(x=reorder(var,rei),xend=var,y=0,yend=rei),size=1,show.legend = F)+geom_point(size=1,show.legend = F)+scale_x_discrete("Features")+scale_y_continuous("Relative importance")+coord_flip()+scale_fill_gradient(low="blue",high="red")+scale_color_gradient(low="blue",high="red")+ggtitle("Relative importance")+theme_bw(5)
# 
# # re-train by top 150 variables
# features2=as.vector(res3$var)
# 
# rfmod=h2o.randomForest(x = features2,
#                        y = "vitalstatus",
#                        training_frame = train_tc.hex,
#                        validation_frame = valid_tc.hex,
#                        nfolds=10,
#                        fold_assignment = "Stratified",
#                        balance_classes = TRUE,
#                        ntrees = 100, max_depth = 50,mtries = -1,sample_rate = 0.8,
#                        stopping_metric = "misclassification",
#                        stopping_tolerance = 0.01,
#                        stopping_rounds = 3,
#                        keep_cross_validation_fold_assignment = F,
#                        keep_cross_validation_predictions=F,
#                        score_each_iteration = TRUE,
#                        seed=333)


# # Get the validation accuracy (AUC)
# 
# xval_tc_performance <- h2o.performance(train_test_tc_dnn_each,xval = T)
# xval_tc_performance@metrics$AUC
# 
# # Get the training accuracy (AUC)
# 
# train_tc_performance <- h2o.performance(train_test_tc_dnn_each,train = T)
# train_tc_performance@metrics$AUC
# 
# # Get the testing accuracy(AUC)
# 
# test_tc_performance <- h2o.performance(train_test_tc_dnn_each,valid = T)
# test_tc_performance@metrics$AUC
# 
# # now make a prediction
# 
# predictions <- h2o.predict(train_test_tc_dnn_each, test_tc.hex)
# df.predictions = as.data.frame(predictions)
# predictions_dnn = table(prediction = df.predictions$predict, actual = test_tc$vitalstatus)
# acc_predicted_tc_dnn = round((predictions_dnn[1,1] + predictions_dnn[2,2]) / sum(predictions_dnn) * 100, 2)
# acc_predicted_tc_dnn
# predictions
# ?h2o.predict
# result<-as.data.frame(round(predictions))
# result<-as.data.frame(predictions)
# 
# 
# h2o.accuracy(test_tc_performance )
# h2o.specificity(test_tc_performance)
# h2o.sensitivity(test_tc_performance)
# h2o.confusionMatrix(test_tc_performance, thresholds = 0.962664373663552)
# 
# # Automl
# ?h2o.automl
# aml_tc <- h2o.automl(x = x, y = y,
#                      training_frame = train_tc.hex,
#                      validation_frame = vaild_tc.hex,
#                      max_models = 20,
#                      seed = 1)
# 
# # View the AutoML Leaderboard
# 
# lb_tc <- aml_tc@leaderboard
# print(lb_tc, n = nrow(lb_tc))  # Print all rows instead of default (6 rows)
# 
# # The leader model is stored here
# aml_tc@leader
# 
# # If you need to generate predictions on a test set, you can make
# # predictions directly on the `"H2OAutoML"` object, or on the leader
# # model object directly
# 
# pred_tc <- h2o.predict(aml_tc, test_tc.hex)  # predict(aml, test) also works
# df.pred = as.data.frame(pred_tc)
# 
# predicted_tc = table(prediction = df.pred$predict, actual = test_tc$vitalstatus)
# acc_predicted_tc = round((predicted_tc[1,1] + predicted_tc[2,2]) / sum(predicted_tc) * 100, 2)
# 
# 
# # save the model
# model_path <- h2o.saveModel(object=model, path=getwd(), force=TRUE)
# 
# # load the model
# saved_model <- h2o.loadModel(model_path)
