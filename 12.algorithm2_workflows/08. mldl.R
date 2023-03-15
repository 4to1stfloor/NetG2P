# 1) call input

filepath = "/home/seokwon/"
ref_path = paste0(filepath, "99.reference/")
kindofCancer = "brca"

Cancername = "TCGA-BRCA" # you can change this
main.path_tc = paste0(filepath, "00.data/","30.", Cancername,"/")

ex_path = paste0(filepath,"01.externaldata/", kindofCancer, "/") 
CancerType_in = "brca_metabric" # you can change this
main.path_in = paste0(ex_path,CancerType_in, "/")

CancerType = "TCGA-BRCA"

# plot setting 
library(ggplot2)
library(ggrepel)

my_theme <- function(base_size = 12, base_family = "sans"){
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      panel.grid.major = element_line(color = "grey"),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "aliceblue"),
      strip.background = element_rect(fill = "darkgrey", color = "grey", size = 1),
      strip.text = element_text(face = "bold", size = 12, color = "white"),
      legend.position = "right",
      legend.justification = "top", 
      panel.border = element_rect(color = "grey", fill = NA, size = 0.5)
    )
}

# A - TCGA

data_tc = readRDS(file = paste0(main.path_tc,CancerType,"_pathwaylink_all.rds"))
data_tc = readRDS(file = paste0(main.path_tc,"/",CancerType,"_pathwaylink_all_0or1.rds"))
data_tc = readRDS(file = paste0(main.path_tc,"/",CancerType,"_pathwaylink_all_1or1000000.rds"))
data_tc = readRDS(file = paste0(main.path_tc,"/",CancerType,"_pathwaylink_all_0_000001or1000000.rds"))

# data_tc = data_tc[,-which(dim(data_tc)[1] == colSums(data_tc[,-length(colnames(data_tc))]))]
data_tc = data_tc[,-which(colSums(data_tc[,-length(colnames(data_tc))]) == 0)]
# B - exteral

data_ex = readRDS(paste0(main.path_in, CancerType_in,"_pathwaylink_all_0or1.rds"))
data_ex = readRDS(paste0(main.path_in, CancerType_in,"_pathwaylink_all_1or1000000.rds"))
data_ex = readRDS(paste0(main.path_in,CancerType_in,"_pathwaylink_all_0_000001or1000000.rds"))

# 2) build a model
library(h2o)
localH2O = h2o.init(ip = "localhost", port = 54321, startH2O = TRUE,min_mem_size = "10G",nthreads = 24)
set.seed(1234)

train_val = sample(1:nrow(data_tc), nrow(data_tc)*0.6) 
testtmp_val = -train_val

# (1) TCGA data

train_tc <- data_tc[train_val,]
testtmp_tc <- data_tc[testtmp_val,]
valid_val = sample(1:nrow(testtmp_tc), nrow(testtmp_tc)*0.5) 

test_val = -valid_val
valid_tc <- testtmp_tc[valid_val,]
test_tc <-  testtmp_tc[test_val,]

# (1) Define input (features) and output (response) variables"

features = setdiff(colnames(data_tc),"vitalstatus")
response = colnames(data_tc)[length(colnames(data_tc))]

# (2) exteranldata (test_all_data sample)

features_ex = setdiff(colnames(data_ex),"vitalstatus")
response_ex = colnames(data_ex)[length(colnames(data_ex))]

# (1) - Convert the outcome variable into factor

train_tc$vitalstatus <- as.factor(train_tc$vitalstatus)
valid_tc$vitalstatus <- as.factor(valid_tc$vitalstatus)
test_tc$vitalstatus <- as.factor(test_tc$vitalstatus)

# (2) - Convert the outcome variable into factor

data_ex$vitalstatus <- as.factor(data_ex$vitalstatus)

# (1) - Convert the data into h2o form (TCGA)

train_tc.hex <- as.h2o(x = train_tc, destination_frame = "train_tc.hex")
test_tc.hex <- as.h2o(x = test_tc, destination_frame = "test_tc.hex")
valid_tc.hex <- as.h2o(x = valid_tc, destination_frame = "valid_tc.hex")

# (2) - Convert the data into h2o form (external)

data_ex.hex <- as.h2o(x = data_ex, destination_frame = "data_ex.hex")

##############################################################################
###### 1. model building - Random search grid
##############################################################################

# set hyper-parameter

hyper_params <- list(
  activation = c("Rectifier", "Maxout", "Tanh", "RectifierWithDropout", "MaxoutWithDropout", "TanhWithDropout"), 
  hidden = list(c(5, 5, 5), c(10, 10, 10), c(50, 50, 50), c(100, 100, 100)),
  epochs = c(50, 100, 200),
  l1 = c(0, 0.00001, 0.0001), 
  l2 = c(0, 0.00001, 0.0001),
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

# search Random grid 
remove(dl_grid)
dl_grid <- h2o.grid(algorithm = "deeplearning", 
                    x = features,
                    y = response,
                    grid_id = "dl_grid",
                    training_frame = train_tc.hex,
                    validation_frame = valid_tc.hex,
                    nfolds = 5,                           
                    fold_assignment = "Stratified",
                    hyper_params = hyper_params,
                    search_criteria = search_criteria,
                    seed = 42
)
?h2o.grid()

hyper_params2 <- list(ntrees = seq(50, 500, by = 50),
                      mtries = seq(3, 5, by = 1),
                      max_depth = seq(10, 30, by = 10),
                      min_rows = seq(1, 3, by = 1),
                      nbins = seq(20, 30, by = 10),
                      sample_rate = c(0.55, 0.632, 0.75))

d2_grid <- h2o.grid(algorithm = "randomForest", 
                    x = features,
                    y = response,
                    grid_id = "d2_grid",
                    training_frame = train_tc.hex,
                    validation_frame = valid_tc.hex,
                    nfolds = 5,                           
                    fold_assignment = "Stratified",
                    hyper_params = hyper_params2,
                    search_criteria = search_criteria,
                    seed = 42
)

# performance metrics where smaller is better -> order with decreasing = FALSE
sort_options_1 <- c("mean_per_class_error", "mse", "err")

for (sort_by_1 in sort_options_1) {
  
  grid <- h2o.getGrid("dl_grid", sort_by = sort_by_1, decreasing = FALSE)
  
  model_ids <- grid@model_ids
  best_model <- h2o.getModel(model_ids[[1]])
  
  assign(paste0("best_model_", sort_by_1), best_model)
  
}

# performance metrics where bigger is better -> order with decreasing = TRUE
sort_options_2 <- c("auc", "precision", "accuracy", "recall", "specificity")

for (sort_by_2 in sort_options_2) {
  
  grid <- h2o.getGrid("dl_grid", sort_by = sort_by_2, decreasing = TRUE)
  
  model_ids <- grid@model_ids
  best_model <- h2o.getModel(model_ids[[1]])
  
  assign(paste0("best_model_", sort_by_2), best_model)
  
}

library(tibble)

sort_options <- c("mean_per_class_error", "mse", "err", "auc", "precision", "accuracy", "recall", "specificity")
best_model = NULL
errors_df_comb = NULL
for (sort_by in sort_options) {
  
  best_model <- get(paste0("best_model_", sort_by))
  errors <- h2o.mean_per_class_error(best_model, train = TRUE, valid = TRUE, xval = TRUE)
  
  errors_df <- data.frame(model_id = best_model@model_id, sort = sort_by, errors) %>%
    rownames_to_column(var = "rowname")
  
  if (sort_by == "mean_per_class_error") {
    
    errors_df_comb <- errors_df
    
  } else {
    
    errors_df_comb <- rbind(errors_df_comb, errors_df)
    
  }
}

order <- subset(errors_df_comb, rowname == "xval") %>%
  arrange(errors)

errors_df_comb %>%
  mutate(sort = factor(sort, levels = order$sort)) %>%
  ggplot(aes(x = sort, y = errors, fill = model_id)) +
  facet_grid(rowname ~ ., scales = "free") +
  geom_bar(stat = "identity", alpha = 0.8) +
  scale_fill_brewer(palette = "Set1") +
  my_theme() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.margin = unit(c(0.5, 0, 0, 1), "cm")) +
  labs(x = "")


# model performance

for (sort_by in sort_options) {
  
  best_model <- get(paste0("best_model_", sort_by))
  mse_auc_test <- data.frame(model_id = best_model@model_id,
                             sort = sort_by, 
                             mse = h2o.mse(h2o.performance(best_model, test_tc.hex)),
                             auc = h2o.auc(h2o.performance(best_model, test_tc.hex)))
  
  if (sort_by == "mean_per_class_error") {
    
    mse_auc_test_comb <- mse_auc_test
    
  } else {
    
    mse_auc_test_comb <- rbind(mse_auc_test_comb, mse_auc_test)
    
  }
}

library(tidyr)

mse_auc_test_comb %>%
  tidyr::gather(x, y, mse:auc) %>%
  mutate(sort = factor(sort, levels = order$sort))%>%
  ggplot(aes(x = sort, y = y, fill = model_id)) +
  facet_grid(x ~ ., scales = "free") +
  geom_bar(stat = "identity", alpha = 0.8, position = "dodge") +
  scale_fill_brewer(palette = "Set1") +
  my_theme() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.margin = unit(c(0.5, 0, 0, 1.5), "cm")) +
  labs(x = "", y = "value", fill = "")

# 
best_model = NULL
finalRf_predictions = NULL
finalRf_predictions_comb = NULL
for (sort_by in sort_options) {
  
  best_model <- get(paste0("best_model_", sort_by))
  
  finalRf_predictions <- data.frame(model_id = rep(best_model@model_id, nrow(test_tc.hex)),
                                    sort = rep(sort_by, nrow(test_tc.hex)),
                                    
                                    actual = as.vector(test_tc$vitalstatus), 
                                    as.data.frame(h2o.predict(object = best_model, newdata = test_tc.hex)))
  
  
  
  finalRf_predictions$predict_stringent <- ifelse(finalRf_predictions$Dead > 0.8, "Dead", 
                                                  ifelse(finalRf_predictions$Alive > 0.8, "Alive", "uncertain"))
  finalRf_predictions$accurate_stringent <- ifelse(finalRf_predictions$actual == finalRf_predictions$predict_stringent, "yes", 
                                                   ifelse(finalRf_predictions$predict_stringent == "uncertain", "na", "no"))
  finalRf_predictions$cal_acc <- ifelse(finalRf_predictions$Alive > finalRf_predictions$Dead , "Alive", "Dead")
  finalRf_predictions$accurate <- ifelse(finalRf_predictions$actual == finalRf_predictions$cal_acc, "yes", "no")
  
  
  if (sort_by == "mean_per_class_error") {
    
    finalRf_predictions_comb <- finalRf_predictions
    
  } else {
    
    finalRf_predictions_comb <- rbind(finalRf_predictions_comb, finalRf_predictions)
    
  }
  
  tmp_dl_best = h2o.performance(best_model, test_tc.hex)
  dl_best_data = h2o.confusionMatrix(tmp_dl_best)
  print(dl_best_data)
  print(best_model@model_id)
}

finalRf_predictions_comb %>%
  ggplot(aes(x = actual, fill = accurate)) +
  geom_bar(position = "dodge") +
  scale_fill_brewer(palette = "Set1") +
  my_theme() +
  facet_wrap(~ sort, ncol = 4) +
  labs(fill = "Were\npredictions\naccurate?",
       title = "Default predictions")

finalRf_predictions_comb %>%
  subset(accurate_stringent != "na") %>%
  ggplot(aes(x = actual, fill = accurate_stringent)) +
  geom_bar(position = "dodge") +
  scale_fill_brewer(palette = "Set1") +
  my_theme() +
  facet_wrap(~ sort, ncol = 4) +
  labs(fill = "Were\npredictions\naccurate?",
       title = "Stringent predictions")


result_sum = data.frame()

for (model_id in unique(finalRf_predictions_comb$model_id)) {
  tmp.finaldl = finalRf_predictions_comb[which(finalRf_predictions_comb$model_id == model_id),]
  for (sort_id in unique(tmp.finaldl$sort)) {
    tmp2.finaldl = tmp.finaldl[which(tmp.finaldl$sort == sort_id),]
    acc_tmp = sum(tmp2.finaldl$accurate == "yes") / nrow(test_tc) * 100
    print(acc_tmp)
    tmp_df = data.frame(model_id = model_id, sort = sort_id, acc = acc_tmp)
    
    result_sum = rbind(result_sum,tmp_df)
    tmp_df = NULL
    
  } 
  
  acc_tmp = NULL
  tmp2.finaldl =NULL
  tmp.finaldl = NULL
  
} 
sum(finalRf_predictions_comb$actual == "Alive") / length(finalRf_predictions_comb$actual)
result_sum
test_model = get(paste0("best_model_accuracy"))
h2o.make_leaderboard(dl_grid, test_tc.hex)
h2o.performance(test_model,test_tc.hex)
## randomforest 
sort_options_1 <- c("mean_per_class_error", "mse", "err")

for (sort_by_1 in sort_options_1) {
  
  grid <- h2o.getGrid("d2_grid", sort_by = sort_by_1, decreasing = FALSE)
  
  model_ids <- grid@model_ids
  best_model <- h2o.getModel(model_ids[[1]])
  
  assign(paste0("best_model_", sort_by_1), best_model)
  
}

# performance metrics where bigger is better -> order with decreasing = TRUE
sort_options_2 <- c("auc", "precision", "accuracy", "recall", "specificity")

for (sort_by_2 in sort_options_2) {
  
  grid <- h2o.getGrid("d2_grid", sort_by = sort_by_2, decreasing = TRUE)
  
  model_ids <- grid@model_ids
  best_model <- h2o.getModel(model_ids[[1]])
  
  assign(paste0("best_model_", sort_by_2), best_model)
  
}

library(tibble)

sort_options <- c("mean_per_class_error", "mse", "err", "auc", "precision", "accuracy", "recall", "specificity")
best_model = NULL
errors_df_comb = NULL
for (sort_by in sort_options) {
  
  best_model <- get(paste0("best_model_", sort_by))
  errors <- h2o.mean_per_class_error(best_model, train = TRUE, valid = TRUE, xval = TRUE)
  
  errors_df <- data.frame(model_id = best_model@model_id, sort = sort_by, errors) %>%
    rownames_to_column(var = "rowname")
  
  if (sort_by == "mean_per_class_error") {
    
    errors_df_comb <- errors_df
    
  } else {
    
    errors_df_comb <- rbind(errors_df_comb, errors_df)
    
  }
}

order <- subset(errors_df_comb, rowname == "xval") %>%
  arrange(errors)

errors_df_comb %>%
  mutate(sort = factor(sort, levels = order$sort)) %>%
  ggplot(aes(x = sort, y = errors, fill = model_id)) +
  facet_grid(rowname ~ ., scales = "free") +
  geom_bar(stat = "identity", alpha = 0.8) +
  scale_fill_brewer(palette = "Set1") +
  my_theme() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.margin = unit(c(0.5, 0, 0, 1), "cm")) +
  labs(x = "")


# model performance

for (sort_by in sort_options) {
  
  best_model <- get(paste0("best_model_", sort_by))
  mse_auc_test <- data.frame(model_id = best_model@model_id,
                             sort = sort_by, 
                             mse = h2o.mse(h2o.performance(best_model, test_tc.hex)),
                             auc = h2o.auc(h2o.performance(best_model, test_tc.hex)))
  
  if (sort_by == "mean_per_class_error") {
    
    mse_auc_test_comb <- mse_auc_test
    
  } else {
    
    mse_auc_test_comb <- rbind(mse_auc_test_comb, mse_auc_test)
    
  }
}

library(tidyr)

mse_auc_test_comb %>%
  tidyr::gather(x, y, mse:auc) %>%
  mutate(sort = factor(sort, levels = order$sort))%>%
  ggplot(aes(x = sort, y = y, fill = model_id)) +
  facet_grid(x ~ ., scales = "free") +
  geom_bar(stat = "identity", alpha = 0.8, position = "dodge") +
  scale_fill_brewer(palette = "Set1") +
  my_theme() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.margin = unit(c(0.5, 0, 0, 1.5), "cm")) +
  labs(x = "", y = "value", fill = "")

# 
best_model = NULL
finalRf_predictions = NULL
finalRf_predictions_comb = NULL
for (sort_by in sort_options) {
  
  best_model <- get(paste0("best_model_", sort_by))
  
  finalRf_predictions <- data.frame(model_id = rep(best_model@model_id, nrow(test_tc.hex)),
                                    sort = rep(sort_by, nrow(test_tc.hex)),
                                    
                                    actual = as.vector(test_tc$vitalstatus), 
                                    as.data.frame(h2o.predict(object = best_model, newdata = test_tc.hex)))
  
  
  
  finalRf_predictions$predict_stringent <- ifelse(finalRf_predictions$Dead > 0.8, "Dead", 
                                                  ifelse(finalRf_predictions$Alive > 0.8, "Alive", "uncertain"))
  finalRf_predictions$accurate_stringent <- ifelse(finalRf_predictions$actual == finalRf_predictions$predict_stringent, "yes", 
                                                   ifelse(finalRf_predictions$predict_stringent == "uncertain", "na", "no"))
  finalRf_predictions$cal_acc <- ifelse(finalRf_predictions$Alive > finalRf_predictions$Dead , "Alive", "Dead")
  finalRf_predictions$accurate <- ifelse(finalRf_predictions$actual == finalRf_predictions$cal_acc, "yes", "no")
  
  
  if (sort_by == "mean_per_class_error") {
    
    finalRf_predictions_comb <- finalRf_predictions
    
  } else {
    
    finalRf_predictions_comb <- rbind(finalRf_predictions_comb, finalRf_predictions)
    
  }
  
  tmp_dl_best = h2o.performance(best_model, test_tc.hex)
  dl_best_data = h2o.confusionMatrix(tmp_dl_best)
  print(dl_best_data)
  print(best_model@model_id)
}

finalRf_predictions_comb %>%
  ggplot(aes(x = actual, fill = accurate)) +
  geom_bar(position = "dodge") +
  scale_fill_brewer(palette = "Set1") +
  my_theme() +
  facet_wrap(~ sort, ncol = 4) +
  labs(fill = "Were\npredictions\naccurate?",
       title = "Default predictions")

finalRf_predictions_comb %>%
  subset(accurate_stringent != "na") %>%
  ggplot(aes(x = actual, fill = accurate_stringent)) +
  geom_bar(position = "dodge") +
  scale_fill_brewer(palette = "Set1") +
  my_theme() +
  facet_wrap(~ sort, ncol = 4) +
  labs(fill = "Were\npredictions\naccurate?",
       title = "Stringent predictions")


result_sum = data.frame()

for (model_id in unique(finalRf_predictions_comb$model_id)) {
  tmp.finaldl = finalRf_predictions_comb[which(finalRf_predictions_comb$model_id == model_id),]
  for (sort_id in unique(tmp.finaldl$sort)) {
    tmp2.finaldl = tmp.finaldl[which(tmp.finaldl$sort == sort_id),]
    acc_tmp = sum(tmp2.finaldl$accurate == "yes") / nrow(test_tc) * 100
    print(acc_tmp)
    tmp_df = data.frame(model_id = model_id, sort = sort_id, acc = acc_tmp)
    
    result_sum = rbind(result_sum,tmp_df)
    tmp_df = NULL
    
  } 
  
  acc_tmp = NULL
  tmp2.finaldl =NULL
  tmp.finaldl = NULL
  
} 

result_sum

sum(finalRf_predictions_comb$actual == "Alive") / length(finalRf_predictions_comb$actual)
result_sum
test_model = get(paste0("best_model_accuracy"))
h2o.make_leaderboard(d2_grid, test_tc.hex)
h2o.make_leaderboard(dl_grid, test_tc.hex)
h2o.make_leaderboard(d2_grid, data_ex.hex)
h2o.make_leaderboard(dl_grid, data_ex.hex)

test_model = h2o.getModel("d2_grid_model_13")
h2o.performance(test_model,test_tc.hex)
h2o.performance(test_model,data_ex.hex)

# Use same data as above

# GBM hyperparamters (bigger grid than above)
gbm_params2 <- list(learn_rate = seq(0.01, 1, 0.01),
                    max_depth = seq(2, 30, 1),
                    sample_rate = seq(0, 2, 0.1),
                    col_sample_rate = seq(0.1, 1.0, 0.1))

gbm_params2 <- list(learn_rate = seq(0.01, 1, 0.01),
                    max_depth = seq(2, 30, 1))

search_criteria2 <- list(strategy = "RandomDiscrete", max_models = 36, seed = 1)

# Train and validate a random grid of GBMs
?h2o.grid()
gbm_grid2 <- h2o.grid("gbm", x = features, y = response,
                      grid_id = "gbm_grid2",
                      training_frame = train_tc.hex,
                      validation_frame = valid_tc.hex,
                      ntrees = 100,
                      seed = 12,
                      nfolds=5,
                      keep_cross_validation_predictions=TRUE,
                      hyper_params = gbm_params2,
                      search_criteria = search_criteria2)

gbm_gridperf2 <- h2o.getGrid(grid_id = "gbm_grid2",
                             sort_by = "auc",
                             decreasing = TRUE)
gbm_gridperf2 <- h2o.getGrid(grid_id = "gbm_grid2",
                             sort_by = "mse",
                             decreasing = FALSE)
print(gbm_gridperf2)
?h2o.getGrid

# Grab the top GBM model, chosen by validation AUC
best_gbm2 <- h2o.getModel(gbm_gridperf2@model_ids[[1]])

# Now let's evaluate the model performance on a test set
# so we get an honest estimate of top model performance
print(gbm_gridperf2)
View(gbm_gridperf2)
# Grab the top GBM model, chosen by validation AUC
best_gbm1 <- h2o.getModel(gbm_gridperf2@model_ids[[1]])
h2o.make_leaderboard(gbm_grid2, test_tc.hex)
h2o.make_leaderboard(gbm_grid2, data_ex.hex)
test.model = h2o.getModel("gbm_grid2_model_9")
h2o.performance(model = test.model, newdata = test_tc.hex)
h2o.performance(model = test.model, newdata = data_ex.hex)
# Now let's evaluate the model performance on a test set
# so we get an honest estimate of top model performance
?h2o.performance
best_xgb_perf2 <- h2o.performance(model = best_xgb1,
                                  newdata = test_tc.hex)
h2o.varimp_plot(test.model)

plot(h2o.performance(model = test.model, newdata = test_tc.hex) , type = "roc")
plot(h2o.performance(model = test.model, newdata = data_ex.hex) , type = "roc")

h2o.gains_lift_plot(test.model)
h2o.confusionMatrix(test.model)
h2o.confusionMatrix(object = test.model,test_tc.hex)
h2o.confusionMatrix(object = test.model,data_ex.hex)
h2o.auc(h2o.performance(test.model,data_ex.hex))
df_merge_confusion = rbind(h2o.confusionMatrix(test.model),h2o.confusionMatrix(object = test.model,test_tc.hex),h2o.confusionMatrix(object = test.model,data_ex.hex))
write.csv(df_merge_confusion,paste0(main.path_in, "/gbm_randomgrid",h2o.auc(h2o.performance(test.model,data_ex.hex)) ,"_confusionmatrix.csv"))
h2o.saveModel(test.model, path = paste0(main.path_in,"/gbm_randomgrid_0.664"))

# xbm

# Train & Cross-validate a (shallow) XGB-GBM
?h2o.xgboost

hyper_xgb_params <- list(ntrees = seq(10, 1000, 1),
                     learn_rate = seq(0.0001, 0.2, 0.0001),
                     max_depth = seq(1, 20, 1),
                     sample_rate = seq(0.5, 1.0, 0.0001),
                     col_sample_rate = seq(0.2, 1.0, 0.0001))

search_criteria2 <- list(strategy = "RandomDiscrete", max_models = 36, seed = 1)

xgb_grid1 <- h2o.grid("xgboost", x = features, y = response,
                      grid_id = "xgb_grid1",
                      training_frame = train_tc.hex,
                      validation_frame = valid_tc.hex,
                      seed = 1,
                      nfolds=5,
                      keep_cross_validation_predictions=TRUE,
                      hyper_params = hyper_xgb_params,
                      search_criteria = search_criteria2)

xgb_gridperf1 <- h2o.getGrid(grid_id = "xgb_grid1",
                             sort_by = "mse",
                             decreasing = FALSE)

xgb_gridperf1 <- h2o.getGrid(grid_id = "xgb_grid1",
                             sort_by = "auc",
                             decreasing = TRUE)

# Grab the top GBM model, chosen by validation AUC
best_xgb1 <- h2o.getModel(xgb_gridperf1@model_ids[[1]])
h2o.make_leaderboard(xgb_grid1, test_tc.hex)
h2o.make_leaderboard(xgb_grid1, data_ex.hex)
test.model = h2o.getModel("xgb_grid1_model_13")
h2o.performance(model = test.model, newdata = test_tc.hex)
h2o.performance(model = test.model, newdata = data_ex.hex)
# Now let's evaluate the model performance on a test set
# so we get an honest estimate of top model performance
?h2o.performance
best_xgb_perf2 <- h2o.performance(model = best_xgb1,
                                  newdata = test_tc.hex)
h2o.varimp_plot(test.model)

plot(h2o.performance(model = test.model, newdata = test_tc.hex) , type = "roc")
plot(h2o.performance(model = test.model, newdata = data_ex.hex) , type = "roc")

h2o.gains_lift_plot(test.model)



# # Train & Cross-validate another (deeper) XGB-GBM
# my_xgb2 <- h2o.xgboost(x = features,
#                        y = response,
#                        training_frame = train_tc.hex,
#                        validation_frame = valid_tc.hex,
#                        distribution = "bernoulli",
#                        ntrees = 50,
#                        max_depth = 9,
#                        min_rows = 1,
#                        learn_rate = 0.1,
#                        sample_rate = 0.7,
#                        col_sample_rate = 0.8,
#                        nfolds = 5,
#                        fold_assignment = "AUTO",
#                        keep_cross_validation_predictions = TRUE,
#                        reg_lambda = 0.7,
#                        seed = 1,
#                        )
# 


# Train a stacked ensemble using the H2O and XGBoost models from above
xgb_gridperf1 <- h2o.getGrid(grid_id = "xgb_grid1",
                             sort_by = "auc",
                             decreasing = TRUE)
d1_gridperf1 <- h2o.getGrid(grid_id = "dl_grid",
                             sort_by = "auc",
                             decreasing = TRUE)
d2_gridperf1 <- h2o.getGrid(grid_id = "d2_grid",
                             sort_by = "auc",
                             decreasing = TRUE)
gbm_gridperf1 <- h2o.getGrid(grid_id = "gbm_grid2",
                             sort_by = "auc",
                             decreasing = TRUE)

h2o.make_leaderboard(xgb_grid1, test_tc.hex)
h2o.make_leaderboard(xgb_grid1, data_ex.hex)

xgb_model_auc = h2o.getModel("xgb_grid1_model_13")
xgb_model_auc2 = h2o.getModel("xgb_grid1_model_34")
base_models <- list(xgb_model_auc@model_id, xgb_model_auc2@model_id)

h2o.make_leaderboard(gbm_grid2, test_tc.hex)
h2o.make_leaderboard(gbm_grid2, data_ex.hex)

gbm_model_auc = h2o.getModel("gbm_grid2_model_5")
gbm_model_auc2 = h2o.getModel("gbm_grid2_model_9")
base_models <- list(gbm_model_auc@model_id, gbm_model_auc2@model_id)

h2o.make_leaderboard(dl_grid, test_tc.hex)
h2o.make_leaderboard(dl_grid, data_ex.hex)

dl_model_auc = h2o.getModel("dl_grid_model_52")
dl_model_auc2 = h2o.getModel("dl_grid_model_50")
base_models <- list(dl_model_auc@model_id, dl_model_auc2@model_id)

h2o.make_leaderboard(d2_grid, test_tc.hex)
h2o.make_leaderboard(d2_grid, data_ex.hex)

d2_model_auc = h2o.getModel("d2_grid_model_35")
d2_model_auc2 = h2o.getModel("d2_grid_model_24")
base_models <- list(d2_model_auc@model_id, d2_model_auc2@model_id)

ensemble <- h2o.stackedEnsemble(x = features,
                                y = response,
                                training_frame = train_tc.hex,
                                validation_frame = valid_tc.hex,
                                base_models = base_models)

?h2o.stackedEnsemble
# Eval ensemble performance on a test set
perf <- h2o.performance(ensemble, newdata = test_tc.hex)
h2o.performance(ensemble, newdata = data_ex.hex)


# random foreast

set.seed(12345)
rfseed=sample(10000,100,replace=F)

res=NULL
pdf=NULL
train_tc$P9P40
for (i in 1:100){
  
  rfmod=h2o.randomForest(x = features,
                         y = "vitalstatus",
                         training_frame = train_tc.hex,
                         validation_frame = valid_tc.hex,
                         nfolds=10,
                         fold_assignment = "Stratified",
                         balance_classes = TRUE,
                         ntrees = 100, max_depth = 50,mtries = -1,sample_rate = 0.7,
                         stopping_metric = "misclassification",
                         stopping_tolerance = 0.01,
                         stopping_rounds = 3,
                         keep_cross_validation_fold_assignment = F,
                         keep_cross_validation_predictions=F,
                         score_each_iteration = TRUE,
                         seed=rfseed[i])
  
  pdftemp=predict(rfmod,test_tc.hex)%>%as_tibble()%>%mutate(.,Id=row.names(test_tc),Truth=test_tc$vitalstatus,Accuracy=ifelse(test_tc$vitalstatus ==.$predict, "Correct", "Wrong"),Model=i)
  
  vimp=rfmod@model$variable_importances%>%as_tibble()%>%mutate(.,Seed=rfseed[i],Model=i)
  
  res=rbind(res,vimp)
  
  pdf=rbind(pdf,pdftemp)
  
}

res=as_tibble(res)
res$Model=as.factor(res$Model)
pdf$Id=as.factor(pdf$Id)

pdf%>%ggplot(aes(x=Model,y=reorder(Id,Dead),fill=Dead))+geom_tile()+theme_bw()+scale_fill_gradient2(low="green",mid="gold",high="red",midpoint = 0.5)+facet_wrap(~Truth,shrink=TRUE,scale="free",ncol=1)+theme_bw(8)+theme(axis.text.x = element_text(angle =45, hjust = 1))

pdf%>%ggplot(aes(x=Model,y=reorder(Id,1-Dead),fill=Accuracy))+geom_tile()+theme_bw()+scale_fill_manual(values=c("skyblue","red"))+theme_bw(8)+theme(axis.text.x = element_text(angle =45, hjust = 1))+facet_wrap(~Truth,shrink=TRUE,scale="free",ncol=1)+scale_y_discrete("Idx")

# average predicted probability and averaged accuracy of 100 models
pdf2=NULL
getwd()
library(Hmisc)


write.csv(pdf, paste0(main.path_tc,"randomforest_1st_results3.csv"))
for (i in 1:length(rownames(test_tc))){
  Idt=pdf$Id[i]
  dtemp=subset(pdf,Id==Idt)
  sum=Hmisc::describe(dtemp)
  
  predprob=sum$Dead$counts[[5]]
  acc=sum$Accuracy$values[[2]][1]/100
  wrong=sum$Accuracy$values[[2]][2]/100
  
  pdftemp=cbind.data.frame(Id=Idt,predprob,acc,wrong)
  pdf2=rbind(pdf2,pdftemp)
}

pdf2$class=test_tc$vitalstatus
pdf2$predprob=as.numeric(as.character(pdf2$predprob))

pdf2$Id=as.factor(pdf2$Id)

pdf2$wrong[is.na(pdf2$wrong)] <- 0

pdf2%>%ggplot(aes(x=reorder(Id,predprob),y=predprob,fill=predprob))+geom_bar(alpha=0.8,stat="identity")+theme_bw(8)+scale_fill_gradient2(low="gold",mid="red",high="purple",midpoint = 0.4)+facet_wrap(~class,shrink = TRUE,ncol=1,scales = "free")+geom_hline(yintercept = 0.5,linetype=2,size=1)

pdf2%>%gather(wrong,acc,key="Accuracy",value="Prediction")%>%ggplot(aes(x=reorder(Id,-predprob),y=Prediction,fill=Accuracy))+geom_bar(alpha=0.7,position="fill",stat="Identity")+theme_bw(8)+facet_wrap(~class,shrink = TRUE,ncol=1,scales = "free")+scale_fill_manual(values=c("skyblue","red"))+scale_x_discrete("Id")

# relative importance of variables (more important a feature is, more saturated in color that feature would be)

res%>%ggplot(aes(x=reorder(Model,relative_importance),y=reorder(variable,-relative_importance),fill=relative_importance,color=relative_importance))+geom_tile(show.legend = F)+theme_bw()+scale_fill_gradient2(low="white",mid="red",high="purple",midpoint = 50)+scale_color_gradient2(low="white",mid="red",high="purple",midpoint = 50)+theme_bw(5)+theme(axis.text.x = element_text(angle =45, hjust = 1))+scale_x_discrete("Models")+scale_y_discrete("Features")


res%>%ggplot(aes(x=reorder(Model,scaled_importance),y=reorder(variable,-scaled_importance),fill=scaled_importance,color=scaled_importance))+geom_tile(show.legend = F)+theme_bw()+scale_fill_gradient2(low="white",mid="gold",high="red",midpoint = 0.5)+scale_color_gradient2(low="white",mid="gold",high="red",midpoint =0.5)+theme_bw(5)+theme(axis.text.x = element_text(angle =45, hjust = 1))+scale_x_discrete("Models")+scale_y_discrete("Features")

# extract averaged importance score

res2=NULL

for (i in 1:length(colnames(test_tc[,-length(test_tc)]))){
  dtemp=subset(res,variable==features[i])
  sum=psych::describe(dtemp)
  rei=sum[[5]][2]
  sci=sum[[5]][3]
  pct=sum[[5]][4]
  
  rtemp=cbind.data.frame(rei,sci,pct,var=dtemp[[1]][1])
  res2=rbind(res2,rtemp)
  
}

res2%>%ggplot(aes(x=reorder(var,rei),y=rei,fill=rei,color=rei))+geom_segment(aes(x=reorder(var,rei),xend=var,y=0,yend=rei),size=1,show.legend = F)+geom_point(size=1,show.legend = F)+scale_x_discrete("Features")+scale_y_continuous("Relative importance")+coord_flip()+scale_fill_gradient(low="blue",high="red")+scale_color_gradient(low="blue",high="red")+ggtitle("Relative importance")+theme_bw(5)

# extract top150 variables
res3=res2%>%arrange(desc(rei))%>%.[c(1:150),]

res3%>%ggplot(aes(x=reorder(var,rei),y=rei,fill=rei,color=rei))+geom_segment(aes(x=reorder(var,rei),xend=var,y=0,yend=rei),size=1,show.legend = F)+geom_point(size=1,show.legend = F)+scale_x_discrete("Features")+scale_y_continuous("Relative importance")+coord_flip()+scale_fill_gradient(low="blue",high="red")+scale_color_gradient(low="blue",high="red")+ggtitle("Relative importance")+theme_bw(5)

# re-train by top 150 variables
features2=as.vector(res3$var)

rfmod=h2o.randomForest(x = features2,
                       y = "vitalstatus",
                       training_frame = train_tc.hex,
                       validation_frame = valid_tc.hex,
                       nfolds=10,
                       fold_assignment = "Stratified",
                       balance_classes = TRUE,
                       ntrees = 100, max_depth = 50,mtries = -1,sample_rate = 0.8,
                       stopping_metric = "misclassification",
                       stopping_tolerance = 0.01,
                       stopping_rounds = 3,
                       keep_cross_validation_fold_assignment = F,
                       keep_cross_validation_predictions=F,
                       score_each_iteration = TRUE,
                       seed=333)

# lasso
lassomod=h2o.glm(x=features,
                 y="vitalstatus",
                 training_frame=train_tc.hex, 
                 validation_frame = valid_tc.hex,
                 nfolds = 10,seed = 546, 
                 keep_cross_validation_predictions = FALSE,
                 keep_cross_validation_fold_assignment = FALSE, 
                 fold_assignment = "Stratified",
                 family = "binomial",
                 alpha=1,
                 lambda_search = TRUE,nlambdas =20,
                 early_stopping = TRUE,
                 standardize = TRUE,
                 missing_values_handling ="Skip",
                 compute_p_values = FALSE, 
                 remove_collinear_columns = FALSE,
                 intercept = TRUE, 
                 non_negative = FALSE
)

coefs=lassomod@model$coefficients_table%>%as_tibble()%>%subset(.,coefficients!=0)
h2o.std_coef_plot(lassomod)
coefs

# Ridge
ridgemod=h2o.glm(x=features,
                 y="vitalstatus",
                 training_frame=train_tc.hex,
                 validation_frame = valid_tc.hex,
                 nfolds = 10,seed = 999, 
                 keep_cross_validation_predictions = FALSE,
                 keep_cross_validation_fold_assignment = FALSE, 
                 fold_assignment = "Stratified",
                 family = "binomial",
                 alpha=0,
                 lambda_search = TRUE,nlambdas =20,
                 early_stopping = TRUE,
                 standardize = TRUE,
                 missing_values_handling ="Skip",
                 compute_p_values = FALSE, 
                 remove_collinear_columns = FALSE,
                 intercept = TRUE, 
                 non_negative = FALSE
)
h2o.std_coef_plot(ridgemod)
#Elastic net
BiocManager::install("glmnet")
library(glmnet)

trainsetenet=train_tc%>%na.omit()
testsetenet=test_tc%>%na.omit()

trainsetenet[-1]<-lapply(trainsetenet[-1], as.numeric)
testsetenet[-1]<-lapply(testsetenet[-1], as.numeric)

x=as.matrix(trainsetenet[,-1])
y=trainsetenet[,1]
train_tc_elastic_x= as.matrix(train_tc[,-length(colnames(train_tc))])
train_tc_elastic_y= train_tc[,length(colnames(train_tc))]

test_tc_elastic_x= as.matrix(test_tc[,-length(colnames(test_tc))])
test_tc_elastic_x = test_tc[,-length(colnames(test_tc))]
test_tc_elastic_y= test_tc[,length(colnames(test_tc))]

enet = cv.glmnet(train_tc_elastic_x,train_tc_elastic_y,
                 family = "binomial",
                 type.measure = "auc",
                 nfolds = 10,
                 alpha = 0.6)

plot(enet)
enet$lambda.min
coef(enet,s = "lambda.min")%>%as.matrix()%>%as.data.frame()%>%subset(.,`s1`!=0)
predict(enet,as.matrix(test_tc[,-length(colnames(test_tc))]),s = "lambda.min")
library(mlr)

task=makeClassifTask(data=test_tc,target="vitalstatus")

dummylrner=makeLearner("classif.rpart",predict.type ="prob")
dummymod=train(dummylrner,task)
dummypred=predict(dummymod,task)

LASSOPRED=predict(lassomod,test_tc.hex)%>%as.data.frame()
LASSOPRED$cal_predic = ifelse(LASSOPRED$Alive > LASSOPRED$Dead, "Alive", "Dead")
LASSOPRED=LASSOPRED%>%mutate(.,id=c(1:length(rownames(test_tc))),truth=test_tc$vitalstatus,prob.Alive=.[,2],prob.Dead=.[,3],response=.$cal_predic)%>%.[,-c(1:3)]


RIDGEPRED=predict(ridgemod,test_tc.hex)%>%as.data.frame()
RIDGEPRED=RIDGEPRED%>%mutate(.,id=c(1:length(rownames(test_tc))),truth=test_tc$vitalstatus,prob.Alive=.[,2],prob.Dead=.[,3],response=.$predict)%>%.[,-c(1:3)]
remove(ENETPRED)
remove(enetclass)
enetclass=predict(enet,as.matrix(test_tc[,-length(colnames(test_tc))]),type="class",s = "lambda.min")%>%as.data.frame()
ENETPRED=predict(enet,as.matrix(test_tc[,-length(colnames(test_tc))]) ,s = "lambda.min")%>%as.data.frame()%>%mutate(.,id=c(1:length(rownames(test_tc))),truth=test_tc$vitalstatus,prob.Alive=1-.[,1],prob.Dead=.[,1],response=enetclass[,1])%>%.[,-1]

RFPRED=predict(rfmod,test_tc.hex)%>%as.data.frame()
RFPRED$cal_predic = ifelse(RFPRED$Alive > RFPRED$Dead, "Alive", "Dead")
RFPRED=RFPRED%>%mutate(.,id=c(1:length(rownames(test_tc))),truth=test_tc$vitalstatus,prob.Alive=.[,2],prob.Dead=.[,3],response=.$cal_predic)%>%.[,-c(1:3)]

pred1=dummypred
pred1$data=LASSOPRED
pred2=dummypred
pred2$data=ENETPRED
pred3=dummypred
pred3$data=RFPRED
pred4=dummypred
pred4$data=RIDGEPRED%>%na.omit()

measure=list(mlr::bac,mlr::auc,mlr::acc,mlr::f1,mlr::ppv,mlr::npv,mlr::tnr,mlr::tpr)
p1 = NULL
p2 = NULL
p3 = NULL
p4 = NULL

p1=performance(pred1,measure)
p2=performance(pred2,measure)
p3=performance(pred3,measure)
p4=performance(pred4,measure)

lasso=p1%>%as.vector()
enet=p2%>%as.vector()
rf=p3%>%as.vector()
rdg=p4%>%as.vector()

labs=c("BAC","AUC","ACC","F1","PPV","NPV","TPR","TNR")

scores=list("Lasso"=lasso*10,"Elasticnet"=enet*10,"RandomForest"=rf*10,"Ridge"=rdg*10)
scores=list("Elasticnet"=enet*10,"RandomForest"=rf*10,"Ridge"=rdg*10)
scores=list("Lasso"=lasso*10,"RandomForest"=rf*10,"Ridge"=rdg*10)
library(radarchart)

chartJSRadar(scores = scores, labs = labs, maxScale = 10)

measure2=list(ber,brier,fnr,fpr,mmce,fdr)

e1 = NULL
e2 = NULL
e3 = NULL
e4 = NULL

e1=performance(pred1,measure2)
e2=performance(pred2,measure2)
e3=performance(pred3,measure2)
e4=performance(pred4,measure2)

lasso=e1%>%as.vector()
enet=e2%>%as.vector()
rf=e3%>%as.vector()
rdg=e4%>%as.vector()

labs=c("BER","BRIER","FNR","FPR","MMCE","FDR")

scores=list("Lasso"=lasso,"Elasticnet"=enet,"RandomForest"=rf,"Ridge"=rdg)
scores=list("Elasticnet"=enet,"RandomForest"=rf,"Ridge"=rdg)
scores=list("Lasso"=lasso,"RandomForest"=rf,"Ridge"=rdg)
chartJSRadar(scores = scores, labs = labs, maxScale =0.8)


### fin 

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
