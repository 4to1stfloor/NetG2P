
library(survival)
library(gbm)
# set random state
set.seed(2891)

cancer_ex_data = readRDS("/home/seokwon/nas/00.data/target_cancer_folder/TCGA-CESC.rds")


# 2) transfer data for model build
# A - data split

data_tc_ori_alive = cancer_ex_data[cancer_ex_data$vital_status == "Alive",]
data_tc_ori_dead = cancer_ex_data[cancer_ex_data$vital_status == "Dead",]

train_alive = sample(1:nrow(data_tc_ori_alive), nrow(data_tc_ori_alive)*0.6) 
train_dead = sample(1:nrow(data_tc_ori_dead), nrow(data_tc_ori_dead)*0.6) 
testtmp_val_alive = -train_alive
testtmp_val_dead = -train_dead

train_tc_alive = data_tc_ori_alive[train_alive,]
train_tc_dead = data_tc_ori_dead[train_dead,]

train_tc = rbind(train_tc_alive,train_tc_dead)
train_tc$time = train_tc$duration
train_tc$vital_status = NULL
train_tc$duration = NULL

testtmp_tc_alive = data_tc_ori_alive[testtmp_val_alive,]
testtmp_tc_dead = data_tc_ori_dead[testtmp_val_dead,]

valid_val_alive = sample(1:nrow(testtmp_tc_alive), nrow(testtmp_tc_alive)*0.5) 
valid_val_dead = sample(1:nrow(testtmp_tc_dead), nrow(testtmp_tc_dead)*0.5) 

test_val_alive = -valid_val_alive
test_val_dead = -valid_val_dead

valid_tc_alive = testtmp_tc_alive[valid_val_alive,]
valid_tc_dead = testtmp_tc_dead[valid_val_dead,]
valid_tc = rbind(valid_tc_alive,valid_tc_dead)
valid_tc$time = valid_tc$duration
valid_tc$vital_status = NULL
valid_tc$duration = NULL
train_tc = rbind(train_tc , valid_tc)

test_tc_alive = testtmp_tc_alive[test_val_alive,]
test_tc_dead = testtmp_tc_dead[test_val_dead,]
test_tc = rbind(test_tc_alive,test_tc_dead)

train_tc <- train_tc[, c(2:ncol(train_tc), 1)]

# using training set fits gbm model
model <- gbm(Surv(time, status) ~ .,
             distribution='coxph',
             data=train_tc,
             n.trees=3000,
             shrinkage=0.005,
             interaction.depth=1,
             n.minobsinnode=5,
             # cv.folds=5,
             verbose = T,
             n.cores = 20,)

# values of loss function on training set for each tree
print(model$train.error[2901:3000])

best.iter <- gbm.perf(model, plot.it = TRUE, method = 'cv')

# return a vector of prediction on n.trees indicting log hazard scale.f(x)
pred.train <- predict(model, data_train, n.trees = best.iter)
pred.test <- predict(model, data_test, n.trees = best.iter)

# CI
Hmisc::rcorr.cens(-pred.train, Surv(data_train$time, data_train$status))
Hmisc::rcorr.cens(-pred.test, Surv(data_test$time, data_train$status))


# Sepecify time of interest
time.interest <- sort(unique(data_train$time[data_train$status==1]))
# Estimate the cumulative baseline hazard function using training data
basehaz.cum <- basehaz.gbm(data_train$time,data_train$status, pred.train, t.eval = time.interest, cumulative = TRUE)
surf.i <- exp(-exp(pred.test[1])*basehaz.cum)
print(surf.i)

specif.time <- time.interest[10]
cat("Survival Rate of all at time", specif.time, "\n")
surv.rate <- exp(-exp(pred.test)*basehaz.cum[10])
print(surv.rate)
?basehaz.gbm

res_test <- data_test
# predicted outcome for test set
res_test$pred = -pred.test
res_test$survival_rate <- surv.rate
# Save data
write.csv(res_test, file = "result_gbm.csv")
