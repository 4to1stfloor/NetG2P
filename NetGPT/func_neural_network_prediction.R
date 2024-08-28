library(survminer) 
library(ggplot2) 
library(pheatmap)
library(survival) 
library(data.table) 
library(RColorBrewer) 
library(NbClust)
library(cluster) 
library(factoextra) 
library(dplyr)
library(stringr)

CancerType = "TCGA-COAD"

filepath ="C:/Users/user/Desktop/data/"

surv.pvals = readRDS(paste0(filepath, CancerType,"/surv_pvals_54_KEGG_net_GC_20220117_filt5.rds"))
phyper = readRDS(paste0(filepath, CancerType, "/phyper_enrichment_54_KEGG_net_GC_20220117.rds"))

mut.mat = readRDS(paste0(filepath , CancerType,"/", CancerType,"_mut_count_filt_data.rds"))
mut.count = as.data.frame(colSums(mut.mat))
colnames(mut.count) = "counts"

########
## clinical data upload
####

cli = fread("C:/Users/user/Desktop/data/99.reference/all_clin_indexed.csv")
cli_surv = cli[cli$project == CancerType,
               c("submitter_id",
                 "vital_status",
                 "days_to_death",
                 "days_to_last_follow_up")]

cli_surv$deceased = cli_surv$vital_status == "Dead"
cli_surv$overall_survival = ifelse(cli_surv$deceased,
                                   cli_surv$days_to_death,
                                   cli_surv$days_to_last_follow_up)

cli_surv = cli_surv[!is.na(cli_surv$vital_status),]

cli_surv$status[which(cli_surv$vital_status == "Dead")] = 1
cli_surv$status[which(cli_surv$vital_status == "Alive")] = 0

cli_surv$status = as.numeric(cli_surv$status)
cli_surv
########################
#### 01A filtering #####
########################

mut.count$patient_id = substr(rownames(mut.count), 1, 16)
mut.count$forfilt =  rownames(mut.count)
mut.count.filt = mut.count %>% filter(str_detect(forfilt, "-01A-")) 
mut.filt = mut.count.filt[which(mut.count.filt$counts < 1000),]
mut.filt$forfilt = NULL

# dim(mut.filt)
# colnames(phyper)
mut.filt
phyper
int.id = intersect(colnames(phyper), mut.filt$patient_id)
idx.id = match(int.id, colnames(phyper))
phyper.filt = phyper[,idx.id]

surv = surv.pvals
surv$hallmark = rownames(surv)
surv.sig = surv$hallmark[which(surv[,1] < 0.05)]
surv.sig

##########################################################
#### below code make pvals for 0 or 1 to make binary matrix
###########################################################

# for (n in 1:length(surv.sig)) {
#   phyper.each = phyper.filt[which(rownames(phyper.filt) == surv.sig[n]), ]
#   phyper.dat.tmp = as.data.frame(cbind(id = substr(names(phyper.each), 1, 12), pvals = as.numeric(phyper.each)))
#   phyper.dat.tmp$group = ""
#   phyper.dat.tmp$group[which(phyper.dat.tmp$pvals <= 0.05)] = 0
#   phyper.dat.tmp$group[which(phyper.dat.tmp$pvals > 0.05)] = 1
#   if (n == 1) {
#     phyper.dat.bm= as.numeric(phyper.dat.tmp$group)
#   } else {
#     phyper.dat.bm = rbind(phyper.dat.bm, as.numeric(phyper.dat.tmp$group))
#   }
# }
# 
# rownames(phyper.dat.bm) = surv.sig
# colnames(phyper.dat.bm) = substr(colnames(phyper.filt), 1,12)
# 
# dim(phyper.dat.bm)

for (n in 1:length(surv.sig)) {
  phyper.each = phyper.filt[which(rownames(phyper.filt) == surv.sig[n]), ]
  phyper.dat = as.data.frame(cbind(id = substr(names(phyper.each), 1, 12), pvals = as.numeric(phyper.each)))
  
  if (n == 1) {
    phyper.dat.2= as.numeric(phyper.dat$pvals)
  } else {
    phyper.dat.2 = rbind(phyper.dat.2, as.numeric(phyper.dat$pvals))
  }
}

phyper.dat.2
#phyper.dat.2 = t(as.matrix(phyper.dat.2))
rownames(phyper.dat.2) = surv.sig
colnames(phyper.dat.2) = substr(colnames(phyper.filt),1,12)


##########################################################################
#################### test ################################################
##########################################################################



######
#####
#####
df.phyper.2 = as.data.frame(phyper.dat.2)

df.phyper.all= data.frame(matrix(nrow = length(rownames(phyper.dat.2))))
for (all_sub in cli_surv$submitter_id) {
  if (all_sub %in% colnames(phyper.dat.2)) {
    df.phyper.all = cbind(df.phyper.all,df.phyper.2[all_sub])
  }
  
}  
df.phyper.all= df.phyper.all[,-1]
df.phyper.all_t = t(df.phyper.all)
df.phyper.all_t = as.data.frame(df.phyper.all_t)
df.phyper.all_t$submitter_id = rownames(df.phyper.all_t)
df.phyper.all_t

df.phyper.merge = as.data.frame(merge(cli_surv , df.phyper.all_t, by = "submitter_id"))
rownames(df.phyper.merge) = df.phyper.merge$submitter_id

#######
### pre-treated for nn (status 0 or 1)
#####

df.phyper.merge2 = df.phyper.merge[,c(-1,-2,-3,-4,-5,-6)]
df.phyper.merge2= df.phyper.merge2 %>% relocate("status", .after = "P52-P53")
colnames(df.phyper.merge2)= gsub("[[:punct:]]", "", colnames(df.phyper.merge2))

df.phyper.merge2
# Training and Test Data

#######
### pre-treated for nn (status 0 or 1)
#####

trainset <- df.phyper.merge2[1:236, ]
testset <- df.phyper.merge2[237:length(df.phyper.merge2), ]


##############################
###### for neural net 
###############################
df.phyper.merge_vital = df.phyper.merge[,c(-1,-3,-4,-5,-6,-7)]
df.phyper.merge_vital2= df.phyper.merge_vital %>% relocate("vital_status", .after = "P52-P53")
colnames(df.phyper.merge_vital2)= gsub("[[:punct:]]", "", colnames(df.phyper.merge_vital2))

df.phyper.merge_vital2
# Training and Test Data
trainset_v <- df.phyper.merge_vital2[1:236, ]
testset_v <- df.phyper.merge_vital2[237:length(df.phyper.merge2), ]
library(neuralnet)

nn.cli <- neuralnet(status~., data=trainset, hidden=c(10,2), linear.output=FALSE, threshold=0.001)
nn.cli$result.matrix
plot(nn.cli)

# Test the resulting output for neuralnet function 
temp_test <- testset[,-23]

head(temp_test)
nn.results <- compute(nn.cli, temp_test)
predicted_result = nn.results$net.result
cor(predicted_result , testset$status)

results <- data.frame(actual = testset$status, prediction = nn.results$net.result)
results

roundedresults<-sapply(results,round,digits=0)
roundedresultsdf=data.frame(roundedresults)
attach(roundedresultsdf)
table(actual,prediction)

#4. NEURAL NETWORK
nn.cli2 <- neuralnet(status~., data=trainset, hidden=c(10,5), linear.output=FALSE, threshold=0.001)
nn.cli2$result.matrix
plot(nn.cli2)

nncli2.results <- compute(nn.cli2, temp_test)
results.nncli2 <- data.frame(actual = testset$status, prediction = nncli2.results$net.result)
results.nncli2

# table
roundedresults.cli2<-sapply(results.nncli2,round,digits=0)
roundedresultsdf.cli2=data.frame(roundedresults.cli2)
attach(roundedresultsdf.cli2)
table(actual,prediction)


##################################
##### for nnet 
################################
library(nnet)
require(RSNNS)
require(clusterGeneration)

# names(trainset_y)<-'Y1'

trainset_x = trainset[,-length(colnames(trainset))]
trainset_y = as.data.frame(trainset[,length(colnames(trainset))])
colnames(trainset_y) = "status"
rownames(trainset_y) = rownames(trainset)
# make model
cli_nn <-nnet(trainset_x,trainset_y,data=trainset,size=25,linout=T, rang = 0.1,decay = 5e-4, maxit = 200)
cli_nn <-nnet(trainset_x,trainset_y,data=trainset,size=20,linout=T)
cli_mlp = mlp(trainset_vx, trainset_vy, size = 8 , linOut = T)
####
###### validate model for nn
####


testset_x = testset[,-length(colnames(testset))]

results <- data.frame(actual = testset$status, prediction = predict(cli_nn, testset_x))

roundedresults<-sapply(results,round,digits=0)
roundedresultsdf=data.frame(roundedresults)
roundedresultsdf= roundedresultsdf %>% mutate(
  actual = ifelse( actual <0 , 0 ,actual),
  prediction = ifelse(prediction <0 , 0 , prediction )
)
attach(roundedresultsdf)
result.table= table(actual,prediction)
write.csv(result.table, "results_table_for_ann_cli.csv")

mod2<-lm(status~.,data=cbind(trainset_y[,'status',drop=F],trainset_x))

#import sensitivity analysis function
source_url('https://gist.githubusercontent.com/fawda123/6860630/raw/b8bf4a6c88d6b392b1bfa6ef24759ae98f31877c/lek_fun.r')

lek.fun(cli_nn)
lek.fun(cli_mlp, exp.in = trainset_vx)
summary(mod2)
###
###### for calculate weight or draw plot
####


df.phyper.merge$Alive = ifelse(df.phyper.merge$vital_status == "Alive",
                          1,0)
df.phyper.merge$Dead = ifelse(df.phyper.merge$vital_status == "Dead",
                               1,0)
df.phyper.merge3 = df.phyper.merge[,c(-1,-2,-3,-4,-5,-6,-7)]
df.phyper.merge3
df.phyper.merge_vital2= df.phyper.merge3
trainset_v <- df.phyper.merge_vital2[1:236, ]
testset_v <- df.phyper.merge_vital2[237:length(df.phyper.merge_vital2), ]

trainset_vx = trainset_v[,-c(length(colnames(trainset_v))-1,length(colnames(trainset_v)))]
trainset_vy = trainset_v[,c(length(colnames(trainset_v))-1,length(colnames(trainset_v)))]
testset_vx = testset_v[,-c(length(colnames(testset_v))-1,length(colnames(testset_v)))]
trainset_vy= as.data.frame(trainset_vy)
# make model
cli_nnv <-nnet(trainset_vx,trainset_vy,data=trainset_v,size=20,linout=T, rang = 0.1,decay = 5e-4, maxit = 200)

# validate
library(dplyr)
predic_df = as.data.frame(predict(cli_nnv, testset_vx))
results_alive <- data.frame(actual_alive = testset_v$Alive, prediction_alive = predic_df$Alive)
results_Dead <- data.frame(actual_dead = testset_v$Dead, prediction_dead = predic_df$Dead)

roundedresults_alive<-sapply(results_alive,round,digits=0)
roundedresults_dead<-sapply(results_Dead,round,digits=0)

roundedresultsdf_alive=data.frame(roundedresults_alive)
roundedresultsdf_dead=data.frame(roundedresults_dead)

roundedresultsdf_alive= roundedresultsdf_alive %>% mutate(
  actual_alive = ifelse( actual_alive <0 ,
                         0 ,actual_alive),
  actual_alive = ifelse( actual_alive >1 , 
                         1 ,actual_alive),
  prediction_alive = ifelse(prediction_alive <0 , 
                            0 , prediction_alive ),
  prediction_alive = ifelse(prediction_alive >1 , 
                            1 , prediction_alive )
)

roundedresultsdf_dead= roundedresultsdf_dead %>% mutate(
  actual_dead = ifelse( actual_dead <0 ,
                         0 ,actual_dead),
  actual_dead = ifelse( actual_dead >1 , 
                         1 ,actual_dead),
  prediction_dead = ifelse(prediction_dead <0 , 
                            0 , prediction_dead ),
  prediction_dead = ifelse(prediction_dead >1 , 
                            1 , prediction_dead )
)

attach(roundedresultsdf_alive)
result.alive.table= table(actual_alive,prediction_alive)
write.csv(result.alive.table, "results_alive_table_for_ann_cli.csv")

attach(roundedresultsdf_dead)
result.dead.table= table(actual_dead,prediction_dead)
write.csv(result.dead.table, "results_dead_table_for_ann_cli.csv")


##### plot 

plot.nnet(cli_nnv,pos.col='darkgreen',neg.col='darkblue',alpha.val=0.7,rel.rsc=15,
          circle.cex=10,cex=1.4,
          circle.col='brown')

# find weight
plot.nnet(cli_nnv, wts.only=T)

# plot.nnet(cli_nn,pos.col='darkgreen',neg.col='darkblue',alpha.val=0.7,rel.rsc=15,
#           circle.cex=10,cex=1.4,circle.col='brown',all.in='Sepal W.',all.out='Alive')

plot.nnet(cli_nnv,pos.col='darkgreen',neg.col='darkblue',alpha.val=0.7,rel.rsc=15,
          circle.cex=10,cex=1.4,circle.col='brown',all.out='Alive')
          

##
source('https://gist.githubusercontent.com/fawda123/6860630/raw/b8bf4a6c88d6b392b1bfa6ef24759ae98f31877c/lek_fun.r')

library(NeuralNetTools)
garson(cli_nnv) # it will not work at two output mode
lekprofile(cli_nnv)



#relative importance function
library(devtools)
source_url('https://gist.github.com/fawda123/6206737/raw/2e1bc9cbc48d1a56d2a79dd1d33f414213f5f1b1/gar_fun.r')

#relative importance of input variables for Alive or Dead. you can change the str
rel.imp<-gar.fun('Dead',cli_nnv , bar.plot=F)$rel.imp

#color vector based on relative importance of input values
# the value of 22 meaning is the number of input column
cols<-colorRampPalette(c('green','red'))(22)[rank(rel.imp)]

##
#plotting function
source_url('https://gist.githubusercontent.com/fawda123/7471137/raw/466c1474d0a505ff044412703516c34f1a4684a5/nnet_plot_update.r')

# plot model with new color vector
# separate colors for input vectors using a list for 'circle.col'
# all.out does not automatically matched with the results of rel.imp. 
# So you can manually change the all.out value correlate with rel.imp

plot(cli_nnv,circle.col=list(cols,'lightblue'), all.out = 'Dead')

