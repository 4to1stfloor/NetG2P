
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
library(nnet)
require(RSNNS)
require(clusterGeneration)

# 0) call input

filepath = "/home/seokwon/"
ref_path = paste0(filepath, "99.reference/")
kindofCancer = "brca"

Cancername = "TCGA-BRCA" # you can change this
main.path_tc = paste0(filepath, "00.data/","30.", Cancername,"/")

ex_path = paste0(filepath,"01.externaldata/", kindofCancer, "/") 
CancerType_in = "brca_metabric" # you can change this
main.path_in = paste0(ex_path,CancerType_in, "/")

CancerType = "TCGA-BRCA"

# 1) prepare for NN
# (1) for TCGA data upload

#   A. for each pathway
surv.pvals = readRDS(paste0(main.path_tc,"/surv_pvals_54_KEGG_GC_20220117.rds"))
phyper = readRDS(paste0(main.path_tc, "/phyper_enrichment_54_KEGG_GC_20220117.rds"))

#   B. for pathwaylink
surv.pvals = readRDS(paste0(main.path_tc,"/surv_pvals_54_KEGG_net_GC_20220117_filt5.rds"))
phyper = readRDS(paste0(main.path_tc, "/phyper_enrichment_54_KEGG_net_GC_20220117.rds"))

mut.mat = readRDS(paste0(main.path_tc,"/", CancerType,"_mut_count_filt_data.rds"))
mut.count = as.data.frame(colSums(mut.mat))
colnames(mut.count) = "counts"

# (1) - clinical data upload for kidney

cli = fread(paste0(ref_path , "all_clin_indexed.csv"))
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

cli_surv$status = NA
cli_surv$status[which(cli_surv$vital_status == "Dead")] = 1
cli_surv$status[which(cli_surv$vital_status == "Alive")] = 0

cli_surv$status = as.numeric(cli_surv$status)
cli_surv


# (1) - filtering only 01A

mut.count$patient_id = substr(rownames(mut.count), 1, 16)
mut.count$forfilt =  rownames(mut.count)
mut.count.filt = mut.count %>% filter(str_detect(forfilt, "-01A-")) 
mut.filt = mut.count.filt[which(mut.count.filt$counts < 1000),]
mut.filt$forfilt = NULL

int.id = intersect(colnames(phyper), mut.filt$patient_id)
idx.id = match(int.id, colnames(phyper))
phyper.filt = phyper[,idx.id]

#  transfer data 1 (>0.05) or 100 (<0.05)

phyper.filt_t = t(phyper.filt)
phyper.filt_t_num = phyper.filt_t[,-which(colSums(phyper.filt_t) == dim(phyper.filt_t)[1])]

phyper.filt_t_num.bi = ifelse(phyper.filt_t_num >0.05 , 0 , 1)
phyper.filt_t_num.bi = ifelse(phyper.filt_t_num >0.05 , 1 , 1000000)

rownames(phyper.filt_t_num.bi) = substr(rownames(phyper.filt_t_num.bi), 1,12)

df.phyper.filt_t_log_num = as.data.frame(phyper.filt_t_num.bi)
df.phyper.filt_t_log_num$submitter_id = rownames(df.phyper.filt_t_log_num)
df.phyper.filt_t_log_num

df.phyper.merge_log_all = as.data.frame(merge(cli_surv , df.phyper.filt_t_log_num, by = "submitter_id"))
rownames(df.phyper.merge_log_all) = df.phyper.merge_log_all$submitter_id

df.phyper.merge_log_all

df.phyper.merge_log_all3 = df.phyper.merge_log_all[,c(-1,-3,-4,-5,-6,-7)]
df.phyper.merge_log_all3
df.phyper.merge_vital_log_all2= df.phyper.merge_log_all3

colnames(df.phyper.merge_vital_log_all2)= gsub("[[:punct:]]", "", colnames(df.phyper.merge_vital_log_all2))
df.phyper.merge_vital_log_all2= df.phyper.merge_vital_log_all2 %>% relocate("vitalstatus", .after = "P53P54")
# df.phyper.merge_vital_log_all2= df.phyper.merge_vital_log_all2 %>% relocate("vitalstatus", .after = "P54")

saveRDS(df.phyper.merge_vital_log_all2, file = paste0(main.path_tc,"/",CancerType,"_pathwaylink_all_1or0.rds"))

saveRDS(df.phyper.merge_vital_log_all2, file = paste0(main.path_tc,"/",CancerType,"_pathwaylink_all_1or1000000.rds"))

# (1) - transform -log p.val and transverse

phyper.filt_t = t(phyper.filt)
phyper.filt_t_num= phyper.filt_t[,-which(colSums(phyper.filt_t) == dim(phyper.filt_t)[1])]

phyper.filt_t_log = -log(phyper.filt_t_num)

phyper.filt_t_log = t(phyper.filt_t_log)
colnames(phyper.filt_t_log) = substr(colnames(phyper.filt_t_log), 1,12)
df.phyper.filt_t_log = as.data.frame(phyper.filt_t_log)
df.phyper.filt_t_log_num = data.frame(matrix(nrow = length(rownames(df.phyper.filt_t_log))))
for (all_sub in cli_surv$submitter_id) {
  if (all_sub %in% colnames(df.phyper.filt_t_log)) {
    # df.phyper.all.link_num = cbind(df.phyper.all.link_num,df.phyper.filt_num[,all_sub])
    tmp.phyer.filt_num = as.data.frame(df.phyper.filt_t_log[,all_sub])
    colnames(tmp.phyer.filt_num) = all_sub
    rownames(tmp.phyer.filt_num) = rownames(df.phyper.filt_t_log)
    # print(tmp.phyer.filt_num)
    df.phyper.filt_t_log_num = cbind(df.phyper.filt_t_log_num,tmp.phyer.filt_num)
    tmp.phyer.filt_num = NULL
  }
  
} 
df.phyper.filt_t_log_num
df.phyper.filt_t_log_num= df.phyper.filt_t_log_num[,-1]

df.phyper.filt_t_log_num = t(df.phyper.filt_t_log_num)

df.phyper.filt_t_log_num = as.data.frame(df.phyper.filt_t_log_num)
df.phyper.filt_t_log_num$submitter_id = rownames(df.phyper.filt_t_log_num)
df.phyper.filt_t_log_num

df.phyper.merge_log_all = as.data.frame(merge(cli_surv , df.phyper.filt_t_log_num, by = "submitter_id"))
rownames(df.phyper.merge_log_all) = df.phyper.merge_log_all$submitter_id

df.phyper.merge_log_all

df.phyper.merge_log_all3 = df.phyper.merge_log_all[,c(-1,-3,-4,-5,-6,-7)]
df.phyper.merge_log_all3
df.phyper.merge_vital_log_all2= df.phyper.merge_log_all3

colnames(df.phyper.merge_vital_log_all2)= gsub("[[:punct:]]", "", colnames(df.phyper.merge_vital_log_all2))
df.phyper.merge_vital_log_all2= df.phyper.merge_vital_log_all2 %>% relocate("vitalstatus", .after = "P53P54")
# df.phyper.merge_vital_log_all2= df.phyper.merge_vital_log_all2 %>% relocate("vitalstatus", .after = "P54")

saveRDS(df.phyper.merge_vital_log_all2, file = paste0(main.path_tc,"/",CancerType,"_pathwaylink_all.rds"))

# (2) - for external data

main.path_in 

#   A. for each pathway
surv.pvals_in = readRDS(paste0(main.path_in,"/surv_pvals_54_KEGG_GC_20220117.rds"))
phyper_in = readRDS(paste0(main.path_in, "/phyper_enrichment_54_KEGG_GC_20220117.rds"))

#   B. for pathwaylink
surv.pvals_in = readRDS(paste0(main.path_in,"/surv_pvals_54_KEGG_net_GC_20220117_filt5.rds"))
phyper_in = readRDS(paste0(main.path_in, "/phyper_enrichment_54_KEGG_net_GC_20220117.rds"))

mut.mat_in = readRDS(paste0(main.path_in, CancerType_in,"_mut_count_filt_data.rds"))
mut.count_in = as.data.frame(colSums(mut.mat_in))
colnames(mut.count_in) = "counts"


# (2) - clinical data upload

cli_ex = fread(paste0(main.path_in,CancerType_in,"_clinical_data.tsv")) 
cli_ex
cli_ex_surv = cli_ex[,
                     c("Patient ID",
                       "Overall Survival (Months)",
                       "Overall Survival Status")]

cli_ex_surv$vital_status = NA
cli_ex_surv$vital_status[which(cli_ex_surv$`Overall Survival Status` == "1:DECEASED")] = "Dead"
cli_ex_surv$vital_status[which(cli_ex_surv$`Overall Survival Status` == "0:LIVING")] = "Alive"
cli_ex_surv$status = NA
cli_ex_surv$status[which(cli_ex_surv$`Overall Survival Status` == "1:DECEASED")] = 1
cli_ex_surv$status[which(cli_ex_surv$`Overall Survival Status` == "0:LIVING")] = 0

cli_ex_surv = cli_ex_surv[!is.na(cli_ex_surv$`Overall Survival Status`),]
cli_ex_surv$overall_survival = cli_ex_surv$`Overall Survival (Months)`*30
cli_ex_surv$status = as.numeric(cli_ex_surv$status)
cli_ex_surv

# (2) - filtering only 01A 

mut.count_in$patient_id = rownames(mut.count_in)
mut.filt_in = mut.count_in[which(mut.count_in$counts < 1000),]

int.id_in = intersect(colnames(phyper_in), mut.filt_in$patient_id)
idx.id_in = match(int.id_in, colnames(phyper_in))
phyper.filt_in = phyper_in[,idx.id_in]

# (2) -  transfer data 1 (>0.05) or 100 (<0.05)

phyper.filt_in_t = t(phyper.filt_in)
phyper.filt_in_t_num= phyper.filt_in_t[,-which(colSums(phyper.filt_in_t) == dim(phyper.filt_in_t)[1])]

phyper.filt_in_t_num.bi = ifelse(phyper.filt_in_t_num >0.05 , 0 , 1)
phyper.filt_in_t_num.bi = ifelse(phyper.filt_t_num >0.05 , 1 , 1000000)

rownames(phyper.filt_in_t_num.bi) = substr(rownames(phyper.filt_in_t_num.bi), 1,16)

df.phyper.filt_in_t_log_num = as.data.frame(phyper.filt_in_t_num.bi)
df.phyper.filt_in_t_log_num$submitter_id = rownames(df.phyper.filt_in_t_log_num)
df.phyper.filt_in_t_log_num

cli_ex_surv$submitter_id = cli_ex_surv$`Patient ID`

df.phyper.merge_in_log_all = as.data.frame(merge(cli_ex_surv , df.phyper.filt_in_t_log_num, by = "submitter_id"))
rownames(df.phyper.merge_in_log_all) = df.phyper.merge_in_log_all$submitter_id

df.phyper.merge_log_in3 = df.phyper.merge_in_log_all[,c(-1,-2,-3,-4,-6,-7)]
df.phyper.merge_log_in3
df.phyper.merge_vital_log_in2= df.phyper.merge_log_in3

colnames(df.phyper.merge_vital_log_in2)= gsub("[[:punct:]]", "", colnames(df.phyper.merge_vital_log_in2))
df.phyper.merge_vital_log_in2= df.phyper.merge_vital_log_in2 %>% relocate("vitalstatus", .after = "P53P54")
# df.phyper.merge_vital_log_all2= df.phyper.merge_vital_log_all2 %>% relocate("vitalstatus", .after = "P54")

saveRDS(df.phyper.merge_vital_log_in2, file = paste0(main.path_in,CancerType_in,"_pathwaylink_all_0or1.rds"))
saveRDS(df.phyper.merge_vital_log_in2, file = paste0(main.path_in,CancerType_in,"_pathwaylink_all_1or1000000.rds"))

# (2) - transform -log p.val and transverse

phyper.filt_in_t = t(phyper.filt_in)
phyper.filt_in_t_num= phyper.filt_in_t[,-which(colSums(phyper.filt_in_t) == dim(phyper.filt_in_t)[1])]

phyper.filt_in_t_log = -log(phyper.filt_in_t_num)
phyper.filt_in_t_log = t(phyper.filt_in_t_log)

df.phyper.filt_in_t_log = as.data.frame(phyper.filt_in_t_log)
df.phyper.filt_in_t_log_num = data.frame(matrix(nrow = length(rownames(df.phyper.filt_in_t_log))))
for (all_sub in cli_sclc_uco_surv$`Patient ID`) {
  if (all_sub %in% colnames(df.phyper.filt_in_t_log)) {
    # df.phyper.all.link_num = cbind(df.phyper.all.link_num,df.phyper.filt_num[,all_sub])
    tmp.phyer.filt_num = as.data.frame(df.phyper.filt_in_t_log[,all_sub])
    colnames(tmp.phyer.filt_num) = all_sub
    rownames(tmp.phyer.filt_num) = rownames(df.phyper.filt_in_t_log)
    # print(tmp.phyer.filt_num)
    df.phyper.filt_in_t_log_num = cbind(df.phyper.filt_in_t_log_num,tmp.phyer.filt_num)
    tmp.phyer.filt_num = NULL
  }
  
} 
df.phyper.filt_in_t_log_num
df.phyper.filt_in_t_log_num= df.phyper.filt_in_t_log_num[,-1]

df.phyper.filt_in_t_log_num = t(df.phyper.filt_in_t_log_num)

df.phyper.filt_in_t_log_num = as.data.frame(df.phyper.filt_in_t_log_num)
df.phyper.filt_in_t_log_num$submitter_id = rownames(df.phyper.filt_in_t_log_num)
df.phyper.filt_in_t_log_num

cli_sclc_uco_surv$submitter_id = cli_sclc_uco_surv$`Patient ID`

df.phyper.merge_in_log_all = as.data.frame(merge(cli_sclc_uco_surv , df.phyper.filt_in_t_log_num, by = "submitter_id"))
rownames(df.phyper.merge_in_log_all) = df.phyper.merge_in_log_all$submitter_id

# (2) - Alive or Dead

df.phyper.merge_in_log_all_filt = df.phyper.merge_in_log_all

df.phyper.merge_log_in3 = df.phyper.merge_in_log_all_filt[,c(-1,-2,-3,-4,-6,-7)]
df.phyper.merge_log_in3
df.phyper.merge_vital_log_in2= df.phyper.merge_log_in3

colnames(df.phyper.merge_vital_log_in2)= gsub("[[:punct:]]", "", colnames(df.phyper.merge_vital_log_in2)) # '-' does not allow
df.phyper.merge_vital_log_in2= df.phyper.merge_vital_log_in2 %>% relocate("vitalstatus", .after = "P53P54")
df.phyper.merge_vital_log_in2= df.phyper.merge_vital_log_in2 %>% relocate("vitalstatus", .after = "P54")
df.phyper.merge_vital_log_in2 # this will be input data for validate

saveRDS(df.phyper.merge_vital_log_in2, file = paste0(main.path_in, CancerType_in,"_pathwaylink_all.rds"))