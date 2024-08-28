# call library

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

# call input or pathway

filepath = "/home/seokwon/nas/"
Cancerlist = dir(paste0(filepath, "/00.data/filtered_TCGA/"))

num_CancerType = "04.TCGA-CESC"
for (num_CancerType in Cancerlist) {
  
  # 1) prepare for NN
  # (1) for TCGA data upload
  
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/")
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  Cancername = gsub("TCGA-" , "" , CancerType)
  
  phyper_net = readRDS(paste0(main.path_tc, num_CancerType, "/phyper_cellline_enrichment_54_KEGG_net_GC_20220117.rds"))
  phyper_each = readRDS(paste0(main.path_tc, num_CancerType,"/phyper_cellline_enrichment_54_KEGG_GC_20220117.rds"))
  
  if (all.equal(colnames(phyper_net) , colnames(phyper_each))) {
    phper_dual = rbind(phyper_net, phyper_each)
  }
  
  phper_each_log = as.data.frame(t(-log(phyper_each)))
  phper_net_log = as.data.frame(t(-log(phyper_net)))
  phper_dual_log = as.data.frame(t(-log(phper_dual)))
  
  saveRDS(phper_each_log , paste0(main.path_tc, num_CancerType, "/",Cancername, "_cellline_pathwayeach_all_log.rds"))
  
  saveRDS(phper_net_log , paste0(main.path_tc, num_CancerType, "/",Cancername, "_cellline_pathwaylink_all_log.rds"))
  
  saveRDS(phper_dual_log , paste0(main.path_tc, num_CancerType, "/",Cancername, "_cellline_dual_all_log.rds"))
  
}

# 
# # (2) - for external data
# 
# main.path_in 
# 
# #   A. for each pathway
# surv.pvals_in = readRDS(paste0(main.path_in,"/surv_pvals_54_KEGG_GC_20220117.rds"))
# phyper_in = readRDS(paste0(main.path_in, "/phyper_enrichment_54_KEGG_GC_20220117.rds"))
# 
# #   B. for pathwaylink
# surv.pvals_in = readRDS(paste0(main.path_in,"/surv_pvals_54_KEGG_net_GC_20220117_filt5.rds"))
# phyper_in = readRDS(paste0(main.path_in, "/phyper_enrichment_54_KEGG_net_GC_20220117.rds"))
# 
# mut.mat_in = readRDS(paste0(main.path_in, CancerType_in,"_mut_count_filt_data.rds"))
# mut.count_in = as.data.frame(colSums(mut.mat_in))
# colnames(mut.count_in) = "counts"
# 
# 
# # (2) - clinical data upload
# 
# cli_ex = fread(paste0(main.path_in,CancerType_in,"_clinical_data.tsv")) 
# cli_ex$`Patient ID`
# cli_ex_surv = cli_ex[,
#                      c("Sample ID",
#                        "DFS months from neo",
#                        "DFS status from neo")]
# 
# cli_ex_surv$vital_status = NA
# cli_ex_surv$vital_status[which(cli_ex_surv$`DFS status from neo` == "1:Recurred/Progressed")] = "Dead"
# cli_ex_surv$vital_status[which(cli_ex_surv$`DFS status from neo` == "0:DiseaseFree")] = "Alive"
# 
# cli_ex_surv$overall_survival = cli_ex_surv$`DFS months from neo`*30
# cli_ex_surv = cli_ex_surv[!is.na(cli_ex_surv$`DFS status from neo`),]
# cli_ex_surv$status = as.numeric(cli_ex_surv$status)
# cli_ex_surv$submitter_id = substr(cli_ex_surv$`Sample ID` , 1,13)
# 
# cli_ex_surv_filt = cli_ex_surv
# while (sum(duplicated(cli_ex_surv_filt$submitter_id)) != 0 ) {
#   print(sum(duplicated(cli_ex_surv_filt$submitter_id)))
#   dupl = which(duplicated(cli_ex_surv_filt$submitter_id))[1]
#   dupl_id = cli_ex_surv_filt[dupl,]$submitter_id
#   tmp.df = cli_ex_surv_filt[which(cli_ex_surv_filt$submitter_id == dupl_id),]
#   tmp.df.2 = tmp.df[-1,]
#   cli_ex_surv_filt = cli_ex_surv_filt[-which(cli_ex_surv_filt$submitter_id == dupl_id),]
#   cli_ex_surv_filt = rbind(cli_ex_surv_filt,tmp.df.2)
#   tmp.df =NULL
#   tmp.df.2 = NULL
# }
# 
# 
# 
# sum(cli_ex_surv_filt$vital_status == "Dead")
# sum(cli_ex_surv_filt$vital_status == "Alive")
# 
# # (2) - filtering only 01A 
# 
# mut.count_in$patient_id = rownames(mut.count_in)
# mut.filt_in = mut.count_in[which(mut.count_in$counts < 1000),]
# 
# int.id_in = intersect(colnames(phyper_in), mut.filt_in$patient_id)
# idx.id_in = match(int.id_in, colnames(phyper_in))
# phyper.filt_in = phyper_in[,idx.id_in]
# 
# # (2) -  transfer data 1 (>0.05) or 100 (<0.05)
# 
# phyper.filt_in_t = t(phyper.filt_in)
# phyper.filt_in_t_num= phyper.filt_in_t[,-which(colSums(phyper.filt_in_t) == dim(phyper.filt_in_t)[1])]
# phyper.filt_in_t_num.bi = ifelse(phyper.filt_in_t_num >0.05 , 1 , 1000000)
# 
# phyper.filt_in_t_bi = t(phyper.filt_in_t_num.bi)
# 
# df.phyper.filt_in_t_bi = as.data.frame(phyper.filt_in_t_bi)
# df.phyper.filt_in_t_bi_num = data.frame(matrix(nrow = length(rownames(df.phyper.filt_in_t_bi))))
# # colnames(df.phyper.filt_in_t_bi) = substr(gsub(pattern = "_", replacement = "-", x = colnames(df.phyper.filt_in_t_bi)), 1,9) 
# colnames(df.phyper.filt_in_t_bi) = gsub(pattern = "_", replacement = "-", x = colnames(df.phyper.filt_in_t_bi))
# for (all_sub in cli_ex_surv_filt$submitter_id) {
#   if (all_sub %in% colnames(df.phyper.filt_in_t_bi)) {
#     # df.phyper.all.link_num = cbind(df.phyper.all.link_num,df.phyper.filt_num[,all_sub])
#     tmp.phyer.filt_num = as.data.frame(df.phyper.filt_in_t_bi[,all_sub])
#     colnames(tmp.phyer.filt_num) = all_sub
#     rownames(tmp.phyer.filt_num) = rownames(df.phyper.filt_in_t_bi)
#     # print(tmp.phyer.filt_num)
#     df.phyper.filt_in_t_bi_num = cbind(df.phyper.filt_in_t_bi_num,tmp.phyer.filt_num)
#     tmp.phyer.filt_num = NULL
#   }
#   
# } 
# df.phyper.filt_in_t_bi_num
# df.phyper.filt_in_t_bi_num= df.phyper.filt_in_t_bi_num[,-1]
# 
# df.phyper.filt_in_t_bi_num = t(df.phyper.filt_in_t_bi_num)
# 
# df.phyper.filt_in_t_bi_num = as.data.frame(df.phyper.filt_in_t_bi_num)
# df.phyper.filt_in_t_bi_num$submitter_id = rownames(df.phyper.filt_in_t_bi_num)
# df.phyper.filt_in_t_bi_num
# 
# 
# df.phyper.merge_in_bi_all = as.data.frame(merge(cli_ex_surv_filt , df.phyper.filt_in_t_bi_num, by = "submitter_id"))
# rownames(df.phyper.merge_in_bi_all) = df.phyper.merge_in_bi_all$submitter_id
# 
# df.phyper.merge_in_bi_all
# 
# df.phyper.merge_bi_in3 = df.phyper.merge_in_bi_all[,c(-1,-2,-3,-4,-6,-7)]
# df.phyper.merge_bi_in3
# df.phyper.merge_vital_bi_in2= df.phyper.merge_bi_in3
# 
# colnames(df.phyper.merge_vital_bi_in2)= gsub("[[:punct:]]", "", colnames(df.phyper.merge_vital_bi_in2))
# df.phyper.merge_vital_bi_in2= df.phyper.merge_vital_bi_in2 %>% relocate("vitalstatus", .after = "P53P54")
# # df.phyper.merge_vital_log_all2= df.phyper.merge_vital_log_all2 %>% relocate("vitalstatus", .after = "P54")
# 
# saveRDS(df.phyper.merge_vital_bi_in2, file = paste0(main.path_in,CancerType_in,"_pathwaylink_all_1or1000000.rds"))
