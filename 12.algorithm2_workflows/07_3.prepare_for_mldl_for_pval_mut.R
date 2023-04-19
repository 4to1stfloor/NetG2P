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

filepath = paste0("/home/seokwon/nas/")
ref_path = paste0(filepath, "/99.reference/")
type = "each"
bi_num_mode = "log"

Cancerlist = dir(paste0(filepath, "/00.data/filtered_TCGA/"))

for (num_CancerType in Cancerlist) {
  # 1) prepare for ml
  # (1) for TCGA data upload
  
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  
  if (CancerType == "Total-TCGA") {
    break
  }
  
  #   A. for each pathway
  if (type == "each") {
    surv.pvals = readRDS(paste0(main.path_tc,"/surv_pvals_54_mut_KEGG_GC_20220117.rds"))
    phyper = readRDS(paste0(main.path_tc, "/phyper_enrichment_mut_54_KEGG_GC_20220117.rds"))
  } else if (type == "link") {
    #   B. for pathwaylink
    surv.pvals = readRDS(paste0(main.path_tc,"/surv_pvals_54_mut_KEGG_net_GC_20220117.rds"))
    phyper = readRDS(paste0(main.path_tc, "/phyper_enrichment_mut_54_KEGG_net_GC_20220117.rds"))
  }
  
  mut.mat = readRDS(paste0(main.path_tc,"/", CancerType,"_mut_count_filt_data.rds"))
  mut.count = as.data.frame(colSums(mut.mat))
  colnames(mut.count) = "counts"
  
  # (1) - clinical data upload
  
  cli = fread(paste0(ref_path , "all_clin_indexed.csv"))
  if (CancerType == "TCGA-COADREAD") {
    cli_surv = cli[cli$project %in% c("TCGA-COAD","TCGA-READ"),
                   c("submitter_id",
                     "vital_status",
                     "days_to_death",
                     "days_to_last_follow_up")]
  } else if (CancerType == "TCGA-KIDNEY") {
    cli_surv = cli[cli$project %in% c("TCGA-KIRP","TCGA-KIRC","TCGA-KICH"),
                   c("submitter_id",
                     "vital_status",
                     "days_to_death",
                     "days_to_last_follow_up")]
  } else {
    cli_surv = cli[cli$project == CancerType,
                   c("submitter_id",
                     "vital_status",
                     "days_to_death",
                     "days_to_last_follow_up")]
  }
  
  
  cli_surv$deceased = cli_surv$vital_status == "Dead"
  cli_surv$overall_survival = ifelse(cli_surv$deceased,
                                     cli_surv$days_to_death,
                                     cli_surv$days_to_last_follow_up)
  
  cli_surv = cli_surv[!is.na(cli_surv$vital_status),]
  
  cli_surv$status = NA
  cli_surv$status[which(cli_surv$vital_status == "Dead")] = 1
  cli_surv$status[which(cli_surv$vital_status == "Alive")] = 0
  
  cli_surv$status = as.numeric(cli_surv$status)
  
  # (1) - filtering only 01A

  if (CancerType %in% c("TCGA-COADREAD","TCGA-KIDNEY")) {
    mut.count$patient_id = substr(rownames(mut.count), 1, 12) 
  } else {
    mut.count$patient_id = substr(rownames(mut.count), 1, 16) 
  }
  
  mut.count$forfilt =  rownames(mut.count)
  mut.count.filt = mut.count %>% filter(str_detect(forfilt, "-01A-")) 
  mut.filt = mut.count.filt[which(mut.count.filt$counts < 1000),]
  mut.filt$forfilt = NULL
  
  if (CancerType %in% c("TCGA-COADREAD","TCGA-KIDNEY")) {
    int.id = intersect(substr(colnames(phyper),1,12), mut.filt$patient_id) 
    idx.id = match(int.id, substr(colnames(phyper),1,12))
  } else {
    int.id = intersect(substr(colnames(phyper), 1, 16), mut.filt$patient_id) 
    idx.id = match(int.id, substr(colnames(phyper), 1, 16))
  }
  if (length(idx.id) == 0) {
    print("id does not matched")
  }
  
  phyper_filt = phyper[,idx.id]
  
  #  transfer data 1 (>0.05) or 100 (<0.05)
  
  if (length(which(colSums(phyper_filt) == dim(phyper_filt)[1])) !=0) {
    phyper.filt_num= phyper_filt[,-which(colSums(phyper_filt) == dim(phyper_filt)[1])]
  }else {
    phyper.filt_num = phyper_filt
  }
  phyper.filt_t_num = t(phyper.filt_num)

  if (bi_num_mode == "log") {
    phyper.filt_t_num.log = -log(phyper.filt_t_num)
  } else if (bi_num_mode == "bi") {
    phyper.filt_t_num.log = ifelse(phyper.filt_t_num < 0.05, 1 , 0)
  } 
  
  rownames(phyper.filt_t_num.log) = substr(rownames(phyper.filt_t_num.log), 1,12)
  
  df.phyper.filt_t_log_num = as.data.frame(phyper.filt_t_num.log)
  df.phyper.filt_t_log_num$submitter_id = rownames(df.phyper.filt_t_log_num)
  df.phyper.filt_t_log_num
  
  df.phyper.merge_log_all = as.data.frame(merge(cli_surv , df.phyper.filt_t_log_num, by = "submitter_id"))
  rownames(df.phyper.merge_log_all) = df.phyper.merge_log_all$submitter_id
  
  df.phyper.merge_log_all
  
  df.phyper.merge_log_all3 = df.phyper.merge_log_all[,c(-1,-3,-4,-5,-6,-7)]
  df.phyper.merge_log_all3
  df.phyper.merge_vital_log_all2= df.phyper.merge_log_all3
  
  colnames(df.phyper.merge_vital_log_all2)= gsub("[[:punct:]]", "", colnames(df.phyper.merge_vital_log_all2))
  if (type == "link") {
    df.phyper.merge_vital_log_all2= df.phyper.merge_vital_log_all2 %>% relocate("vitalstatus", .after = "P53P54")
    saveRDS(df.phyper.merge_vital_log_all2, file = paste0(main.path_tc,"/",CancerType,"_mut_pathwaylink_all_log.rds"))
    
  } else if (type == "each") {
    df.phyper.merge_vital_log_all2= df.phyper.merge_vital_log_all2 %>% relocate("vitalstatus", .after = "P54")
    saveRDS(df.phyper.merge_vital_log_all2, file = paste0(main.path_tc,"/",CancerType,"_mut_pathwayeach_all_log.rds"))
  }
  
}
