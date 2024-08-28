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
all_clinical = read.csv(file = paste0(ref_path,"all_clin_indexed.csv"))

Cancerlist = Cancerlist[c(-2,-3,-4,-7,-17)]
type = "each"
mode = "dual"
bi_num_mode = "log"

num_CancerType = "04.TCGA-CESC"
for (num_CancerType in Cancerlist) {

  main.path_tc = paste0(filepath, "00.data/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  
  if (mode == "dual" && bi_num_mode == "1") {
    data_tc_link = readRDS(file = paste0(main.path_tc,"/",CancerType,"_pathwaylink_all_0or1.rds"))
    data_tc_each = readRDS(file = paste0(main.path_tc,"/",CancerType,"_pathwayeach_all_0or1.rds"))
    
    if (all.equal(rownames(data_tc_link), rownames(data_tc_each))) {
      data_tc_link$vitalstatus = NULL
      data_tc = cbind(data_tc_link,data_tc_each)
    }
  } else if (mode == "dual" && bi_num_mode == "log") {
    data_tc_link = readRDS(file = paste0(main.path_tc,"/",CancerType,"_pathwaylink_all_log.rds"))
    data_tc_each = readRDS(file = paste0(main.path_tc,"/",CancerType,"_pathwayeach_all_log.rds"))
    
    if (all.equal(rownames(data_tc_link), rownames(data_tc_each))) {
      data_tc_link$vitalstatus = NULL
      data_tc = cbind(data_tc_link,data_tc_each)
    }
    
  } else {
    
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
    } else if (type == "link" && bi_num_mode == "log") {
      data_tc = readRDS(file = paste0(main.path_tc,"/",CancerType,"_pathwaylink_all_log.rds"))
    } else if (type == "each" && bi_num_mode == "log") {
      data_tc = readRDS(file = paste0(main.path_tc,"/",CancerType,"_pathwayeach_all_log.rds"))
    }
    
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
  
  cli = fread(paste0(ref_path , "all_clin_indexed.csv"))
  if (CancerType == "TCGA-COADREAD") {
    
    cli_coad = cli[cli$project == "TCGA-COAD",
                   c("submitter_id",
                     "vital_status",
                     "days_to_death",
                     "days_to_last_follow_up")]
    
    cli_read = cli[cli$project == "TCGA-READ",
                   c("submitter_id",
                     "vital_status",
                     "days_to_death",
                     "days_to_last_follow_up")]
    cli_surv = rbind(cli_coad,cli_read)
    
  } else if (CancerType == "TCGA-KIDNEY") {
    
    cli_kich = cli[cli$project == "TCGA-KICH",
                   c("submitter_id",
                     "vital_status",
                     "days_to_death",
                     "days_to_last_follow_up")]
    
    cli_kirc = cli[cli$project == "TCGA-KIRC",
                   c("submitter_id",
                     "vital_status",
                     "days_to_death",
                     "days_to_last_follow_up")]
    cli_kirp = cli[cli$project == "TCGA-KIRP",
                   c("submitter_id",
                     "vital_status",
                     "days_to_death",
                     "days_to_last_follow_up")]
    cli_surv = rbind(cli_kich,cli_kirc,cli_kirp)
    
  }  else {
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
  cli_surv
  
  data_tc_ori$duration = NA
  for (pat_id in rownames(data_tc_ori)) {
    data_tc_ori[pat_id,]$duration = cli_surv[which(cli_surv$submitter_id == pat_id),]$overall_survival
    data_tc_ori[pat_id,]$status = cli_surv[which(cli_surv$submitter_id == pat_id),]$status
  } 
  data_tc_ori$vitalstatus = NULL

  write.csv(data_tc_ori, file = paste0(main.path_tc, "/", CancerType, "_dual_add_duration_log.csv"))

  
  # fill_cli = read.csv(file = paste0(ref_path,CancerType,"_fill.csv"))
  # fill_cli = fill_cli[,c("submitter_id",
  #                        "vital_status",
  #                        "days_to_death",
  #                        "days_to_last_follow_up",
  #                        "ajcc_pathologic_stage",
  #                        "icd_10_code")]
  # 
  # fill_cli$deceased = fill_cli$vital_status == "Dead"
  # 
  # fill_cli$overall_survival = ifelse(fill_cli$deceased,
  #                                    fill_cli$days_to_death,
  #                                    fill_cli$days_to_last_follow_up)
  # 
  # fill_cli = fill_cli[!is.na(fill_cli$vital_status),]
  # 
  # fill_cli$status = NA
  # fill_cli$status[which(fill_cli$vital_status == "Dead")] = 1
  # fill_cli$status[which(fill_cli$vital_status == "Alive")] = 0
  # 
  # fill_cli$status = as.numeric(fill_cli$status)
  # 
  # fill_cli_df = as.data.frame(get_dummies.(fill_cli , c("ajcc_pathologic_stage","icd_10_code")))
  

  
}



### fin