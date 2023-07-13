library(data.table) 
library(scran)
filepath = "/home/seokwon/nas/"
ref_path = paste0(filepath, "99.reference/")
# Cancerlist = dir(paste0(filepath, "/00.data/total/"))
Cancerlist = dir(paste0(filepath, "/00.data/filtered_TCGA/"))
# Cancerlist = c("34.TCGA-COADREAD" ,"35.TCGA-KIDNEY"  )
sce_path = "/mnt/gluster_server/data/raw/TCGA_data/00.data/"
cli = fread(paste0(ref_path , "all_clin_indexed.csv"))
# num_CancerType = "01.TCGA-READ"
for (num_CancerType in Cancerlist) {
  
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  
  cancer_exp = readRDS(file = paste0(main.path_tc,"/",CancerType,"_exp_TPM_mat_filt_log_gene_symbol.rds"))
  
  cancer_exp = as.data.frame(t(cancer_exp))
  
  if (length(grep("^.*-01A-", rownames(cancer_exp), value = TRUE)) == 0) {
    next
  } else {
    rows_to_keep <- grep("^.*-01A-", rownames(cancer_exp), value = TRUE)
  }

  cancer_exp_filt <- cancer_exp[rows_to_keep, ]
  
  if (sum(duplicated(substr(rownames(cancer_exp_filt),1,12))) != 0) {
    
    for (dup_pat in rownames(cancer_exp_filt[duplicated(substr(rownames(cancer_exp_filt),1,12)),])) {
      tmp_df = cancer_exp_filt[grep(paste0(substr(dup_pat,1,12),"-*"), rownames(cancer_exp_filt)),]
      first = rownames(tmp_df)[1]
      second = rownames(tmp_df)[2]
      tmp_df_mean = as.data.frame(t(colMeans(tmp_df)))
      rownames(tmp_df_mean) = substr(first,1,12)
      
      cancer_exp_filt <- subset(cancer_exp_filt, rownames(cancer_exp_filt) != first)
      cancer_exp_filt <- subset(cancer_exp_filt, rownames(cancer_exp_filt) != second)
      
      cancer_exp_filt = rbind(cancer_exp_filt, tmp_df_mean)
      remove(tmp_df, tmp_df_mean)
    }
  }
  
  dup_names = grep("\\.1$",rownames(cancer_exp_filt),value = TRUE)
  
  if (CancerType %in% c("TCGA-COADREAD","TCGA-KIDNEY")) {
    dup_names_ori = substr(dup_names, 1, 28)
  }else {
    dup_names_ori = substr(dup_names, 1, 16)
  }
  
  
  if (length(dup_names) !=0) {
    for (dup_patients in 1:length(dup_names)) {
      first = dup_names[dup_patients]
      second = dup_names_ori[dup_patients]
      tmp_dup = as.data.frame(cancer_exp_filt[,which(names(cancer_exp_filt) == first)])
      tmp_dup2 = data.frame(lapply(tmp_dup, as.numeric))
      rownames(tmp_dup2) = rownames(cancer_exp_filt)
      colnames(tmp_dup2) = "first"
      
      tmp_dup3 = as.data.frame(cancer_exp_filt[,which(names(cancer_exp_filt) == second)])
      tmp_dup4 = data.frame(lapply(tmp_dup3, as.numeric))
      rownames(tmp_dup4) = rownames(cancer_exp_filt)
      colnames(tmp_dup4) = "second"
      averages <- as.data.frame(rowMeans(cbind(tmp_dup2, tmp_dup4)))
      colnames(averages) = second
      cancer_exp_filt[,first] =NULL
      cancer_exp_filt[,second] = NULL
      cancer_exp_filt = cbind(cancer_exp_filt,averages)
      
    }
  } 
  
  if (sum(duplicated(rownames(cancer_exp_filt))) == 0 ) {
    rownames(cancer_exp_filt) = substr(rownames(cancer_exp_filt), 1, 12)
  }
  
  
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
  cli_surv = cli_surv[!(cli_surv$vital_status == "Not Reported"),]
  
  cli_surv$status = NA
  cli_surv$status[which(cli_surv$vital_status == "Dead")] = 1
  cli_surv$status[which(cli_surv$vital_status == "Alive")] = 0
  
  cli_surv$status = as.numeric(cli_surv$status)

  cancer_exp_filt = as.data.frame(t(cancer_exp_filt))
  cancer_exp_num = data.frame(matrix(nrow = length(rownames(cancer_exp_filt))))
  
  for (all_sub in cli_surv$submitter_id) {
    if (all_sub %in% colnames(cancer_exp_filt)) {
      # df.phyper.all.link_num = cbind(df.phyper.all.link_num,df.phyper.filt_num[,all_sub])
      tmp.phyer.filt_num = as.data.frame(cancer_exp_filt[,all_sub])
      colnames(tmp.phyer.filt_num) = all_sub
      rownames(tmp.phyer.filt_num) = rownames(cancer_exp_filt)
      # print(tmp.phyer.filt_num)
      cancer_exp_num = cbind(cancer_exp_num,tmp.phyer.filt_num)
      tmp.phyer.filt_num = NULL
      }
    }
    
  
  cancer_exp_num= cancer_exp_num[,-1]
  cancer_exp_num_tmp = as.data.frame(t(cancer_exp_num))
  cancer_exp_num_tmp$submitter_id = rownames(cancer_exp_num_tmp)
  
  cancer_exp_num_filt = as.data.frame(merge(cli_surv , cancer_exp_num_tmp, by = "submitter_id"))
  
  rownames(cancer_exp_num_filt) = cancer_exp_num_filt$submitter_id
    
  cancer_exp_num_filt = cancer_exp_num_filt[,c(-1,-2,-3,-4,-5)]

  # hvgs 
  
  sce_for_hvgs = readRDS(paste0(sce_path,num_CancerType,"/", CancerType,".sce_raw_snv_cnv.rds"))

  sce_for_hvgs = logNormCounts(sce_for_hvgs)
  assays_to_remove <- c("counts", "snv", "cnv")
  for (assay_name in assays_to_remove) {
    assay(sce_for_hvgs, assay_name) <- NULL
  }
  dec.sce <- modelGeneVar(sce_for_hvgs)
  
  hvgs = getTopHVGs(dec.sce, n=nrow(cancer_exp_num_filt) * 3)
  
  filt_genes = hvgs[which(hvgs %in% colnames(cancer_exp_num_filt))[1:nrow(cancer_exp_num_filt)]]
  
  cancer_exp_hvgs = cancer_exp_num_filt[,filt_genes]
  cancer_exp_hvgs_fin = cbind(cancer_exp_hvgs, cancer_exp_num_filt[,c("overall_survival", "status")])
  
  saveRDS(cancer_exp_hvgs_fin, paste0("~/nas/hvgs/",CancerType,"_hvgs_exp.rds"))
  write.csv(cancer_exp_hvgs_fin, paste0("~/nas/hvgs/",CancerType,"_hvgs_exp.csv"))
                  
  # random 1
  cancer_exp_num_decreasing = cancer_exp_num_filt[,order(colSums(cancer_exp_num_filt), decreasing = T)][,1:nrow(cancer_exp_num_filt) * 6]
  
  set.seed(13524)
  cancer_exp_rannum1_fin = cancer_exp_num_decreasing[,sample(x=1:ncol(cancer_exp_num_decreasing) ,size = nrow(cancer_exp_num_filt), replace=F)]
  cancer_exp_rannum1_fin = cbind(cancer_exp_rannum1_fin, cancer_exp_num_filt[,c("overall_survival", "status")])
  
  saveRDS(cancer_exp_rannum1_fin, paste0("~/nas/randomly_selected_gene_cancer1/",CancerType,"_random_exp.rds"))
  write.csv(cancer_exp_rannum1_fin, paste0("~/nas/randomly_selected_gene_cancer1/",CancerType,"_random_exp.csv"))
  
  # random 2
  
  set.seed(23524)
  cancer_exp_rannum2_fin = cancer_exp_num_decreasing[,sample(x=1:ncol(cancer_exp_num_decreasing) ,size = nrow(cancer_exp_num_filt), replace=F)]
  cancer_exp_rannum2_fin = cbind(cancer_exp_rannum2_fin, cancer_exp_num_filt[,c("overall_survival", "status")])
  
  saveRDS(cancer_exp_rannum2_fin, paste0("~/nas/randomly_selected_gene_cancer2/",CancerType,"_random_exp.rds"))
  write.csv(cancer_exp_rannum2_fin, paste0("~/nas/randomly_selected_gene_cancer2/",CancerType,"_random_exp.csv"))
  
}
