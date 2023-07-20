library(readxl)
library(tidygraph)
library(data.table) 
library(scran)
filepath = "/home/seokwon/nas/"
ref_path = paste0(filepath, "99.reference/")
# Cancerlist = dir(paste0(filepath, "/00.data/total/"))
Cancerlist = dir(paste0(filepath, "/00.data/filtered_TCGA/"))
# Cancerlist = c("34.TCGA-COADREAD" ,"35.TCGA-KIDNEY"  )
sce_path = "/mnt/gluster_server/data/raw/TCGA_data/00.data/"

# num_CancerType = "04.TCGA-CESC"
for (num_CancerType in Cancerlist) {
  
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  
  # call input
  
  short_long = readRDS(paste0(filepath,"04.Results/short_long/",CancerType,"_best_features_short_long.rds"))
  
  short_long_filt = short_long %>% select(-vitalstatus, -duration, -status)
  
  write.csv(short_long_filt, paste0(filepath,"04.Results/short_long/",CancerType,"_critical_features_short_long.csv"))
  
  
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

  cancer_exp_filt_sl = cancer_exp_filt[rownames(short_long_filt),]
  
  cancer_exp_num_filt = as.data.frame(t(cancer_exp_filt_sl))
  
  # hvgs 
  
  sce_for_hvgs = readRDS(paste0(sce_path,num_CancerType,"/", CancerType,".sce_raw_snv_cnv.rds"))
  
  sce_for_hvgs = logNormCounts(sce_for_hvgs)
  assays_to_remove <- c("counts", "snv", "cnv")
  for (assay_name in assays_to_remove) {
    assay(sce_for_hvgs, assay_name) <- NULL
  }
  dec.sce <- modelGeneVar(sce_for_hvgs)
  
  hvgs = getTopHVGs(dec.sce, n=nrow(cancer_exp_num_filt) * 3)
  
  filt_genes = hvgs[which(hvgs %in% rownames(cancer_exp_num_filt))[1:ncol(cancer_exp_num_filt)]]
  
  cancer_exp_hvgs = cancer_exp_num_filt[filt_genes,]
  cancer_exp_hvgs_fin = as.data.frame(t(cancer_exp_hvgs))
  
  if (all.equal(rownames(cancer_exp_hvgs_fin) , rownames(short_long_filt))) {
    cancer_exp_hvgs_fin$cluster = short_long_filt$cluster
  }
  
  saveRDS(cancer_exp_hvgs_fin, paste0("~/nas/hvgs/",CancerType,"_short_long_hvgs_exp.rds"))
  write.csv(cancer_exp_hvgs_fin, paste0("~/nas/hvgs/",CancerType,"_short_long_hvgs_exp.csv"))
  
  # random 1

  hvgs_for_ran = getTopHVGs(dec.sce, n=nrow(cancer_exp_num_filt) * 6)
  filt_genes_for_ran = hvgs_for_ran[which(hvgs_for_ran %in% rownames(cancer_exp_num_filt))]
  filt_genes_wo_hvgs = filt_genes_for_ran[!( filt_genes_for_ran %in%  filt_genes)]
  
  tmp_cancer_exp_num_filt = as.data.frame(t(cancer_exp_num_filt[filt_genes_wo_hvgs,]))
  tmp_cancer_exp_num_filt <- lapply(tmp_cancer_exp_num_filt, as.numeric)
  tmp_cancer_exp_num_filt <- data.frame(tmp_cancer_exp_num_filt, check.names = FALSE)
  
  # Assuming tmp_cancer_exp_num_filt is a list, convert it to a dataframe while preserving row names
  tmp_cancer_exp_num_filt <- as.data.frame(t(cancer_exp_num_filt[filt_genes_wo_hvgs,]))
  
  # Convert each column to numeric
  tmp_cancer_exp_num_filt <- lapply(tmp_cancer_exp_num_filt, as.numeric)
  
  # Create a new dataframe with the converted columns and preserve row names
  tmp_cancer_exp_num_filt <- data.frame(tmp_cancer_exp_num_filt, check.names = FALSE)
  
  # Restore the original row names
  rownames(tmp_cancer_exp_num_filt) <- rownames(as.data.frame(t(cancer_exp_num_filt[filt_genes_wo_hvgs,])))

  set.seed(13524)
  cancer_exp_rannum1_fin = tmp_cancer_exp_num_filt[,order(colSums(tmp_cancer_exp_num_filt), decreasing = T)][,sample(x=1:ncol(tmp_cancer_exp_num_filt) ,size = nrow(tmp_cancer_exp_num_filt), replace=F)]
  
  if (all.equal(rownames(cancer_exp_rannum1_fin) , rownames(short_long_filt))) {
    cancer_exp_rannum1_fin$cluster = short_long_filt$cluster
  }
  
  saveRDS(cancer_exp_rannum1_fin, paste0("~/nas/randomly_selected_gene_cancer1/",CancerType,"_short_long_random_exp.rds"))
  write.csv(cancer_exp_rannum1_fin, paste0("~/nas/randomly_selected_gene_cancer1/",CancerType,"_short_long_random_exp.csv"))
  
  # random 2
  
  set.seed(23524)
  cancer_exp_rannum2_fin = tmp_cancer_exp_num_filt[,order(colSums(tmp_cancer_exp_num_filt), decreasing = T)][,sample(x=1:ncol(tmp_cancer_exp_num_filt) ,size = nrow(tmp_cancer_exp_num_filt), replace=F)]
  
  if (all.equal(rownames(cancer_exp_rannum2_fin) , rownames(short_long_filt))) {
    cancer_exp_rannum2_fin$cluster = short_long_filt$cluster
  }
  
  saveRDS(cancer_exp_rannum2_fin, paste0("~/nas/randomly_selected_gene_cancer2/",CancerType,"_short_long_random_exp.rds"))
  write.csv(cancer_exp_rannum2_fin, paste0("~/nas/randomly_selected_gene_cancer2/",CancerType,"_short_long_random_exp.csv"))
  
}
