library(org.Hs.eg.db)
library(data.table)
library(scater)
library(scran)
library(pheatmap)
library(devtools)
library(data.table)
library(scuttle)
library(scran)
install_github("Sage-Bionetworks/CMSclassifier")
library(CMSclassifier)
# zscore normalization

tcga.calc.zscore = function(sce, target.genes){
  message("Calculating z-score with respect to all diploid cells. Version 2023.05.03")
  common.genes = intersect(rownames(sce), target.genes)
  if (length(common.genes) == 0) {
    stop("None of the target genes found in sce. Please check your nomenclature")
  } else if (length(common.genes) != length(target.genes)) {
    message("Some of the genes from query does not exist in this cancer type. It will result in NAs")
    message("Missing genes are: ", paste(setdiff(target.genes, common.genes), collapse = ", "))
  }
  sce.sub = subset(sce, rownames(sce) %in% common.genes,)
  #i am not checking assay names.
  count.mat = assay(sce.sub, 1)
  cnv.mat = assay(sce.sub, 3)
  z.mat = matrix(data = NA, nrow = nrow(count.mat), ncol = ncol(count.mat))
  colnames(z.mat) = colnames(count.mat)
  rownames(z.mat) = rownames(count.mat)
  for (i in 1:nrow(count.mat)) {
    idx.di = which(cnv.mat[i,] == 2)
    query.mean = mean(count.mat[i, idx.di], na.rm = T)
    query.sd = sd(count.mat[i, idx.di], na.rm = T)
    z.mat[i,] = (count.mat[i,] - query.mean)/query.sd
  }
  return(z.mat)
}
filepath = "/home/seokwon/nas/"
ref_path = paste0(filepath, "99.reference/")

Cancerlist = dir(paste0(filepath, "/00.data/filtered_TCGA/"))
sce_path = "/mnt/gluster_server/data/raw/TCGA_data/00.data/"
cli = fread("~/nas/99.reference/all_clin_indexed.csv")
setwd("~/nas/04.Results/compare_surv/")

# for all
for (num_CancerType in Cancerlist) {
  
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  
  # call input
  
  sce = readRDS(paste0(sce_path,num_CancerType,"/", CancerType,".sce_raw_snv_cnv.rds"))
  sce_for_hvgs = sce 
  
  sce_for_hvgs = logNormCounts(sce_for_hvgs)
  # sce_for_hvgs = subset(sce_for_hvgs, rownames(sce_for_hvgs) %in% V(g.ppi.conn)$name,)
  
  assays_to_remove <- c("counts", "snv", "cnv")
  for (assay_name in assays_to_remove) {
    assay(sce_for_hvgs, assay_name) <- NULL
  }
  
  
  dec.sce <- modelGeneVar(sce_for_hvgs)
  hvgs = getTopHVGs(dec.sce, n=1000)
  sce_exp = tcga.calc.zscore(sce = sce, hvgs)
  #loading RDPN output and function 2 output 
  sce_exp = as.data.frame(sce_exp)
  sce_exp = sce_exp[complete.cases(sce_exp),]
  
  sce_count = assay(sce, 1)

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
  
  sce_exp_filt_t_df = as.data.frame(sce_exp)
  
  if (CancerType %in% c("TCGA-COADREAD","TCGA-KIDNEY")) {
    cols_to_keep <- grep("^.*-01A-", colnames(sce_exp_filt_t_df), value = TRUE)
  } else {
    cols_to_keep <- grep("^.*-01A$", colnames(sce_exp_filt_t_df), value = TRUE)
  }
  
  sce_exp_filt_t_df_filt <- sce_exp_filt_t_df[, cols_to_keep]
  
  for (dup_pat in rownames(sce_exp_filt_t_df_filt[duplicated(substr(rownames(sce_exp_filt_t_df_filt),1,12)),])) {
    tmp_df = sce_exp_filt_t_df_filt[grep(paste0(substr(dup_pat,1,12),"-*"), rownames(sce_exp_filt_t_df_filt)),]
    first = rownames(tmp_df)[1]
    second = rownames(tmp_df)[2]
    tmp_df_mean = as.data.frame(t(colMeans(tmp_df)))
    rownames(tmp_df_mean) = substr(first,1,12)
    
    sce_exp_filt_t_df_filt <- subset(sce_exp_filt_t_df_filt, rownames(sce_exp_filt_t_df_filt) != first)
    sce_exp_filt_t_df_filt <- subset(sce_exp_filt_t_df_filt, rownames(sce_exp_filt_t_df_filt) != second)
    
    sce_exp_filt_t_df_filt = rbind(sce_exp_filt_t_df_filt, tmp_df_mean)
    remove(tmp_df, tmp_df_mean)
  }
  
  
  dup_names = grep("\\.1$",names(sce_exp_filt_t_df_filt),value = TRUE)
  if (CancerType %in% c("TCGA-COADREAD","TCGA-KIDNEY")) {
    dup_names_ori = substr(dup_names, 1, 28)
  }else {
    dup_names_ori = substr(dup_names, 1, 16)
  }
  
  
  if (length(dup_names) !=0) {
    for (dup_patients in 1:length(dup_names)) {
      first = dup_names[dup_patients]
      second = dup_names_ori[dup_patients]
      tmp_dup = as.data.frame(sce_exp_filt_t_df_filt[,which(names(sce_exp_filt_t_df_filt) == first)])
      tmp_dup2 = data.frame(lapply(tmp_dup, as.numeric))
      rownames(tmp_dup2) = rownames(sce_exp_filt_t_df_filt)
      colnames(tmp_dup2) = "first"
      
      tmp_dup3 = as.data.frame(sce_exp_filt_t_df_filt[,which(names(sce_exp_filt_t_df_filt) == second)])
      tmp_dup4 = data.frame(lapply(tmp_dup3, as.numeric))
      rownames(tmp_dup4) = rownames(sce_exp_filt_t_df_filt)
      colnames(tmp_dup4) = "second"
      averages <- as.data.frame(rowMeans(cbind(tmp_dup2, tmp_dup4)))
      colnames(averages) = second
      sce_exp_filt_t_df_filt[,first] =NULL
      sce_exp_filt_t_df_filt[,second] = NULL
      sce_exp_filt_t_df_filt = cbind(sce_exp_filt_t_df_filt,averages)
      
    }
  } 
  
  if (sum(duplicated(names(sce_exp_filt_t_df_filt))) == 0 ) {
    colnames(sce_exp_filt_t_df_filt) = substr(colnames(sce_exp_filt_t_df_filt), 1, 12)
  }

  sce_exp_filt2_df = sce_exp_filt_t_df_filt[,substr(colnames(sce_exp_filt_t_df_filt),1,12) %in% cli_surv$submitter_id]
  
  if (ncol(sce_exp_filt_t_df_filt) != length(cli_surv$submitter_id) ) {
    cli_surv_filt  = cli_surv[cli_surv$submitter_id %in% substr(colnames(sce_exp_filt2_df),1,12) ,]
  }
  
  if (CancerType %in% c("TCGA-COADREAD","TCGA-KIDNEY")) {
    
    sce_exp_filt2_df = sce_exp_filt_t_df_filt[,intersect(substr(colnames(sce_exp_filt_t_df_filt),1,12) , cli_surv$submitter_id)]
    
  } 
  
  if (all.equal(ncol(sce_exp_filt2_df),  length(cli_surv_filt$submitter_id) )) {
      sce_exp_filt2_df = as.data.frame(t(sce_exp_filt2_df))
      sce_exp_filt2_df$vital_status = cli_surv_filt$vital_status
      
      }
    
  
  
  
  sce_count_log_with_sym_replaced <- sce_exp_filt2_df
  sce_count_log_with_sym_replaced[is.na(sce_exp_filt2_df) ] <- 0
  
  tmp_vital = sce_exp_filt2_df$vital_status
  sce_exp_filt3_df = sce_exp_filt2_df[,-ncol(sce_exp_filt2_df)]
  sce_exp_filt3_df[sce_exp_filt3_df > 4] <- 4
  sce_exp_filt3_df[sce_exp_filt3_df < -4] <- -4
  sce_exp_filt4_df = cbind(sce_exp_filt3_df,vital_status = tmp_vital)
  
  sce_exp_filt4_df = sce_exp_filt4_df[order(sce_exp_filt4_df$vital_status),]
  annotation_df <- data.frame(vital_status = sce_exp_filt4_df$vital_status)
  rownames(annotation_df) <- rownames(sce_exp_filt4_df)
  
  # Create a named color vector for the unique values of vital_status
  vital_status_colors <- c("Alive" = "yellow", "Dead" = "black")
  names(vital_status_colors) <- unique(annotation_df$vital_status)
  
  png(filename = paste0(CancerType,"_exp_colum_clusterT_surv.png"),
      width = 25, height = 25,  units = "cm" ,pointsize = 12,
      bg = "white", res = 1200, family = "")
  
  tmp = pheatmap::pheatmap(as.matrix(t(sce_exp_filt4_df[,-ncol(sce_exp_filt4_df)])),
                            annotation_col = annotation_df,
                            annotation_colors = list(vital_status = vital_status_colors),
                            cluster_cols = T)
  print(tmp)
  
  dev.off()
  
  png(filename = paste0(CancerType,"_exp_column_clusterF_surv.png"),
      width = 25, height = 25,  units = "cm" ,pointsize = 12,
      bg = "white", res = 1200, family = "")
  
  tmp = pheatmap::pheatmap(as.matrix(t(sce_exp_filt4_df[,-ncol(sce_exp_filt4_df)])),
                           annotation_col = annotation_df,
                           annotation_colors = list(vital_status = vital_status_colors),
                           cluster_cols = F)
  print(tmp)
  
  dev.off()
  
  png(filename = paste0(CancerType,"_exp_cluster_surv_complex.png"),
      width = 25, height = 25,  units = "cm" ,pointsize = 12,
      bg = "white", res = 1200, family = "")
  tmp2 = ComplexHeatmap::pheatmap(as.matrix(t(sce_exp_filt4_df[,-ncol(sce_exp_filt4_df)])),
                                 column_split = annotation_df$vital_status,
                                 annotation_col = annotation_df,
                                 annotation_colors = list(vital_status = vital_status_colors),
                                 cluster_cols = T)

  print(tmp2)

  dev.off()
  
  remove(tmp,tmp2,sce_exp_filt4_df,sce_exp_filt3_df,sce_exp_filt2_df,cli_surv_filt,sce_exp_filt_t_df_filt,cli_surv, sce_for_hvgs, sce)
}

