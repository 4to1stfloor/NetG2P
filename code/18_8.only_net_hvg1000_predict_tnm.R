library(devtools)
library(CMSclassifier)
library(org.Hs.eg.db)
library(data.table)
library(scater)
library(scran)
library(pheatmap)
library(maftools)
library(RColorBrewer)

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
setwd("~/nas/04.Results/subtypes/tmn/")

# for all
for (num_CancerType in Cancerlist) {
  
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  
  if (CancerType == "TCGA-COADREAD") {
    cli_surv = cli[cli$project %in% c("TCGA-COAD","TCGA-READ"),
                   c("submitter_id",
                     "vital_status",
                     "days_to_death",
                     "days_to_last_follow_up",
                     "ajcc_pathologic_t",
                     "ajcc_pathologic_n",
                     "ajcc_pathologic_m")]
    
  } else if (CancerType == "TCGA-KIDNEY") {
    cli_surv = cli[cli$project %in% c("TCGA-KIRP","TCGA-KIRC","TCGA-KICH"),
                   c("submitter_id",
                     "vital_status",
                     "days_to_death",
                     "days_to_last_follow_up",
                     "ajcc_pathologic_t",
                     "ajcc_pathologic_n",
                     "ajcc_pathologic_m")]
    
  } else {
    cli_surv = cli[cli$project == CancerType,
                   c("submitter_id",
                     "vital_status",
                     "days_to_death",
                     "days_to_last_follow_up",
                     "ajcc_pathologic_t",
                     "ajcc_pathologic_n",
                     "ajcc_pathologic_m")]
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
  
  if (sum(!is.na(cli_surv$ajcc_pathologic_t)) == 0 || sum(!is.na(cli_surv$ajcc_pathologic_n)) == 0 || sum(!is.na(cli_surv$ajcc_pathologic_m)) == 0) {
    print(paste0( "there are not tmn data in ",CancerType))
    next
  }
  
  cli_surv_filtered = cli_surv[which(!is.na(cli_surv$ajcc_pathologic_t)),]
  cli_surv_filtered = cli_surv_filtered[which(!is.na(cli_surv_filtered$ajcc_pathologic_n)),]
  cli_surv_filtered = cli_surv_filtered[which(!is.na(cli_surv_filtered$ajcc_pathologic_m)),]
  
  cli_surv_filtered$ajcc_pathologic_t <- ifelse(grepl("^T1", cli_surv_filtered$ajcc_pathologic_t), "T1", cli_surv_filtered$ajcc_pathologic_t)
  cli_surv_filtered$ajcc_pathologic_t <- ifelse(grepl("^T2", cli_surv_filtered$ajcc_pathologic_t), "T2", cli_surv_filtered$ajcc_pathologic_t)
  cli_surv_filtered$ajcc_pathologic_t <- ifelse(grepl("^T3", cli_surv_filtered$ajcc_pathologic_t), "T3", cli_surv_filtered$ajcc_pathologic_t)
  cli_surv_filtered$ajcc_pathologic_t <- ifelse(grepl("^T4", cli_surv_filtered$ajcc_pathologic_t), "T4", cli_surv_filtered$ajcc_pathologic_t)
  
  cli_surv_filtered$ajcc_pathologic_n <- ifelse(grepl("^N1", cli_surv_filtered$ajcc_pathologic_n), "N1", cli_surv_filtered$ajcc_pathologic_n)
  cli_surv_filtered$ajcc_pathologic_n <- ifelse(grepl("^N2", cli_surv_filtered$ajcc_pathologic_n), "N2", cli_surv_filtered$ajcc_pathologic_n)
  cli_surv_filtered$ajcc_pathologic_n <- ifelse(grepl("^N3", cli_surv_filtered$ajcc_pathologic_n), "N3", cli_surv_filtered$ajcc_pathologic_n)
  
  cli_surv_filtered$ajcc_pathologic_m <- ifelse(grepl("^M1", cli_surv_filtered$ajcc_pathologic_m), "M1", cli_surv_filtered$ajcc_pathologic_m)
  
  # call input
  
  sce = readRDS(paste0(sce_path,num_CancerType,"/", CancerType,".sce_raw_snv_cnv.rds"))
  net = readRDS(paste0(main.path_tc, "/", list.files(main.path_tc, pattern = "^net_prop_total")))
  
  sce_for_hvgs = sce 
  
  sce_for_hvgs = logNormCounts(sce_for_hvgs)
  # sce_for_hvgs = subset(sce_for_hvgs, rownames(sce_for_hvgs) %in% V(g.ppi.conn)$name,)
  
  assays_to_remove <- c("counts", "snv", "cnv")
  for (assay_name in assays_to_remove) {
    assay(sce_for_hvgs, assay_name) <- NULL
  }
  
  dec.sce <- modelGeneVar(sce_for_hvgs)
  hvgs = getTopHVGs(dec.sce, n=1000)
  
  net_df = as.data.frame(net)
  
  cols_to_keep <- grep("^.*-01A-", colnames(net_df), value = TRUE)
  
  net_filt_df <- net_df[, cols_to_keep]
  
  if (sum(duplicated(substr(rownames(net_filt_df),1,12))) != 0) {
    
    for (dup_pat in rownames(net_filt_df[duplicated(substr(rownames(net_filt_df),1,12)),])) {
      tmp_df = net_filt_df[grep(paste0(substr(dup_pat,1,12),"-*"), rownames(net_filt_df)),]
      first = rownames(tmp_df)[1]
      second = rownames(tmp_df)[2]
      tmp_df_mean = as.data.frame(t(colMeans(tmp_df)))
      rownames(tmp_df_mean) = substr(first,1,12)
      
      net_filt_df <- subset(net_filt_df, rownames(net_filt_df) != first)
      net_filt_df <- subset(net_filt_df, rownames(net_filt_df) != second)
      
      net_filt_df = rbind(net_filt_df, tmp_df_mean)
      remove(tmp_df, tmp_df_mean)
    }
  }
  
  dup_names = grep("\\.1$",names(net_filt_df),value = TRUE)
  
  if (CancerType %in% c("TCGA-COADREAD","TCGA-KIDNEY")) {
    dup_names_ori = substr(dup_names, 1, 28)
  }else {
    dup_names_ori = substr(dup_names, 1, 16)
  }
  
  
  if (length(dup_names) !=0) {
    for (dup_patients in 1:length(dup_names)) {
      first = dup_names[dup_patients]
      second = dup_names_ori[dup_patients]
      tmp_dup = as.data.frame(net_filt_df[,which(names(net_filt_df) == first)])
      tmp_dup2 = data.frame(lapply(tmp_dup, as.numeric))
      rownames(tmp_dup2) = rownames(net_filt_df)
      colnames(tmp_dup2) = "first"
      
      tmp_dup3 = as.data.frame(net_filt_df[,which(names(net_filt_df) == second)])
      tmp_dup4 = data.frame(lapply(tmp_dup3, as.numeric))
      rownames(tmp_dup4) = rownames(net_filt_df)
      colnames(tmp_dup4) = "second"
      averages <- as.data.frame(rowMeans(cbind(tmp_dup2, tmp_dup4)))
      colnames(averages) = second
      net_filt_df[,first] =NULL
      net_filt_df[,second] = NULL
      net_filt_df = cbind(net_filt_df,averages)
      
    }
  } 
  
  if (sum(duplicated(names(net_filt_df))) == 0 ) {
    colnames(net_filt_df) = substr(colnames(net_filt_df), 1, 12)
  }
  
  net_filt2_df = net_filt_df[,substr(colnames(net_filt_df),1,12) %in% cli_surv_filtered$submitter_id]
  
  if (ncol(net_filt2_df) != length(cli_surv_filtered$submitter_id) ) {
    tmp_pat = intersect(cli_surv_filtered$submitter_id,substr(colnames(net_filt2_df),1,12))
    cli_surv_filt  = cli_surv_filtered[which(cli_surv_filtered$submitter_id %in% tmp_pat),]
    net_filt2_df  = net_filt2_df[,tmp_pat]
  }
  
  if (CancerType %in% c("TCGA-COADREAD","TCGA-KIDNEY")) {
    
    net_filt2_df = net_filt2_df[,intersect(substr(colnames(net_filt2_df),1,12) , cli_surv_filtered$submitter_id)]
    
  } 
  
  if (all.equal(ncol(net_filt2_df),  length(cli_surv_filt$submitter_id) )) {
    net_filt2_df = as.data.frame(t(net_filt2_df))
    net_filt2_df$pat_t = cli_surv_filt$ajcc_pathologic_t 
    net_filt2_df$pat_n = cli_surv_filt$ajcc_pathologic_n 
    net_filt2_df$pat_m = cli_surv_filt$ajcc_pathologic_m 
  }
  
  net_filt3_df = net_filt2_df[,which(colnames(net_filt2_df) %in% hvgs)]
  # net_filt3_replaced <- net_filt3_df
  # net_filt3_replaced[is.na(net_filt3_replaced) ] <- 0
  
  tmp_tnm = net_filt2_df[,c("pat_t","pat_n","pat_m")]

  net_filt3_df[net_filt3_df > 4] <- 4
  net_filt3_df[net_filt3_df < -4] <- -4
  net_filt4_df_df = cbind(net_filt3_df,tmp_tnm)
  
  net_filt4_df_df = net_filt4_df_df[order(net_filt4_df_df$pat_t),]
  annotation_df <- data.frame(pat_t = net_filt4_df_df$pat_t,
                              pat_n = net_filt4_df_df$pat_n,
                              pat_m = net_filt4_df_df$pat_m)
  rownames(annotation_df) <- rownames(net_filt4_df_df)
  
  # Create a named color vector for the unique values of vital_status
  
  num_pat_t = brewer.pal(length(unique(net_filt4_df_df$pat_t)), "RdGy")
  col_pat_t = setNames(num_pat_t, unique(net_filt4_df_df$pat_t))
  
  num_pat_n = brewer.pal(length(unique(net_filt4_df_df$pat_n)), "Spectral")
  col_pat_n = setNames(num_pat_n, unique(net_filt4_df_df$pat_n))
  
  num_pat_m = brewer.pal(length(unique(net_filt4_df_df$pat_m)), "Set3")
  col_pat_m = setNames(num_pat_m, unique(net_filt4_df_df$pat_m))
  
  tnm_colors = list(pat_t = col_pat_t , pat_n = col_pat_n, pat_m = col_pat_m)
  
  #
  png(filename = paste0(CancerType,"_net_colum_clusterT_tnm.png"),
      width = 25, height = 25,  units = "cm" ,pointsize = 12,
      bg = "white", res = 1200, family = "")
  
  tmp = pheatmap::pheatmap(as.matrix(t(net_filt4_df_df[,which(!colnames(net_filt4_df_df) %in% c("pat_t","pat_n","pat_m"))])),
                           annotation_col = annotation_df,
                           annotation_colors = tnm_colors,
                           cluster_cols = T)
  print(tmp)
  
  dev.off()
  
  png(filename = paste0(CancerType,"_net_column_clusterF_tnm.png"),
      width = 25, height = 25,  units = "cm" ,pointsize = 12,
      bg = "white", res = 1200, family = "")
  
  tmp2 = pheatmap::pheatmap(as.matrix(t(net_filt4_df_df[,which(!colnames(net_filt4_df_df) %in% c("pat_t","pat_n","pat_m"))])),
                            annotation_col = annotation_df,
                            annotation_colors = tnm_colors,
                            cluster_cols = F)
  print(tmp2)
  
  dev.off()
  
  png(filename = paste0(CancerType,"_net_cluster_t_complex.png"),
      width = 25, height = 25,  units = "cm" ,pointsize = 12,
      bg = "white", res = 1200, family = "")
  tmp3 = ComplexHeatmap::pheatmap(as.matrix(t(net_filt4_df_df[,which(!colnames(net_filt4_df_df) %in% c("pat_t","pat_n","pat_m"))])),
                                  column_split = annotation_df$pat_t,
                                  annotation_col = annotation_df,
                                  annotation_colors = tnm_colors,
                                  cluster_cols = T)
  print(tmp3)
  
  dev.off()
  
  png(filename = paste0(CancerType,"_net_cluster_n_complex.png"),
      width = 25, height = 25,  units = "cm" ,pointsize = 12,
      bg = "white", res = 1200, family = "")
  
  tmp4 = ComplexHeatmap::pheatmap(as.matrix(t(net_filt4_df_df[,which(!colnames(net_filt4_df_df) %in% c("pat_t","pat_n","pat_m"))])),
                                  column_split = annotation_df$pat_n,
                                  annotation_col = annotation_df,
                                  annotation_colors = tnm_colors,
                                  cluster_cols = T)
  
  print(tmp4)
  
  dev.off()
  
  png(filename = paste0(CancerType,"_net_cluster_m_complex.png"),
      width = 25, height = 25,  units = "cm" ,pointsize = 12,
      bg = "white", res = 1200, family = "")
  tmp5 = ComplexHeatmap::pheatmap(as.matrix(t(net_filt4_df_df[,which(!colnames(net_filt4_df_df) %in% c("pat_t","pat_n","pat_m"))])),
                                  column_split = annotation_df$pat_m,
                                  annotation_col = annotation_df,
                                  annotation_colors = tnm_colors,
                                  cluster_cols = T)
  
  print(tmp5)
  
  dev.off()
  
}
