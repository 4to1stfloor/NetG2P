
library(devtools)
library(CMSclassifier)
library(org.Hs.eg.db)
library(data.table)
library(scater)
library(scran)
library(pheatmap)
library(maftools)
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
  mut = readRDS(paste0(sce_path,num_CancerType,"/", CancerType,"_maf.rds"))
  mut_filtered <- subsetMaf(maf = mut, query = "Variant_Classification != 'Nonsense_Mutation'")
  
  mut_count = mutCountMatrix(mut_filtered)
  sce_for_hvgs = sce 
  
  sce_for_hvgs = logNormCounts(sce_for_hvgs)
  # sce_for_hvgs = subset(sce_for_hvgs, rownames(sce_for_hvgs) %in% V(g.ppi.conn)$name,)
  
  assays_to_remove <- c("counts", "snv", "cnv")
  for (assay_name in assays_to_remove) {
    assay(sce_for_hvgs, assay_name) <- NULL
  }
  
  dec.sce <- modelGeneVar(sce_for_hvgs)
  hvgs = getTopHVGs(dec.sce, n=1000)
  
 
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
  
  mut_count = as.data.frame(mut_count)
  
  # if (CancerType %in% c("TCGA-COADREAD","TCGA-KIDNEY")) {
  #   cols_to_keep <- grep("^.*-01A-", colnames(mut_count), value = TRUE)
  # } else {
  #   cols_to_keep <- grep("^.*-01A$", colnames(mut_count), value = TRUE)
  # }
  cols_to_keep <- grep("^.*-01A-", colnames(mut_count), value = TRUE)
 
  mut_count_filt_df <- mut_count[, cols_to_keep]
  
  if (sum(duplicated(substr(rownames(mut_count_filt_df),1,12))) != 0) {
    
    for (dup_pat in rownames(mut_count_filt_df[duplicated(substr(rownames(mut_count_filt_df),1,12)),])) {
      tmp_df = mut_count_filt_df[grep(paste0(substr(dup_pat,1,12),"-*"), rownames(mut_count_filt_df)),]
      first = rownames(tmp_df)[1]
      second = rownames(tmp_df)[2]
      tmp_df_mean = as.data.frame(t(colMeans(tmp_df)))
      rownames(tmp_df_mean) = substr(first,1,12)
      
      mut_count_filt_df <- subset(mut_count_filt_df, rownames(mut_count_filt_df) != first)
      mut_count_filt_df <- subset(mut_count_filt_df, rownames(mut_count_filt_df) != second)
      
      mut_count_filt_df = rbind(mut_count_filt_df, tmp_df_mean)
      remove(tmp_df, tmp_df_mean)
    }
  }
  
  dup_names = grep("\\.1$",names(mut_count_filt_df),value = TRUE)
  
  if (CancerType %in% c("TCGA-COADREAD","TCGA-KIDNEY")) {
    dup_names_ori = substr(dup_names, 1, 28)
  }else {
    dup_names_ori = substr(dup_names, 1, 16)
  }
  
  
  if (length(dup_names) !=0) {
    for (dup_patients in 1:length(dup_names)) {
      first = dup_names[dup_patients]
      second = dup_names_ori[dup_patients]
      tmp_dup = as.data.frame(mut_count_filt_df[,which(names(mut_count_filt_df) == first)])
      tmp_dup2 = data.frame(lapply(tmp_dup, as.numeric))
      rownames(tmp_dup2) = rownames(mut_count_filt_df)
      colnames(tmp_dup2) = "first"
      
      tmp_dup3 = as.data.frame(mut_count_filt_df[,which(names(mut_count_filt_df) == second)])
      tmp_dup4 = data.frame(lapply(tmp_dup3, as.numeric))
      rownames(tmp_dup4) = rownames(mut_count_filt_df)
      colnames(tmp_dup4) = "second"
      averages <- as.data.frame(rowMeans(cbind(tmp_dup2, tmp_dup4)))
      colnames(averages) = second
      mut_count_filt_df[,first] =NULL
      mut_count_filt_df[,second] = NULL
      mut_count_filt_df = cbind(mut_count_filt_df,averages)
      
    }
  } 
  
  if (sum(duplicated(names(mut_count_filt_df))) == 0 ) {
    colnames(mut_count_filt_df) = substr(colnames(mut_count_filt_df), 1, 12)
  }

  mut_count_filt2_df = mut_count_filt_df[,substr(colnames(mut_count_filt_df),1,12) %in% cli_surv$submitter_id]
  
  if (ncol(mut_count_filt2_df) != length(cli_surv$submitter_id) ) {
    tmp_pat = intersect(cli_surv$submitter_id,substr(colnames(mut_count_filt2_df),1,12))
    cli_surv_filt  = cli_surv[which(cli_surv$submitter_id %in% tmp_pat),]
    mut_count_filt2_df  = mut_count_filt2_df[,tmp_pat]
  }
  
  if (CancerType %in% c("TCGA-COADREAD","TCGA-KIDNEY")) {
    
    mut_count_filt2_df = mut_count_filt2_df[,intersect(substr(colnames(mut_count_filt2_df),1,12) , cli_surv$submitter_id)]
    
  } 
  
  if (all.equal(ncol(mut_count_filt2_df),  length(cli_surv_filt$submitter_id) )) {
    mut_count_filt2_df = as.data.frame(t(mut_count_filt2_df))
    mut_count_filt2_df$vital_status = cli_surv_filt$vital_status
      
      }
  
  mut_filtered@clinical.data$vital_status = NA
  mut_filtered@clinical.data$submitter_id = substr(mut_filtered@clinical.data$Tumor_Sample_Barcode,1,12)
  mut_filtered@clinical.data$vital_status <- as.character(mut_filtered@clinical.data$vital_status)
  
  for (maf_patients in substr(mut_filtered@clinical.data$Tumor_Sample_Barcode,1,12)) {
    if (maf_patients  %in% cli_surv$submitter_id) {
      mut_filtered@clinical.data[which(mut_filtered@clinical.data$submitter_id == maf_patients),]$vital_status = 
        cli_surv[which(cli_surv$submitter_id == maf_patients),]$vital_status
    } 
    
  }
  
  # # maf.data_2@clinical.data = na.omit(maf.data_2@clinical.data)
  # mut_filtered@clinical.data$vital_status = ifelse(mut_filtered@clinical.data$vital_status == FALSE, 0 , 1)
  # mut_filtered@clinical.data$vital_status = as.factor(mut_filtered@clinical.data$vital_status)
  genes = intersect(hvgs ,rownames(mut_count))[1:250]
  
  # for count above median
  setwd("~/nas/04.Results/compare_surv/")
  png(filename = paste0(CancerType,"_hvgs_top250_oncoplot.png"),
      width = 25, height = 25,  units = "cm" ,pointsize = 12,
      bg = "white", res = 1200, family = "")
  
  plot_count = oncoplot(maf = mut_filtered, genes = genes,
                        clinicalFeatures = "vital_status",
                        sortByAnnotation = TRUE)
  print(plot_count)
  dev.off()
  
  remove(mut_count_filt_df,mut_count_filt2_df,mut_filtered,cli_surv_filt,cli_surv, sce_for_hvgs, sce)
}

#   
# # for complex_heatmap divide by a specific column annotation 
# col_colors = list(CMS = c("black", "yellow", "green", "red"))
# col_colors = list(CMS = c("black", "yellow", "green"))
# col_ann = subset(sce_count_t[order(sce_count_t$CMS),],select = CMS)
# names(col_colors$CMS) = unique(col_ann$CMS)
# sce_count_t = sce_count_t[rownames(col_ann),]
# cheatmap_mat = as.matrix(t(sce_count_t[,-length(sce_count_t)]))
# 
# library(ComplexHeatmap)
# ComplexHeatmap::pheatmap(cheatmap_mat, 
#                          column_split = col_ann$CMS,
#                          labels_col = "",
#                          show_rownames = T, 
#                          show_colnames = F, 
#                          annotation_col = col_ann,
#                          annotation_colors = col_colors,
#                          # clustering_method = "average",
#                          cluster_cols = T,
#                          cluster_rows = T)
# 
# 


