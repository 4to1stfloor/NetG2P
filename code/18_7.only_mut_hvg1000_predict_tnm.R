
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
  mut = readRDS(paste0(sce_path,num_CancerType,"/", CancerType,"_maf.rds"))
  mut_filtered <- subsetMaf(maf = mut, query = "Variant_Classification != 'Nonsense_Mutation'")

  sce_for_hvgs = sce 
  
  sce_for_hvgs = logNormCounts(sce_for_hvgs)
  # sce_for_hvgs = subset(sce_for_hvgs, rownames(sce_for_hvgs) %in% V(g.ppi.conn)$name,)
  
  assays_to_remove <- c("counts", "snv", "cnv")
  for (assay_name in assays_to_remove) {
    assay(sce_for_hvgs, assay_name) <- NULL
  }
  
  dec.sce <- modelGeneVar(sce_for_hvgs)
  hvgs = getTopHVGs(dec.sce, n=1000)

  mut_filtered@clinical.data$pat_t = NA
  mut_filtered@clinical.data$pat_n = NA
  mut_filtered@clinical.data$pat_m = NA

  mut_filtered@clinical.data$pat_t <- as.character(mut_filtered@clinical.data$pat_t)
  mut_filtered@clinical.data$pat_n <- as.character(mut_filtered@clinical.data$pat_n)
  mut_filtered@clinical.data$pat_m <- as.character(mut_filtered@clinical.data$pat_m)
  
  mut_filtered@clinical.data$submitter_id = substr(mut_filtered@clinical.data$Tumor_Sample_Barcode,1,12)
  
  for (maf_patients in substr(mut_filtered@clinical.data$Tumor_Sample_Barcode,1,12)) {
    if (maf_patients  %in% cli_surv_filtered$submitter_id) {
      mut_filtered@clinical.data[which(mut_filtered@clinical.data$submitter_id == maf_patients),]$pat_t = 
        cli_surv_filtered[which(cli_surv_filtered$submitter_id == maf_patients),]$ajcc_pathologic_t
      
      mut_filtered@clinical.data[which(mut_filtered@clinical.data$submitter_id == maf_patients),]$pat_n = 
        cli_surv_filtered[which(cli_surv_filtered$submitter_id == maf_patients),]$ajcc_pathologic_n 
      
      mut_filtered@clinical.data[which(mut_filtered@clinical.data$submitter_id == maf_patients),]$pat_m = 
        cli_surv_filtered[which(cli_surv_filtered$submitter_id == maf_patients),]$ajcc_pathologic_m 
    } 
    
  }
  
  mut_filtered@clinical.data = mut_filtered@clinical.data[which(!is.na(mut_filtered@clinical.data$pat_t)),]
  mut_filtered@clinical.data = mut_filtered@clinical.data[which(!is.na(mut_filtered@clinical.data$pat_n)),]
  mut_filtered@clinical.data = mut_filtered@clinical.data[which(!is.na(mut_filtered@clinical.data$pat_m)),]
  
  # # maf.data_2@clinical.data = na.omit(maf.data_2@clinical.data)
  # mut_filtered@clinical.data$vital_status = ifelse(mut_filtered@clinical.data$vital_status == FALSE, 0 , 1)
  # mut_filtered@clinical.data$vital_status = as.factor(mut_filtered@clinical.data$vital_status)
  genes = intersect(hvgs ,mut_filtered@gene.summary$Hugo_Symbol)[1:250]
  
  # for count above median

  png(filename = paste0(CancerType,"_mut_top250_t_oncoplot.png"),
      width = 25, height = 25,  units = "cm" ,pointsize = 12,
      bg = "white", res = 1200, family = "")
  
  plot_count1 = oncoplot(maf = mut_filtered, genes = genes,
                        clinicalFeatures = "pat_t",
                        sortByAnnotation = TRUE)
  print(plot_count1)
  dev.off()
  
  png(filename = paste0(CancerType,"_mut_top250_n_oncoplot.png"),
      width = 25, height = 25,  units = "cm" ,pointsize = 12,
      bg = "white", res = 1200, family = "")
  
  plot_count2 = oncoplot(maf = mut_filtered, genes = genes,
                        clinicalFeatures = "pat_n",
                        sortByAnnotation = TRUE)
  print(plot_count2)
  dev.off()
  
  png(filename = paste0(CancerType,"_mut_top250_m_oncoplot.png"),
      width = 25, height = 25,  units = "cm" ,pointsize = 12,
      bg = "white", res = 1200, family = "")
  
  plot_count3 = oncoplot(maf = mut_filtered, genes = genes,
                        clinicalFeatures = "pat_m",
                        sortByAnnotation = TRUE)
  print(plot_count3)
  dev.off()
  
  remove(mut_count_filt_df,mut_count_filt2_df,mut_filtered,cli_surv_filt,cli_surv, sce_for_hvgs, sce)
}
