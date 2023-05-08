library(TCGAbiolinks)
library(SummarizedExperiment)
library(SingleCellExperiment)
library(maftools)
library(dplyr)
library(DT)
library(biomaRt)
library(data.table)
# maf file load

filepath = "/home/seokwon/nas/"
ref_path = paste0(filepath, "99.reference/")
Cancerlist = dir(paste0(filepath, "/00.data/filtered_TCGA/"))
Cancerlist = Cancerlist[-3]
cli = fread(paste0(ref_path , "all_clin_indexed.csv"))

total_spe_pval = read.csv(paste0(filepath,"/04.Results/Total_results_survpval.csv"))
dual_total_genes = readRDS(paste0(ref_path, "/KEGG_dual_total_genes.rds"))

# for all
for (num_CancerType in Cancerlist) {
  
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  
  best_features = readRDS(paste0(main.path_tc,"/h2o_bias_pval_dual_cut_50/network/",CancerType,"_best_features_links.rds"))
  best_features$pathwaylinks = ifelse(best_features$from == best_features$to,best_features$from,  paste0(best_features$from, best_features$to))
  best_features = best_features[1:total_spe_pval[which(total_spe_pval$CancerType == CancerType),]$num_of_features,]

  # only has links
  
  total_dual_genes <- dual_total_genes[which(dual_total_genes$Pathway %in% best_features$pathwaylinks),]$Genes

  # load maf file
  mafdata = readRDS(paste0("/mnt/gluster_server/data/raw/TCGA_data/00.data/",num_CancerType,"/",CancerType,"_maf.rds"))
  mafdata <- subsetMaf(maf = mafdata, query = "Variant_Classification != 'Nonsense_Mutation'")
  
  # add vital_status to maffile 
  
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
  
  cli_surv = cli_surv[!is.na(cli_surv$overall_survival),]
  cli_surv = cli_surv[!is.na(cli_surv$vital_status),]
  cli_surv = cli_surv[!which(cli_surv$vital_status == "Not Reported"),]
  cli_surv = cli_surv[!which(cli_surv$vital_status == "NA"),]
  
  cli_surv$status = NA
  cli_surv$status[which(cli_surv$vital_status == "Dead")] = 1
  cli_surv$status[which(cli_surv$vital_status == "Alive")] = 0
  
  cli_surv$status = as.numeric(cli_surv$status)
  cli_surv$submitter_id = as.character(cli_surv$submitter_id)
  cli_surv$vital_status = as.character(cli_surv$vital_status)
  
  mafdata@clinical.data$vital_status = NA
  mafdata@clinical.data$vital_status = as.character(mafdata@clinical.data$vital_status)
  
  mafdata@clinical.data$submitter_id = substr(mafdata@clinical.data$Tumor_Sample_Barcode,1,12)
  
  for (maf_patients in mafdata@clinical.data$submitter_id ) {
    if (maf_patients  %in% cli_surv$submitter_id) {
      mafdata@clinical.data[which(mafdata@clinical.data$submitter_id == maf_patients),]$vital_status = cli_surv[which(cli_surv$submitter_id == maf_patients),]$vital_status
    } 
    
  } 

  mafdata@clinical.data = mafdata@clinical.data[!is.na(mafdata@clinical.data$vital_status),]
  # maf.data_2@clinical.data = na.omit(maf.data_2@clinical.data)
  # mafdata@clinical.data$vital_status = ifelse(mafdata@clinical.data$vital_status == FALSE, 0 , 1)
  mafdata@clinical.data$vital_status = as.factor(mafdata@clinical.data$vital_status)
  
  # formation of count genes
  uni_genes = unique(total_dual_genes)
  
  bf_genes_count = data.frame(name_genes = uni_genes)
  bf_genes_count$count = NA
  
  for (genes in uni_genes) {
    bf_genes_count[which(bf_genes_count$name_genes == genes),]$count = sum(total_dual_genes %in% genes)
  }
  
  bf_genes_count = bf_genes_count[order(-bf_genes_count$count),]
  row.names(bf_genes_count) <- NULL

  dirpath = paste0(main.path_tc,"/",CancerType,"_analysis_mutation_bf/",CancerType,"_analysis_mutation_from_surv/")
  
  if(!dir.exists(dirpath)){
    dir.create(dirpath)
    print(paste0("Created folder: ", dirpath))
  } else {
    print(paste0("Folder already exists: ", dirpath))
  }

  setwd(dirpath)
  
  write.csv(bf_genes_count, paste0(dirpath,CancerType,"_counts_gene_for_bf.csv"))
  
  # for count above median
  png(filename = paste0(CancerType,"_counts_oncoplot.png"),
      width = 25, height = 25,  units = "cm", pointsize = 12,
      bg = "white", res = 1200, family = "")
  
  plot_count = oncoplot(maf = mafdata, 
                        genes = rownames(mutCountMatrix(mafdata))[rownames(mutCountMatrix(mafdata)) %in% bf_genes_count[which(bf_genes_count$count > median(unique(bf_genes_count$count))),]$name_genes])
  
  print(plot_count)
  dev.off()
  
  # for pathway
  # subset only bestfeatures(pathway + pathwaylink) for each cancer
  
  bf_dual_genes = subset(dual_total_genes, dual_total_genes$Pathway %in% best_features$pathwaylinks)
  
  bf_dual_genes_ordered = data.frame()

  # re-order by importance of features
  for (bf_pathway in best_features$pathwaylinks) {
    bf_dual_genes_ordered = rbind(bf_dual_genes_ordered , bf_dual_genes[which(bf_dual_genes$Pathway == bf_pathway),])
  }
  
  bf_dual_genes_ordered = data.frame(bf_dual_genes_ordered, row.names = NULL)
  bf_dual_genes_ordered_filtered = bf_dual_genes_ordered[which(bf_dual_genes_ordered$Genes %in% mafdata@data$Hugo_Symbol),]
  row.names(bf_dual_genes_ordered_filtered) <- NULL
  uni_dualpathway = unique(bf_dual_genes_ordered_filtered$Pathway)[1:2]
  
  bf_dual_genes_ordered_filtered$Genes %in% mafdata@data
  bf_dual_genes_ordered_filtered$Genes %in% mafdata@data$Hugo_Symbol 
  bf_dual_genes_ordered_filtered$Genes %in% rownames(mutCountMatrix(mafdata))
  
  maf_subset <- subset(mafdata@data, Hugo_Symbol %in% bf_dual_genes_ordered_filtered$Genes)
  maf_subset <- MAF(maf_subset)
  
  # top 2 pathway
  png(filename = paste0(CancerType,"_dual_oncoplot.png"),
      width = 35, height = 50,  units = "cm", pointsize = 12,
      bg = "white", res = 1200, family = "")

  plot_pathwaylink = oncoplot(maf = mafdata,
                              pathways = bf_dual_genes_ordered_filtered[which(bf_dual_genes_ordered_filtered$Pathway %in% uni_dualpathway),],
                              clinicalFeatures = "vital_status",
                              sortByAnnotation = TRUE)
  colnames(maf_subset@data)[2] <- "Hugo_Symbol"
  maf_subset@data$Hugo_Symbol = as.factor(maf_subset@data$Hugo_Symbol)
  oncoplot(maf = maf_subset,
           pathways = bf_dual_genes_ordered_filtered[which(bf_dual_genes_ordered_filtered$Pathway %in% uni_dualpathway),],
           clinicalFeatures = "vital_status",
           sortByAnnotation = TRUE)
  
  print(plot_pathwaylink)
  dev.off()

  remove(bf_dual_genes_ordered,best_features,mafdata,plot_pathwaylink,cli_surv,bf_genes_count)
}
