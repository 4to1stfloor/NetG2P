library(TCGAbiolinks)
library(SummarizedExperiment)
library(SingleCellExperiment)
library(maftools)
library(dplyr)
library(DT)
library(biomaRt)

# maf file load

filepath = "/home/seokwon/nas/"
ref_path = paste0(filepath, "99.reference/")
Cancerlist = dir(paste0(filepath, "/00.data/filtered_TCGA/"))

cli = fread(paste0(ref_path , "all_clin_indexed.csv"))

total_spe_pval = read.csv(paste0(filepath,"/04.Results/Total_results_survpval.csv"))
onco_for_pathway_gene = readRDS(paste0(ref_path, "/KEGG_pathway_shared_each_gene.rds"))
pathwayeach_gene = readRDS(paste0(ref_path, "/Kegg_pathway_genes.rds"))
link_genes = readRDS(paste0(ref_path, "/KEGG_pathway_shared_genes.rds"))
link_genes_mode2 = readRDS("~/nas/99.reference/KEGG_pathway_shared_each_gene.rds")
dual_total_genes = readRDS(paste0(ref_path, "/KEGG_dual_total_genes.rds"))

Cancerlist = Cancerlist[7:12]
# for all
for (num_CancerType in Cancerlist) {
  
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  
  best_features = readRDS(paste0(main.path_tc,"/h2o_bias_pval_dual_cut_50/network/",CancerType,"_best_features_links.rds"))
  best_features$pathwaylinks = ifelse(best_features$from == best_features$to,best_features$from,  paste0(best_features$from, best_features$to) )
  best_spe_features_surv = total_spe_pval[which(total_spe_pval$CancerType == CancerType),]$num_of_features

  # only has links
  
  total_link_genes <- dual_total_genes[which(dual_total_genes$Pathway %in% best_features$pathwaylinks[1:best_spe_features_surv]),]$Genes
  
  # load maf file
  maf.data_2 = readRDS(paste0("/mnt/gluster_server/data/raw/TCGA_data/00.data/",num_CancerType,"/",CancerType,"_maf.rds"))
  
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
  
  cli_surv$status = NA
  cli_surv$status[which(cli_surv$vital_status == "Dead")] = 1
  cli_surv$status[which(cli_surv$vital_status == "Alive")] = 0
  
  cli_surv$status = as.numeric(cli_surv$status)
  
  maf.data_2@clinical.data$vital_status = NA
  maf.data_2@clinical.data$submitter_id = substr(maf.data_2@clinical.data$Tumor_Sample_Barcode,1,12)
  
  for (maf_patients in substr(maf.data_2@clinical.data$Tumor_Sample_Barcode,1,12)) {
    if (maf_patients  %in% cli_surv$submitter_id) {
      maf.data_2@clinical.data[which(maf.data_2@clinical.data$submitter_id == maf_patients),]$vital_status = cli_surv[which(cli_surv$submitter_id == maf_patients),]$status
    } 
    
  }
  # maf.data_2@clinical.data = na.omit(maf.data_2@clinical.data)
  maf.data_2@clinical.data$vital_status = ifelse(maf.data_2@clinical.data$vital_status == FALSE, 0 , 1)
  maf.data_2@clinical.data$vital_status = as.factor(maf.data_2@clinical.data$vital_status)
  
  # formation of count genes
  uni_genes = unique(total_link_genes)
  
  bf_genes_count = data.frame(name_genes = uni_genes)
  bf_genes_count$count = NA
  
  for (genes in uni_genes) {
    bf_genes_count[which(bf_genes_count$name_genes == genes),]$count = sum(total_link_genes %in% genes)
  }
  
  bf_genes_count = bf_genes_count[order(-bf_genes_count$count),]
  row.names(bf_genes_count) <- NULL

  dirpath = paste0(main.path_tc,"/",CancerType,"_analysis_mutation_bf/",CancerType,"_analysis_mutation_from_surv/")
  dir.create(dirpath)
  setwd(dirpath)
  
  write.csv(bf_genes_count, paste0(dirpath,CancerType,"_counts_gene_for_bf.csv"))
  
  # for count above median
  png(filename = paste0(CancerType,"_counts_oncoplot.png"),
      width = 25, height = 25,  units = "cm", pointsize = 12,
      bg = "white", res = 1200, family = "")
  
  plot_count = oncoplot(maf = maf.data_2, genes = bf_genes_count[which(bf_genes_count$count > median(unique(bf_genes_count$count))),]$name_genes)
  print(plot_count)
  dev.off()
  
  # it will be changed -> do not respect surv features 
  # w/o only pathway
  best_features_onlylinks <- best_features[nchar(best_features$pathwaylinks) >= 4, ]
  
  # subset only bestfeatures for each cancer
  bf_links_genes = subset(onco_for_pathway_gene, onco_for_pathway_gene$Pathway %in% best_features_onlylinks$pathwaylinks)
  
  bf_links_genes = subset(dual_total_genes, dual_total_genes$Pathway %in% best_features$pathwaylinks[1:best_spe_features_surv])
  
  bf_links_genes_ordered = data.frame()
  
  # re-order by importance of features
  for (bf_pathway in best_features$pathwaylinks[1:best_spe_features_surv]) {
    bf_links_genes_ordered = rbind(bf_links_genes_ordered , bf_links_genes[which(bf_links_genes$Pathway == bf_pathway),])
  }
  
  bf_links_genes_ordered = data.frame(bf_links_genes_ordered, row.names = NULL)
  
  
  # top 10 pathway
  png(filename = paste0(CancerType,"_pathwaylink_oncoplot.png"),
      width = 35, height = 50,  units = "cm", pointsize = 12,
      bg = "white", res = 1200, family = "")
  
  plot_pathwaylink = oncoplot(maf = maf.data_2, 
                              pathways = bf_links_genes_ordered[which(bf_links_genes_ordered$Pathway %in% best_features$pathwaylinks[1:best_spe_features_surv]),], 
                              clinicalFeatures = "vital_status",
                              sortByAnnotation = TRUE)
  print(plot_pathwaylink)
  dev.off()

  remove(bf_links_genes_ordered,best_features_onlylinks,best_features,maf.data_2,plot_pathwaylink,cli_surv,bf_genes_count)
}

