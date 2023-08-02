library(stringr)
library(data.table)
library(readxl)
library(dplyr)
library(readxl)
library(tidyverse)
library(rrvgo)

filepath = "/home/seokwon/nas/"
ref_path = paste0(filepath, "99.reference/")

link_genes = readRDS(paste0(ref_path, "/KEGG_pathway_shared_genes.rds"))
single_genes = readRDS(paste0(ref_path, "/Kegg_pathway_genes.rds"))
link_genes$n_genes = NULL
link_genes_filtered = link_genes[which(link_genes$shared_genes != ""),]

link_genes_filtered <- link_genes_filtered %>%
  mutate(shared_genes = strsplit(shared_genes, ",")) %>%
  unnest(cols = shared_genes)
link_genes_filtered_df = as.data.frame(link_genes_filtered)
colnames(link_genes_filtered_df) = colnames(single_genes)

Cancerlist = dir(paste0(filepath, "/00.data/filtered_TCGA/"))

surv_total_results = read_xlsx("~/nas/04.Results/Total_results_survpval2.xlsx")

# for fic
fig_path = paste0(filepath,"04.Results/bestfeatures/Genes")
if(!dir.exists(fig_path)){
  dir.create(fig_path)
  print(paste0("Created folder: ", fig_path))
} else {
  print(paste0("Folder already exists: ", fig_path))
}
setwd(fig_path)

# for all
for (num_CancerType in Cancerlist) {
  
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  
  # call input
  
  cancer_bf = read.csv(paste0(filepath,"04.Results/bestfeatures/",CancerType, "_best_features.csv"))
  cancer_bf_cut = cancer_bf[1:surv_total_results[which(surv_total_results$CancerType == CancerType),]$num_of_features,]
  
  write.csv(cancer_bf_cut,paste0(filepath,"04.Results/bestfeatures/",CancerType, "_critical_features.csv"))

  # 일단 unique gene으로 해봄
  total_link_genes = link_genes_filtered_df[which(link_genes_filtered_df$Pathway %in% cancer_bf_cut$variable),]$Genes
  total_single_genes = single_genes[which(single_genes$Pathway %in% cancer_bf_cut$variable),]$Genes
  total_bf_genes = c(total_link_genes,total_single_genes)
  total_bf_genes = data.frame(Genes = total_bf_genes)
  
  write.csv(total_bf_genes,paste0(CancerType, "_genes_from_critical_features.csv"),row.names = F)
  
}  


