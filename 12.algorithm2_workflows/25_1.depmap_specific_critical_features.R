## DepMap case1: density ridgeline plots
## Date: 2023-07-21
## Juyeon Cho

library(tidyverse)
library(ggplot2)
library(ggridges)
library(reshape2)
library(readxl)
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
cancer_specific_critical_features = read_csv("~/nas/04.Results/critical_features/unique_short_long_critical_features_for_SW.csv")
meta_for_TCGA = read.csv(paste0(ref_path, "/Depmap_meta_filt_to_TCGA.csv"))

# for fic
fig_path = "~/nas/04.Results/drug/depmap/"
if(!dir.exists(fig_path)){
  dir.create(fig_path)
  print(paste0("Created folder: ", fig_path))
} else {
  print(paste0("Folder already exists: ", fig_path))
}
setwd(fig_path)

## Read Dataset 
meta.tbl <- read.csv(file = '/mnt/gluster_server/data/raw/DepMap/23Q2/Model.csv', 
                     sep=',', header=T, fill=T)

ge.tbl <- read.csv(file = '/mnt/gluster_server/data/raw/DepMap/23Q2/CRISPRGeneEffect.csv', 
                   sep=',', header=T, fill=T, row.names=1)
# rnai
# ge.tbl <- read.csv(file = '/mnt/gluster_server/data/raw/DepMap/23Q2/', 
#                    sep=',', header=T, fill=T, row.names=1)

## Organize gene name for Gene Effect table
ge.tbl$cell_id <- rownames(ge.tbl)
all_gene <- colnames(ge.tbl)
match_gene <- sub("\\..*","",all_gene)
colnames(ge.tbl) <- match_gene
# num_CancerType = "04.TCGA-CESC"
for (num_CancerType in Cancerlist) {
  
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  ModelType = gsub('TCGA-','', CancerType)
  
  # call input
  
  if (sum(cancer_specific_critical_features$which_cancer == CancerType) == 0) {
    next
  } else {
    
    cancer_critical_features = cancer_specific_critical_features[which(cancer_specific_critical_features$which_cancer == CancerType),]
    cancer_critical_features = cancer_critical_features[which(cancer_critical_features[,CancerType] != "common"),]
    total_link_genes = link_genes_filtered_df[which(link_genes_filtered_df$Pathway %in% cancer_critical_features$features),]$Genes
    total_single_genes = single_genes[which(single_genes$Pathway %in% cancer_critical_features$features),]$Genes
    interested_gene_set = c(total_link_genes,total_single_genes)
  }

  ## Interested gene set
  geneset_oi <- colnames(ge.tbl)[which(colnames(ge.tbl) %in% interested_gene_set)]
  
  #subset(GrowthPattern != 'Organoid' & DepmapModelType == 'COAD') %>%
  #subset(GrowthPattern != 'Organoid' & DepmapModelType == 'COAD' & PrimaryOrMetastasis == 'Primary') %>%
  # subset(GrowthPattern != 'Organoid' & DepmapModelType == 'BLCA' & PrimaryOrMetastasis == 'Primary' & SampleCollectionSite == 'Colon') %>%
  
  ## Select cells you want
  cell_oi <- meta.tbl %>% 
    #subset(GrowthPattern != 'Organoid' & DepmapModelType == 'COAD') %>%
    #subset(GrowthPattern != 'Organoid' & DepmapModelType == 'COAD' & PrimaryOrMetastasis == 'Primary') %>%
    # subset(GrowthPattern != 'Organoid' & DepmapModelType == 'BLCA' & PrimaryOrMetastasis == 'Primary' & SampleCollectionSite == 'Colon') %>%
    subset(GrowthPattern != 'Organoid' & 
             PrimaryOrMetastasis == 'Primary' &
             OncotreeLineage == unique(meta_for_TCGA[which(meta_for_TCGA$TCGA == CancerType),]$oncotree) &
             OncotreeCode %in% unique(meta_for_TCGA[which(meta_for_TCGA$TCGA == CancerType),]$oncotreecode)) %>%
    dplyr::select(ModelID) %>% unlist(use.names=F)
  
  ## Subset dataset
  ge.sub <- ge.tbl[rownames(ge.tbl) %in% cell_oi, c('cell_id', geneset_oi)] 
  ge.interest <- melt(ge.sub, id.vars = 'cell_id')
  ge.interest[which(ge.interest$value == min(ge.interest$value)),]
  
  png(filename = paste0(CancerType,"_depmap_specific_critical_features.png"),
      width = 35, height = 35,  units = "cm" ,pointsize = 12,
      bg = "white", res = 1200, family = "")
  
  ## Kenel Density Curve
  total_out = ggplot(ge.interest, aes(x = value, y = variable, fill = after_stat(x))) +
    geom_density_ridges_gradient(scale = 2, gradient_lwd = .5, color = "black",
                                 jittered_points = TRUE,
                                 position = position_points_jitter(width = 0.05, height = 0),
                                 point_shape = '|', point_size = 3, point_alpha = 1, alpha = 0.7) +
    scale_fill_viridis_c(option = "plasma", name = "Gene Effect") +
    geom_vline(xintercept = -1, color='darkred')+
    geom_vline(xintercept = c(1,0,-2), color = 'darkgrey') +
    geom_vline(xintercept = c(-2.5, -1.5, -0.5, 0.5), color = 'darkgrey', linetype = 'dashed') +
    
    ylab('Gene') +
    xlab('Gene Effect \n(Chronos score, 23Q2)') +
    theme_bw()
  
  print(total_out)
  
  dev.off()
  

  
}
