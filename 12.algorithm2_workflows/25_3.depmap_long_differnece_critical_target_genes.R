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

meta_for_TCGA = read.csv(paste0(ref_path, "/Depmap_meta_filt_to_TCGA.csv"))

## Read Dataset 
meta.tbl <- read.csv(file = paste0(ref_path,'DepMap/Model.csv'), 
                     sep=',', header=T, fill=T)

ge.tbl <- read.csv(file = paste0(ref_path, 'DepMap/CRISPRGeneEffect.csv'), 
                   sep=',', header=T, fill=T, row.names=1)

# for fic
fig_path = "~/nas/04.Results/drug/depmap/"
if(!dir.exists(fig_path)){
  dir.create(fig_path)
  print(paste0("Created folder: ", fig_path))
} else {
  print(paste0("Folder already exists: ", fig_path))
}
setwd(fig_path)

# num_CancerType = "04.TCGA-CESC"
for (num_CancerType in Cancerlist) {
  
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  ModelType = gsub('TCGA-','', CancerType)
  cancer_bf = read.csv(paste0("~/nas/04.Results/short_long/ttest_common/",CancerType,"_critical_features_short_long_common.csv"))
  
  # call input
  
  common_features = cancer_bf %>%
    filter( classification == "common")
  
  short_features = cancer_bf %>%
    filter( classification == "short")
  
  long_features = cancer_bf %>%
    filter( classification == "long")
  
  common_link_genes = link_genes_filtered_df[which(link_genes_filtered_df$Pathway %in% common_features$variable),]$Genes
  common_single_genes = single_genes[which(single_genes$Pathway %in% common_features$variable),]$Genes
  total_common_genes = c(common_link_genes,common_single_genes)
  
  short_link_genes = link_genes_filtered_df[which(link_genes_filtered_df$Pathway %in% short_features$variable),]$Genes
  short_single_genes = single_genes[which(single_genes$Pathway %in% short_features$variable),]$Genes
  total_short_genes = c(short_link_genes,short_single_genes)
  
  long_link_genes = link_genes_filtered_df[which(link_genes_filtered_df$Pathway %in% long_features$variable),]$Genes
  long_single_genes = single_genes[which(single_genes$Pathway %in% long_features$variable),]$Genes
  total_long_genes = c(long_link_genes,long_single_genes)
  
  total_difference_long_genes = setdiff(total_long_genes,total_short_genes)
  total_difference_long_genes = setdiff(total_difference_long_genes , total_common_genes)
  
  if (length(total_difference_long_genes) == 0) {
    next
  }
  ## Interested gene set
  depmap_common_genes <- colnames(ge.tbl)[which(colnames(ge.tbl) %in% total_difference_long_genes)]
  
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
  ge.sub <- ge.tbl[rownames(ge.tbl) %in% cell_oi, c('cell_id', depmap_common_genes)] 
  ge.interest <- melt(ge.sub, id.vars = 'cell_id')
  
  gene.filter.plot <- as.data.frame(matrix(nrow = length(depmap_common_genes), ncol = 2))
  colnames(gene.filter.plot) <- c("Genes", "mean_depmap")
  
  gene.filter.plot$Genes <- depmap_common_genes
  # Mean(depmap score)
  gene.filter.plot$mean_depmap <- sapply(1:length(depmap_common_genes), 
                                         function(x) -mean(dplyr::filter(ge.interest, variable == depmap_common_genes[x])$value) ) 
  
  gene_filtere_edit = gene.filter.plot %>% filter(mean_depmap > 0.5)
  
  ge_edit_specific = ge.interest[which(ge.interest$variable %in% gene_filtere_edit$Genes),]
  
  png(filename = paste0(CancerType,"_depmap_specific_difference_long_critical_features.png"),
      width = 35, height = 35,  units = "cm" ,pointsize = 12,
      bg = "white", res = 1200, family = "")
  
  ## Kenel Density Curve
  total_out = ggplot(ge_edit_specific, aes(x = value, y = variable, fill = after_stat(x))) +
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
