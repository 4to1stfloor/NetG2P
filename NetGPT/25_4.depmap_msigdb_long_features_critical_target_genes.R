library(tidyverse)
library(ggplot2)
library(ggridges)
library(reshape2)
library(readxl)
library(GSA)
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

gmt_file <-GSA.read.gmt(paste(ref_path,
                              "db/c5.go.bp.v2023.1.Hs.symbols.gmt",
                              sep = "/"))

# for fic
fig_path = "~/nas/04.Results/drug/depmap/"
if(!dir.exists(fig_path)){
  dir.create(fig_path)
  print(paste0("Created folder: ", fig_path))
  } else {
    print(paste0("Folder already exists: ", fig_path))
    }
setwd(fig_path)

transform_string <- function(s) {
  return(paste("GOBP_", gsub(" ", "_", toupper(s)), sep = ""))
}


# num_CancerType = "04.TCGA-CESC"
for (num_CancerType in Cancerlist) {
  
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  ModelType = gsub('TCGA-','', CancerType)
  cancer_bf = read.csv(paste0("~/nas/04.Results/short_long/ttest_common/",CancerType,"_critical_features_short_long_common.csv"))
  
  # total_link_genes = link_genes_filtered_df[which(link_genes_filtered_df$Pathway %in% cancer_bf$variable),]$Genes
  # total_single_genes = single_genes[which(single_genes$Pathway %in% cancer_bf$variable),]$Genes
  # total_bf_genes = c(total_link_genes,total_single_genes)
  # total_bf_genes
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

  
  depmap_pathwaylinks_df = link_genes_filtered_df[which(link_genes_filtered_df$Genes %in% unique(ge_edit_specific$variable) ),]
  depmap_pathway_df = single_genes[which(single_genes$Genes %in% unique(ge_edit_specific$variable) ),]
  
  # write_csv(depmap_pathwaylinks_df, paste0(CancerType,"_depmap_pathwaylinks.csv"))
  # write_csv(depmap_pathway_df, paste0(CancerType,"_depmap_pathway.csv"))
  # 
  # for fic
  fig_path = paste0(filepath,"04.Results/GOenrichment_test/treemap_edit/five_category/")
  setwd(fig_path)

  short_long_features = read_csv(paste0(filepath, "04.Results/short_long/ttest_common/",CancerType, "_critical_features_short_long_common.csv"))
  
  for (class_feature in unique(short_long_features$classification)) {
    assign(paste0(class_feature,"_reducedGO"),  readRDS(paste0(CancerType,"_","reducedTerms",class_feature,".rds")))

  }
  

  long_pt_tf <- as.vector(sapply(long_reducedGO$parentTerm, transform_string))
  
  exist_db_long = gmt_file$geneset.names[which(gmt_file$geneset.names %in% unique(long_pt_tf))]
  
  n = 0
  genelist_from_msig = list()
  for (num_bp in which(gmt_file$geneset.names %in% exist_db_long)) {
    n = n+1
    genelist_from_msig[[exist_db_long[n]]] = gmt_file$genesets[[num_bp]]
  }
  
  depmap_to_bpGO <- data.frame(genes = character(0), GO_name = character(0))
  
  for (key in names(genelist_from_msig)) {
    matching_genes <- unique(ge_edit_specific$variable)[unique(ge_edit_specific$variable) %in% genelist_from_msig[[key]]]
    
    if (length(matching_genes) > 0) {
      new_row <- data.frame(genes = paste(matching_genes, collapse = " "), GO_name = key)
      depmap_to_bpGO <- rbind(depmap_to_bpGO, new_row)
    }
  }
  
  fig_path = "~/nas/04.Results/drug/depmap/"
  setwd(fig_path)
  write_csv(depmap_to_bpGO, paste0(CancerType,"_depmap_sig_genes_to_long_msigdb.csv"))

}
