library(ggbiplot)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggfortify)
library(pheatmap)
library(RColorBrewer)
library(Rtsne)
# library(tsne)
# library(umap)
library(tidyr)
library(dplyr)
library(reshape2)
library(readxl)
library(ggridges)
library(ggdist)
library(ggpubr)

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

## Read Dataset 
meta.tbl <- read.csv(file = paste0(ref_path,'DepMap/Model.csv'), 
                     sep=',', header=T, fill=T)

ge.tbl <- read.csv(file = paste0(ref_path, 'DepMap/CRISPRGeneEffect.csv'), 
                   sep=',', header=T, fill=T, row.names=1)

meta_for_TCGA = read.csv(paste0(ref_path, "/Depmap_meta_filt_to_TCGA.csv"))

Cancerlist = dir(paste0(filepath, "/00.data/filtered_TCGA/"))
setwd("~/nas/04.Results/drug/depmap/")
####
Cancerlist = Cancerlist[c(-11,-12)]
# num_CancerType =  "04.TCGA-CESC"

for (num_CancerType in Cancerlist) {
  
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  Cancername = gsub('TCGA-' , '', CancerType)
  
  # call input
  
  gc_TCGA = readRDS(paste0("~/nas/04.Results/short_long/", CancerType,"_critical_features_short_long.rds"))
  
  short_cluster = gc_TCGA %>% filter(cluster == "short") %>% dplyr::select(-vitalstatus, -duration,-status, -cluster)
  long_cluster = gc_TCGA %>% filter(cluster == "long") %>% dplyr::select(-vitalstatus, -duration,-status, -cluster)
  wo_info = gc_TCGA %>% dplyr::select(-vitalstatus, -duration, -status)
  gc_cellline_filt_df = readRDS(paste0("~/nas/04.Results/drug/depmap/gdsc/",CancerType, "_DM_sl_cluster.rds"))

  #### validation 
  cellline_long = gc_cellline_filt_df %>% filter(cluster == "long")
  cellline_short = gc_cellline_filt_df %>% filter(cluster == "short")
  
  depmap_common_link_genes = link_genes_filtered_df %>% 
    filter(Pathway %in% colnames(gc_cellline_filt_df %>% select(-cluster))) %>%
    pull(Genes) %>%
    unique()
  
  depmap_common_each_genes = single_genes %>% 
    filter(Pathway %in% colnames(gc_cellline_filt_df %>% select(-cluster))) %>%
    pull(Genes) %>%
    unique()
  
  depmap_common_genes = unique(c(depmap_common_link_genes,depmap_common_each_genes))
  
  filtered_df = data.frame()
  for (critical_genes in depmap_common_genes) {
    
    ge_sub_short = ge.tbl %>% 
      filter(rownames(ge.tbl) %in% rownames(cellline_short)) %>%
      dplyr::select(cell_id,any_of(critical_genes))
    
    ge_sub_long = ge.tbl %>% 
      filter(rownames(ge.tbl) %in% rownames(cellline_long)) %>%
      dplyr::select(cell_id,any_of(critical_genes))
    
    ge_short_interest <- melt(ge_sub_short, id.vars = 'cell_id')
    ge_short_interest = ge_short_interest[!is.na(ge_short_interest$value),]
    ge_long_interest <- melt(ge_sub_long, id.vars = 'cell_id')
    ge_long_interest = ge_long_interest[!is.na(ge_long_interest$value),]
    
    if (!is.character(ge_long_interest) & !is.logical(ge_short_interest)) {
      if (nrow(ge_long_interest) <=1 | length(ge_long_interest) <=1 ) {
        next
      } else if (nrow(ge_short_interest) <=1 | length(ge_short_interest) <=1) {
        next
      } else {
        t_test_res =  t.test(ge_short_interest$value, ge_long_interest$value)
        
        if (t_test_res$p.value <0.05) {
          
          tmp_test_df = data.frame(genes = critical_genes,
                                   delta_long_to_short =  mean(ge_short_interest$value) - mean(ge_long_interest$value),
                                   sig = "significant",
                                   specific = ifelse(mean(ge_short_interest$value) - mean(ge_long_interest$value) > 0 ,
                                                     "long", "short"))
          
          filtered_df = rbind(filtered_df, tmp_test_df)
          
        } 
      }
    } else {
      next
    }
    
  }
  
  if (nrow(filtered_df) == 0) {
    next
  } else {
    filtered_ordered_df = filtered_df %>% arrange(delta_long_to_short)
  }
  
  tmp_all = filtered_ordered_df %>% filter(sig == "significant")
  
  write.csv(tmp_all , paste0("~/nas/04.Results/drug/depmap/",Cancername, "_depmap_sig_total.csv"))
  
}
  