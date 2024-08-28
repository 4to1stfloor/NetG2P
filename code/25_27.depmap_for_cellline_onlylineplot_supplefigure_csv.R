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
# num_CancerType =  "26.TCGA-LUSC"
# Cancerlist = c("10.TCGA-BLCA","26.TCGA-LUSC", "32.TCGA-UCEC")
for (num_CancerType in Cancerlist) {
  
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  Cancername = gsub('TCGA-' , '', CancerType)
  
  # call input
  
  gc_TCGA = readRDS(paste0("~/nas/04.Results/short_long/", CancerType,"_critical_features_short_long.rds"))
  gc_cellline_filt_df = readRDS(paste0("~/nas/04.Results/drug/depmap/gdsc/",CancerType, "_DM_sl_cluster.rds"))
  
  #### validation 
  cellline_long = gc_cellline_filt_df %>% filter(cluster == "long")
  cellline_short = gc_cellline_filt_df %>% filter(cluster == "short")
  
  depmap_common_link_genes <- unique(link_genes_filtered_df[which(link_genes_filtered_df$Pathway %in% colnames(gc_cellline_filt_df %>% dplyr::select(-cluster))),]$Genes)
  depmap_common_each_genes <- unique(single_genes[which(single_genes$Pathway %in% colnames(gc_cellline_filt_df %>% dplyr::select(-cluster))),]$Genes)
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
        
        if (Cancername %in% c("BLCA","LUSC", "UCEC")) {
          pvalue_cutoff = 0.01
        } else {
          pvalue_cutoff = 0.05
        }
        
        if (t_test_res$p.value < pvalue_cutoff) {
          tmp_test_df = data.frame(genes = critical_genes,
                                   delta_long_to_short =  mean(ge_long_interest$value) - mean(ge_short_interest$value),
                                   mean_of_long = mean(ge_long_interest$value),
                                   mean_of_short = mean(ge_short_interest$value),
                                   pvalue = t_test_res$p.value,
                                   sig = "significant")
          
          filtered_df = rbind(filtered_df, tmp_test_df)
          
        } else {
          tmp_test_df = data.frame(genes = critical_genes,
                                   delta_long_to_short = mean(ge_long_interest$value) - mean(ge_short_interest$value),
                                   mean_of_long = mean(ge_long_interest$value),
                                   mean_of_short = mean(ge_short_interest$value),
                                   pvalue = t_test_res$p.value,
                                   sig = "not")
          
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
  
  library(ggrepel)
  
  filtered_ordered_df = filtered_ordered_df %>%
    mutate( label = case_when(sig == "significant" ~ genes,
                              .default = NA),
            alpha = case_when(sig == "significant" ~ T,
                              .default = F),
            logp = case_when(sig == "significant" ~ -log2(pvalue),
                             .default = 0))
  
  manual_size = ifelse(filtered_ordered_df$logp != 0 , filtered_ordered_df$logp -2, 0.5)
  filtered_ordered_df = filtered_ordered_df %>% 
    mutate( direction = case_when(delta_long_to_short > 0 &
                                    mean_of_long > mean_of_short & 
                                    sig == "significant" ~ "short_sig",
                                  delta_long_to_short > 0 &
                                    mean_of_long < mean_of_short & 
                                    sig == "significant" ~ "long_sig",
                                  
                                  delta_long_to_short < 0 &
                                    mean_of_long > mean_of_short &
                                    sig == "significant" ~ "short_sig",
                                  delta_long_to_short < 0 &
                                    mean_of_long < mean_of_short &
                                    sig == "significant" ~ "long_sig",
                                  .default = "nonsig"))

  tmp_demap_df = filtered_ordered_df %>% 
    select(genes, delta_long_to_short, mean_of_long, mean_of_short, pvalue)
  
  write_csv(tmp_demap_df, paste0(Cancername, "_depmap_total_screen.csv"))

  
}


