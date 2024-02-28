TLDR.sample.finder = function(tissue = NULL, TCGACode = NULL, OncoTreeCode = NULL,
                              snv.gene = NULL) {
  #define database
  meta.cell = readRDS("/mnt/gluster_server/data/reference/TLDR/meta_cells_primary.rds")
  meta.gene = readRDS("/mnt/gluster_server/data/reference/TLDR/meta_genes.rds")
  meta.map = readRDS("/mnt/gluster_server/data/reference/TLDR/meta_TCGA_Oncotree.rds")
  #creating arrays for term validation
  tcga.codes = unique(meta.map$TCGACode)
  oncotree.codes = unique(unlist(strsplit(meta.map$OncotreeCode, ","))) #might have to go deeper and read from meta.cell
  tissue.type = unique(meta.map$tissue)
  #entry check
  tissue.arr = c(is.null(tissue), is.null(TCGACode), is.null(OncoTreeCode))
  if (sum(tissue.arr) == 3) {
    stop("You must provide a category between tissue, TCGA and OncoTreeCode")
  } else if (sum(tissue.arr) < 2) {
    stop("TLDR currently only supports one sample type query.")
  }
  #i can probably compress this, but later
  #validity check
  if (!is.null(tissue)) {
    valid.tissue = intersect(tissue.type, tissue)
    if (length(valid.tissue) == 0) {
      stop("There are no valid tissue types in your query. Valid tissue types are: ", paste(unique(tissue.type), collapse = ", "))
    }
    meta.map.f = meta.map[meta.map$tissue %in% valid.tissue,]
  }
  if (!is.null(TCGACode)) {
    valid.TCGA = intersect(tcga.codes, TCGACode)
    if (length(valid.TCGA) == 0) {
      stop("There are no valid TCGA codes in your query. Valid TCGA codes are: ", paste(tcga.codes, collapse = ", "))
    }
    meta.map.f = meta.map[meta.map$TCGACode %in% valid.TCGA,]
  }
  if (!is.null(OncoTreeCode)) {
    valid.onco = intersect(oncotree.codes %in% OncoTreeCode)
    if (length(valid.onco) == 0) {
      stop("There are no validOncoTreeCodes in your query. If you are unsure, try looking with tissue types or TCGA codes")
    }
    meta.map.f = meta.map[meta.map$OncotreeCode %in% valid.onco,]
  }
  message("Your query resulted in TCGA: ", paste0(meta.map.f$TCGACode, collapse = ", "))
  message("Your query resulted in OncoTreeCode: ", paste0(meta.map.f$OncotreeCode, collapse = ", "))
  message("Your query resulted in lineage: ", paste0(meta.map.f$OncotreeLineage, collapse = ", "))
  message("Searching for cells!")
  meta.cell.t = meta.cell[meta.cell$OncotreeLineage %in% meta.map.f$tissue,]
  if (is.null(tissue)) {
    message("Further filtering by OncoTreeCode...")
    split.onco = unlist(strsplit(meta.map.f$OncotreeCode, ","))
    meta.cell.t$target.cell = meta.cell.t$OncotreeCode %in% split.onco
    meta.cell.t$target.cell[meta.cell.t$OncotreePrimaryDisease == "Non-Cancerous"] = T
    meta.cell.f = meta.cell.t[meta.cell.t$target.cell,]
  } else {
    meta.cell.f = meta.cell.t
  }
  message(sum(meta.cell.f$DepMap), " cell lines in DepMap, ", sum(meta.cell.f$COSMIC), " cell lines in COSMIC, ",sum(meta.cell.f$LINCS), " cell lines in LINCS.")
  message("Returning the queried information as data.frame. Good luck!")
  return(meta.cell.f)
}

###

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
library(reshape2)
library(readxl)
# library(ggridges)
library(ggdist)
library(ggpubr)
library(tidyverse)
library(tidygraph)
library(ggraph)
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
gdsc = readRDS("/mnt/gluster_server/data/reference/GDSC/2024_01_11/GDSC_data_combined.rds")
meta_cell = readRDS("/mnt/gluster_server/data/reference/TLDR/meta_cells_primary.rds")
nih_drug = read_xlsx("~/nas/99.reference/Nih_drug_data.xlsx")

nih_drug_filted = nih_drug %>%
  separate_rows(Drug_Name, sep = "\\(") %>%
  mutate(Drug_Name = str_trim(gsub("\\)", "", Drug_Name))) %>%
  group_by(Tissue_name) %>%
  distinct(Drug_Name, .keep_all = TRUE)

anticancer_drug = read_xlsx("~/nas/99.reference/anticancer_fund_cancerdrugsdb.xlsx")

anticancer_drug_filted = anticancer_drug %>%
  separate_rows(Indications, sep = "\\;") %>%
  mutate(Indications = str_trim(Indications))

anticancer_drug_filted = anticancer_drug_filted[!is.na(anticancer_drug_filted$Indications),]
anticancer_drug_filted = anticancer_drug_filted %>% select(Product,Indications,Targets,CESC,BLCA,STAD,LUAD,LIHC,OV,LUSC,LGG,BRCA,UCEC,COADREAD,KIDNEY, everything())

Cancerlist = dir(paste0(filepath, "/00.data/filtered_TCGA/"))
setwd("~/nas/04.Results/drug/depmap/gdsc/")
####
Cancerlist = Cancerlist[c(-11,-12)]

# TLDR.sample.finder(TCGACode = "TCGA-CESC") %>% filter(DepMap == T & COSMIC == T & GDSC == T)

num_CancerType = "29.TCGA-LGG"

for (num_CancerType in Cancerlist) {
  
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  Cancername = gsub('TCGA-' , '', CancerType)
  
  # call input
  gc_cellline = readRDS(paste0("~/nas/00.data/filtered_TCGA/", num_CancerType, "/",Cancername,"_cellline_dual_all_log.rds"))
  gc_TCGA = readRDS(paste0("~/nas/04.Results/short_long/", CancerType,"_critical_features_short_long.rds"))
  
  filt_cancer_cell = meta_cell %>% 
    filter(DepMap.ID %in% rownames(gc_cellline)) %>%
    filter(DepMap == T & GDSC == T)
  
  setwd("~/nas/04.Results/drug/depmap/gdsc/")
  gdsc_each = gdsc %>% filter(COSMIC.ID %in% filt_cancer_cell$COSMIC.ID)
  
  gc_cellline_filt_df = readRDS(paste0(CancerType, "_DM_sl_cluster.rds"))
  
  # gdsc_each
  # filt_cancer_cell
  
  gc_filtered_cellline_df = gc_cellline_filt_df %>% filter(rownames(.) %in% filt_cancer_cell$DepMap.ID)
  
  gdsc_each_filt = left_join(gdsc_each, filt_cancer_cell %>% select(COSMIC.ID, DepMap.ID), by = "COSMIC.ID")
  
  tmp_df = data.frame(DepMap.ID = rownames(gc_filtered_cellline_df), cluster = gc_filtered_cellline_df$cluster)
  gdsc_w_cluster = left_join(gdsc_each_filt, tmp_df, by = "DepMap.ID")
  
  gdsc_w_cluster = gdsc_w_cluster %>%
    select(DepMap.ID, cluster, DRUG_ID,DRUG_NAME, PUTATIVE_TARGET, PATHWAY_NAME, AUC, RMSE, LN_IC50,Z_SCORE) %>% 
    arrange(DepMap.ID)
  
  tmp_anti = anticancer_drug_filted %>% 
    select(Product, Indications, Targets, 
                                    colnames(anticancer_drug_filted)[colnames(anticancer_drug_filted) == Cancername])
  
  tmp_anti_filt = tmp_anti[which(tmp_anti[,Cancername] == "Y"),]
                                    
  tmp_drug = gdsc_w_cluster %>% filter(DRUG_NAME %in% tmp_anti_filt$Product)

  library(ggpubr)
  library(ggsignif)
  
  # for fic
  fig_path = paste0(filepath,"/04.Results/drug/depmap/gdsc/", Cancername)
  if(!dir.exists(fig_path)){
    dir.create(fig_path)
    print(paste0("Created folder: ", fig_path))
  } else {
    print(paste0("Folder already exists: ", fig_path))
  }
  setwd(fig_path)
  
  total_long_gd = tmp_drug %>% filter(cluster == "long")
  total_short_gd = tmp_drug %>% filter(cluster == "short")
  
  if (nrow(total_long_gd) >1 & nrow(total_short_gd) >1) {
    anno = t.test(total_long_gd$Z_SCORE, total_short_gd$Z_SCORE)$p.value
    
    total_gd = ggplot(tmp_drug , aes( x = cluster , y = Z_SCORE, fill=cluster)) + 
      geom_violin(color ="black") +
      geom_boxplot(width=0.1, color = "black" , fill="white")+
      # scale_color_manual(values="black","black") + 
      scale_fill_manual(values=c("#4DAF4A", "#E41A1C")) +
      
      geom_signif(
        annotation = paste0(formatC(anno, digits = 3)),
        map_signif_level = TRUE, 
        comparisons = list(c("long", "short")),
        # y_position = 4.05, xmin = 1, xmax = 3,
        tip_length = c(0.22, 0.02),
      ) +
      # stat_compare_means(label.y = 10) +
      theme_minimal()
    
    ggsave(filename = paste0(CancerType,"_total_anti_cellline_IC50_zscore.svg"), total_gd)
  } else {
    
    next
  }
}
