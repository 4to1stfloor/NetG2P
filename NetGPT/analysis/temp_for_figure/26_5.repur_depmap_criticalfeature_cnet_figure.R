
library(readxl)
library(tidyverse)
library(svglite)
library(ggplot2)
library(ggsignif)
library(openxlsx)
library(ggpubr)
library(ggsignif)
library(ggraph)
library(tidygraph)

filepath = "/home/seokwon/nas/"
ref_path = paste0(filepath, "99.reference/")

Cancerlist = dir(paste0(filepath, "/00.data/filtered_TCGA/"))
Cancerlist = Cancerlist[c(-11,-12)]

link_genes = readRDS(paste0(ref_path, "/KEGG_pathway_shared_genes.rds"))
single_genes = readRDS(paste0(ref_path, "/Kegg_pathway_genes.rds"))
link_genes$n_genes = NULL
link_genes_filtered = link_genes[which(link_genes$shared_genes != ""),]

link_genes_filtered <- link_genes_filtered %>%
  mutate(shared_genes = strsplit(shared_genes, ",")) %>%
  unnest(cols = shared_genes)
link_genes_filtered_df = as.data.frame(link_genes_filtered)
colnames(link_genes_filtered_df) = colnames(single_genes)

#drug 
anticancer_drug = read_xlsx("~/nas/99.reference/anticancer_fund_cancerdrugsdb.xlsx")

anticancer_drug_filted = anticancer_drug %>%
  separate_rows(Indications, sep = "\\;") %>%
  mutate(Indications = str_trim(Indications))
nih_drug = read_xlsx("~/nas/99.reference/Nih_drug_data.xlsx")

nih_drug_filted = nih_drug %>%
  separate_rows(Drug_Name, sep = "\\(") %>%
  mutate(Drug_Name = str_trim(gsub("\\)", "", Drug_Name))) %>%
  group_by(Tissue_name) %>%
  distinct(Drug_Name, .keep_all = TRUE)

meta_for_TCGA = read.csv(paste0(ref_path, "/Depmap_meta_filt_to_TCGA.csv"))
gdsc = readRDS("/mnt/gluster_server/data/reference/GDSC/2024_01_11/GDSC_data_combined.rds")
meta_cell = readRDS("/mnt/gluster_server/data/reference/TLDR/meta_cells_primary.rds")

criteria_filt = read_xlsx("~/nas/99.reference/DrugCorrection.xlsx")
# num_CancerType = "11.TCGA-STAD"

total_repur = read_xlsx("~/nas/04.Results/drug/depmap/gdsc/total_repurposing_screening_w_target.xlsx")
# Cancerlist %in% 
Cancerlist = Cancerlist[grep(paste0(unique(total_repur$cancer_name) , collapse = "|") , Cancerlist)]
# num_CancerType = "10.TCGA-BLCA"

for (num_CancerType in Cancerlist) {
  
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  Cancername = gsub('TCGA-' , '', CancerType)
  
  # call input
  gc_cellline = readRDS(paste0("~/nas/00.data/filtered_TCGA/", num_CancerType, "/",Cancername,"_cellline_dual_all_log.rds"))
  gc_TCGA = readRDS(paste0("~/nas/04.Results/short_long/", CancerType,"_critical_features_short_long.rds"))
  
  if (num_CancerType == "04.TCGA-CESC") {
    gdsc_sig_cf = data.frame()
    gdsc_specific = data.frame()
  } else {
    gdsc_specific = read.csv(paste0("~/nas/04.Results/drug/depmap/gdsc/", Cancername , "/sig/", Cancername, "_gdsc_screening_sig.csv"))
  }
  
  depmap_sig_genes = read.csv(paste0("~/nas/04.Results/drug/depmap/" , Cancername, "_depmap_sig_w_cf.csv"))
  
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
    select(DepMap.ID, cluster, DRUG_ID,mapped_drug_name  ,DRUG_NAME, PUTATIVE_TARGET, PATHWAY_NAME, AUC, RMSE, LN_IC50,Z_SCORE) %>% 
    arrange(DepMap.ID)
  
  gdsc_w_cluster = left_join(gdsc_w_cluster, criteria_filt, by = c("DRUG_NAME" = "OldName")) %>%
    mutate(DRUG_NAME_new = coalesce(Correction, DRUG_NAME)) %>%
    select(DepMap.ID,cluster , DRUG_ID,DRUG_NAME_new , DRUG_NAME,Correction, everything())
  
  nih_each = nih_drug_filted %>% filter(Cancer_type == Cancername)
  
  nih_each_edit = left_join(nih_each, criteria_filt, by = c("Drug_Name" = "OldName")) %>%
    mutate(main_name_new = coalesce(Correction, Drug_Name ))
  
  tmp_anti = anticancer_drug_filted %>% 
    select(Product, Indications, Targets, 
           colnames(anticancer_drug_filted)[colnames(anticancer_drug_filted) == Cancername])
  
  tmp_anti_filt = tmp_anti[which(tmp_anti[,Cancername] == "Y"),]
  
  anticancer_drug_filted_edit = left_join(tmp_anti_filt, criteria_filt, by = c("Product" = "OldName")) %>%
    mutate(main_name_new = coalesce(Correction, Product )) %>% 
    select(main_name_new , Product,any_of(Cancername), Targets)
  
  gdsc_w_cluster_filt = gdsc_w_cluster %>% 
    filter(DRUG_NAME_new %in% 
             unique(c(unique(nih_each_edit$main_name_new), 
                      # unique(cli_drug_filt_edit$main_name_merge), 
                      unique(anticancer_drug_filted_edit$main_name_new))))
  
  # original sig from depmap
  
  depmap_sig_cf = c(depmap_sig_genes %>%
                      na.omit(pathwaylink) %>%
                      pull(pathwaylink) %>%
                      unique() , 
                    depmap_sig_genes %>%
                      na.omit(Pathway) %>%
                      pull(Pathway) %>%
                      unique())
  
  # original sig from gdsc
  
  if (num_CancerType == "04.TCGA-CESC") {
    gdsc_sig_cf = data.frame()
    gdsc_specific = data.frame()
  } else {
    gdsc_sig_cf = gdsc_specific %>% na.omit(critical_features) %>% pull(critical_features) %>% unique()
  }
  
  # for repurposing 
  repur_anti_cf = total_repur %>% filter(cancer_name == Cancername) %>% pull(critical_features_from_anti) 
  repur_anti_cf = repur_anti_cf[!is.na(repur_anti_cf)]
  
  repur_gdsc_cf = total_repur %>% filter(cancer_name == Cancername) %>% pull(critical_features_from_gdsc)
  repur_gdsc_cf = repur_gdsc_cf[!is.na(repur_gdsc_cf)]
  
  total_repur_cf = c(repur_anti_cf , repur_gdsc_cf)
  
  # network for 
  cf_sl = read.csv(paste0(filepath, "04.Results/short_long/ttest_common/",CancerType,"_critical_features_short_long_common.csv"))
  

  merge_depm_repur = unique(c(depmap_sig_cf,total_repur_cf))
  
  common_features = cf_sl %>% 
    filter(classification != "common") %>%
    filter(variable %in% merge_depm_repur) %>%
    pull(variable)
  
  common_link_genes = link_genes_filtered_df %>% 
    filter(Pathway %in% common_features) %>%
    pull(Genes)
  
  common_single_genes = single_genes %>%
    filter(Pathway %in% common_features) %>%
    pull(Genes)
  
  total_common_genes = unique(c(common_link_genes,common_single_genes))
  
  library(clusterProfiler)
  library(org.Hs.eg.db)
  
  total_common_genes_en = AnnotationDbi::select(org.Hs.eg.db, total_common_genes, 'ENTREZID', 'SYMBOL')[
    which(!is.na(AnnotationDbi::select(org.Hs.eg.db, total_common_genes, 'ENTREZID', 'SYMBOL')$ENTREZID)),]
  total_common_genes_en <- data.frame(total_common_genes_en, row.names = NULL)
  
  ego = enrichGO(total_common_genes_en$ENTREZID,
           OrgDb = org.Hs.eg.db,
           ont = "BP",
           pAdjustMethod = "BH",
           pvalueCutoff = 0.01,
           qvalueCutoff = 0.05,
           readable = T)
  
  egox = setReadable(ego,'org.Hs.eg.db', 'ENTREZID')
  
  
  cnet = cnetplot(egox,
           categorySize="pvalue",
           colorEdge = TRUE)
  setwd("~/nas/04.Results/GOenrichment_test/BP/")
  ggsave(cnet , 
         filename = paste0(Cancername,"_repur_depmap_criticalfeatures_cnetplot.svg"),
         height = 10,
         width = 17)

  # lihc = readRDS(paste0("~/nas/04.Results/short_long/", CancerType,"_critical_features_short_long_with_drug.rds"))
}



