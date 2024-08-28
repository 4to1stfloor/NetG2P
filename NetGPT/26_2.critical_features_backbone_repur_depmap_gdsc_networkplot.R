
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
  
  edges = data.frame(matrix(ncol = 4))
  colnames(edges) = c("from","to","group","variable")
  # critical_features = "P16"
  
  for (critical_features in cf_sl %>% pull(variable)) {
    # nodes = unique(c(nodes , paste0("P",strsplit(critical_features, "P")[[1]][2]) , paste0("P",strsplit(critical_features, "P")[[1]][3])))
    
    if (is.na(strsplit(critical_features, "P")[[1]][3] )) {
      edges[critical_features,"from"] = paste0("P",strsplit(critical_features, "P")[[1]][2])
      edges[critical_features,"to"] = paste0("P",strsplit(critical_features, "P")[[1]][2])
      edges[critical_features,"group"] = cf_sl %>% filter(variable == critical_features) %>% select(classification) %>% pull()
      edges[critical_features,"variable"] = critical_features
      
    } else {
      edges[critical_features,"from"] = paste0("P",strsplit(critical_features, "P")[[1]][2]) 
      edges[critical_features,"to"] = paste0("P",strsplit(critical_features, "P")[[1]][3])
      edges[critical_features,"group"] = cf_sl %>% filter(variable == critical_features) %>% select(classification) %>% pull()
      edges[critical_features,"variable"] = critical_features
      
    }
    
  }
  
  edges = edges[-1,]
  # spe_features_drug
  rownames(edges) = NULL

  edges = left_join(edges, cf_sl , by = "variable") %>% select(from, to, group, variable , classification)

  set.seed(1234)
  backbone = edges %>%
    as_tbl_graph() %N>% 
    left_join(., cf_sl %>% select(variable, classification), by = c("name" = "variable")) %N>%
    mutate(classification = ifelse(is.na(classification), "", classification)) %N>%
    mutate(color = case_when(classification == "short"  ~ "short",
                             classification == "long"  ~ "long",
                             classification == "common"  ~ "common",
                             .default = "")) %E>%
    mutate(color = case_when(classification == "short"  ~ "short",
                             classification == "long"  ~ "long",
                             classification == "common"  ~ "common",
                             .default = "")) %E>%
    mutate(edge_weight = case_when(color == "short" ~ 1,
                                   color == "long" ~ 1,
                                   .default = 0.3)) %N>% 
    mutate(node_alpha = case_when(classification %in% c("short", "long") ~ 1,
                                   .default = 0.5)) %E>%
    mutate(edge_alpha = case_when(classification %in% c("short", "long") ~ 1,
                                  .default = 0.5)) %>% 
    ggraph(layout = "nicely") +
    geom_edge_link(aes(color = color,
                       edge_alpha = edge_alpha,
                       edge_width = edge_weight)) +
    scale_edge_width(range = c(0.5,1.5)) +   
    scale_edge_alpha(range = c(0.3,1)) +
    geom_node_point(aes(color = color,
                        alpha = node_alpha), # 노드 색깔
                    size = 5) +               # 노드 크기
    geom_node_text(aes(label = name),         # 텍스트 표시
                   repel = T,                 # 노드밖 표시
                   size = 5) +  
    # scale_color_manual(values = c("short" = "blue", "common" = "red", "long" = "green")) +
    scale_color_manual(values = c("short" = "#E41A1C", 
                                  "long" = "#4DAF4A",
                                  "common"= "grey")) +
    # scale_color_manual(values = c("short" = "#E41A1C", "long" = "#4DAF4A")) +
    scale_edge_colour_manual(values =c("short" = "#E41A1C", 
                                       "long" = "#4DAF4A",
                                       "common"= "grey")) +
    # geom_segment(aes(x = -3, y = -1.5, xend = 0.1, yend = 0.1),
    #                           arrow = arrow(length = unit(0.1, "cm"))) +
    theme_graph() 
  
  
  set.seed(1234)
  drug_net = edges %>%
    as_tbl_graph() %N>% 
    left_join(., cf_sl %>% select(variable, classification), by = c("name" = "variable")) %N>%
    mutate(classification = ifelse(is.na(classification), "", classification)) %N>% 
    mutate(sp_cf = case_when(name %in% depmap_sig_cf ~ 1,
                             name %in% gdsc_sig_cf ~ 2,
                             name %in% total_repur_cf ~ 3,
                             .default = 0 )) %E>%
    mutate(sp_cf = case_when(variable  %in% depmap_sig_cf ~ 1,
                             variable  %in% gdsc_sig_cf ~ 2,
                             variable  %in% total_repur_cf ~ 3,
                             .default = 0 )) %N>%
    mutate(color = case_when(sp_cf == 1  ~ "depmap",
                             sp_cf == 2  ~ "gdsc",
                             sp_cf == 3  ~ "repur",
                             .default = "unspe")) %E>%
    mutate(color = case_when(sp_cf == 1  ~ "depmap",
                             sp_cf == 2  ~ "gdsc",
                             sp_cf == 3  ~ "repur",
                             .default = "unspe")) %E>%
    mutate(edge_weight = case_when(color == "unspe" ~ 0.3,
                                   .default = 1))  %N>% 
    mutate(node_alpha = case_when(classification %in% c("short", "long") ~ 1,
                                  .default = 0.5)) %E>%
    mutate(edge_alpha = case_when(color == "unspe" ~ 0.5,
                                  .default = 1)) %>% 
    ggraph(layout = "nicely") +
    geom_edge_link(aes(color = color,
                       edge_alpha = edge_alpha,
                       edge_width = edge_weight)) +             # 엣지 명암
    scale_edge_width(range = c(0.5,1.5)) +   
    scale_edge_alpha(range = c(0.3,1)) +
    geom_node_point(aes(color = color,
                        alpha = node_alpha),     # 노드 색깔
                    size = 5) +               # 노드 크기
    geom_node_text(aes(label = name),         # 텍스트 표시
                   repel = T,                 # 노드밖 표시
                   size = 5) +  
    # scale_color_manual(values = c("short" = "blue", "common" = "red", "long" = "green")) +
    scale_color_manual(values = c("unspe" = "grey50", 
                                  "repur" = "#925E9FFF",
                                  "depmap"= "#E41A1C",
                                  "gdsc" = "#AD002AFF")) +
    # scale_color_manual(values = c("short" = "#E41A1C", "long" = "#4DAF4A")) +
    scale_edge_colour_manual(values =c("unspe" = "grey50", 
                                       "repur" = "#925E9FFF",
                                       "depmap"= "#E41A1C",
                                       "gdsc" = "#AD002AFF")) +
    # geom_segment(aes(x = -3, y = -1.5, xend = 0.1, yend = 0.1),
    #                           arrow = arrow(length = unit(0.1, "cm"))) +
    theme_graph()           
  
  setwd("~/nas/04.Results/drug/")
  ggsave(backbone , filename = paste0(Cancername,"_backbone_w_slc.svg"))
  ggsave(drug_net , filename = paste0(Cancername,"_spe_repur_network.svg"))
  
  # lihc = readRDS(paste0("~/nas/04.Results/short_long/", CancerType,"_critical_features_short_long_with_drug.rds"))
}
