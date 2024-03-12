
library(readxl)
library(tidyverse)
library(svglite)
library(ggplot2)
library(ggsignif)
library(openxlsx)
library(ggpubr)
library(ggsignif)

filepath = "/home/seokwon/nas/"
ref_path = paste0(filepath, "99.reference/")

Cancerlist = dir(paste0(filepath, "/00.data/filtered_TCGA/"))
Cancerlist = Cancerlist[c(-11,-12)]

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
num_CancerType = "19.TCGA-LIHC"

for (num_CancerType in Cancerlist) {
  
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  Cancername = gsub('TCGA-' , '', CancerType)
  
  # call input
  gc_cellline = readRDS(paste0("~/nas/00.data/filtered_TCGA/", num_CancerType, "/",Cancername,"_cellline_dual_all_log.rds"))
  gc_TCGA = readRDS(paste0("~/nas/04.Results/short_long/", CancerType,"_critical_features_short_long.rds"))
  cli_drug = read_xlsx(paste0(ref_path,"TCGA_clinical_drug/",Cancername, "_drug_info_update.xlsx")) # original
  
  # cli_drug = cli_drug %>% select(-X)
  cli_drug_filt = cli_drug %>% 
    slice(-1,-2) %>% 
    filter(!pharmaceutical_therapy_drug_name %in% c("[Not Available]", "Unknown")) %>%
    select(-bcr_patient_uuid, -bcr_drug_uuid)

  cli_drug_filt_edit = left_join(cli_drug_filt, criteria_filt, by = c("pharmaceutical_therapy_drug_name" = "OldName")) %>%
    mutate(main_name_merge = coalesce(Correction, pharmaceutical_therapy_drug_name)) %>%
    select(main_name_merge, Correction,pharmaceutical_therapy_drug_name, everything())
  
  filt_cancer_cell = meta_cell %>% 
    filter(DepMap.ID %in% rownames(gc_cellline)) %>%
    filter(DepMap == T & GDSC == T)
  
  gdsc_each = gdsc %>% filter(COSMIC.ID %in% filt_cancer_cell$COSMIC.ID)
  
  gc_cellline_filt_df = readRDS(paste0("~/nas/04.Results/drug/depmap/gdsc/",CancerType, "_DM_sl_cluster.rds"))
  
  # gdsc_each
  # filt_cancer_cell
  
  gc_filtered_cellline_df = gc_cellline_filt_df %>% filter(rownames(.) %in% filt_cancer_cell$DepMap.ID)
  
  gdsc_each_filt = left_join(gdsc_each, filt_cancer_cell %>% select(COSMIC.ID, DepMap.ID), by = "COSMIC.ID")
  
  tmp_df = data.frame(DepMap.ID = rownames(gc_filtered_cellline_df), cluster = gc_filtered_cellline_df$cluster)
  gdsc_w_cluster = left_join(gdsc_each_filt, tmp_df, by = "DepMap.ID")
  
  gdsc_w_cluster = gdsc_w_cluster %>%
    select(DepMap.ID, cluster, DRUG_ID,DRUG_NAME, PUTATIVE_TARGET, PATHWAY_NAME, AUC, RMSE, LN_IC50,Z_SCORE) %>% 
    arrange(DepMap.ID)
  
  gdsc_w_cluster = left_join(gdsc_w_cluster, criteria_filt, by = c("DRUG_NAME" = "OldName")) %>%
    mutate(DRUG_NAME_new = coalesce(Correction, DRUG_NAME)) %>% 
    select(DepMap.ID,cluster , DRUG_ID,DRUG_NAME_new , DRUG_NAME,Correction, everything())
  
  # nih_each = nih_drug_filted %>% filter(Cancer_type != Cancername) # for repurposing. if do not want, change the equal sign
  
  # nih_each_edit = left_join(nih_each, criteria_filt, by = c("Drug_Name" = "OldName")) %>%
  #   mutate(main_name_new = coalesce(Correction, Drug_Name ))
  
  # tmp_anti = anticancer_drug_filted[which(is.na(anticancer_drug_filted[,Cancername])),] 
  # tmp_anti_filt = tmp_anti %>% 
  #   select(Product, Indications, Targets, 
  #          colnames(anticancer_drug_filted)[15:24][colnames(anticancer_drug_filted)[15:24] != Cancername])
  
  # tmp_anti_filt = tmp_anti_filt[rowSums(is.na(tmp_anti_filt[4:12])) != 9,] # for repurposing. if do not want, change the equal sign
  
  # anticancer_drug_filted_edit = left_join(tmp_anti_filt, criteria_filt, by = c("Product" = "OldName")) %>%
  #   mutate(main_name_new = coalesce(Correction, Product )) %>% 
  #   select(main_name_new , Product, Targets, everything())
  
  gdsc_w_cluster_filt = gdsc_w_cluster %>% 
    filter(DRUG_NAME_new %in% unique(cli_drug_filt_edit$main_name_merge))
  

  
  # for fic
  fig_path = paste0(filepath,"/04.Results/drug/depmap/gdsc/", Cancername)
  if(!dir.exists(fig_path)){
    dir.create(fig_path)
    print(paste0("Created folder: ", fig_path))
  } else {
    print(paste0("Folder already exists: ", fig_path))
  }
  setwd(fig_path)
  
  total_long_gd = gdsc_w_cluster_filt %>% filter(cluster == "long")
  total_short_gd = gdsc_w_cluster_filt %>% filter(cluster == "short")
  
  if (nrow(total_long_gd) >1 & nrow(total_short_gd) >1) {
    anno = t.test(total_long_gd$Z_SCORE, total_short_gd$Z_SCORE)$p.value
    
    total_gd = ggplot(gdsc_w_cluster_filt , aes( x = cluster , y = Z_SCORE, fill=cluster)) + 
      geom_violin(color ="black") +
      geom_boxplot(width=0.1, color = "black" , fill="white")+
      # scale_color_manual(values="black","black") + 
      scale_fill_manual(values=c("#4DAF4A", "#E41A1C")) +
      
      geom_signif(
        annotation = paste0(formatC(anno, digits = 3)),
        map_signif_level = TRUE, 
        comparisons = list(c("long", "short")),
        # y_position = 4.05, xmin = 1, xmax = 3,
        # tip_length = c(0.22, 0.02),
      ) +
      ggtitle(Cancername) +
      # stat_compare_means(label.y = 10) +
      theme_minimal()
    
    ggsave(filename = paste0(CancerType,"_total_treated_drug_IC50_zscore.svg"), total_gd)
  } else {
    
    next
  }
}


