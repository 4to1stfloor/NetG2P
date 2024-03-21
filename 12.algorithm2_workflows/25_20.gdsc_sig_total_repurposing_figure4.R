library(readxl)
library(tidyverse)
library(svglite)
library(ggplot2)
library(ggsignif)
library(openxlsx)
library(ggpubr)
library(rstatix)

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
total_repur = read_xlsx("~/nas/04.Results/drug/depmap/gdsc/total_repurposing_screening_w_target.xlsx")
# num_CancerType = "11.TCGA-STAD"
# unique(total_repur$cancer_name)
total_gd_select_arrange = data.frame()
Cancerlist = Cancerlist[c(2,3,7)]
for (num_CancerType in Cancerlist) {
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  Cancername = gsub('TCGA-' , '', CancerType)
  
  # call input
  gc_cellline = readRDS(paste0("~/nas/00.data/filtered_TCGA/", num_CancerType, "/",Cancername,"_cellline_dual_all_log.rds"))
  
  filt_cancer_cell = meta_cell %>% 
    filter(DepMap.ID %in% rownames(gc_cellline)) %>%
    filter(DepMap == T & GDSC == T)
  
  # setwd("~/nas/04.Results/drug/depmap/gdsc/")
  gdsc_each = gdsc %>% filter(COSMIC.ID %in% filt_cancer_cell$COSMIC.ID)
  
  gc_cellline_filt_df = readRDS(paste0("~/nas/04.Results/drug/depmap/gdsc/",CancerType, "_DM_sl_cluster.rds"))
  
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
    mutate(test_drug = coalesce(Correction, mapped_drug_name)) %>% 
    select(DepMap.ID,cluster , DRUG_ID,DRUG_NAME_new , test_drug , DRUG_NAME,Correction, everything())
  
  nih_each = nih_drug_filted %>% filter(Cancer_type != Cancername) # for repurposing. if do not want, change the equal sign
  
  nih_each_edit = left_join(nih_each, criteria_filt, by = c("Drug_Name" = "OldName")) %>%
    mutate(main_name_new = coalesce(Correction, Drug_Name ))
  
  tmp_anti = anticancer_drug_filted[which(is.na(anticancer_drug_filted[,Cancername])),] 
  tmp_anti_filt = tmp_anti %>% 
    select(Product, Indications, Targets, 
           colnames(anticancer_drug_filted)[15:24][colnames(anticancer_drug_filted)[15:24] != Cancername])
  
  tmp_anti_filt = tmp_anti_filt[rowSums(is.na(tmp_anti_filt[4:12])) != 9,] # for repurposing. if do not want, change the equal sign
  
  anticancer_drug_filted_edit = left_join(tmp_anti_filt, criteria_filt, by = c("Product" = "OldName")) %>%
    mutate(main_name_new = coalesce(Correction, Product )) %>% 
    select(main_name_new , Product, Targets, everything())
  
  repur_drug = total_repur %>% filter(cancer_name == Cancername) %>% pull(drug_name)
  
  gdsc_w_cluster_filt = gdsc_w_cluster %>% filter(DRUG_NAME_new %in% repur_drug)
  
  if (all.equal(unique(gdsc_w_cluster_filt$DRUG_NAME_new), repur_drug)) {
    print("right!")
  }else {
    print(Cancername , "Wrong")
  }
  
  gd_select = gdsc_w_cluster_filt %>% select(DRUG_NAME_new,cluster, Z_SCORE)

  gd_select_arrange = gd_select %>%
    arrange(DRUG_NAME_new)
  gd_select_arrange$cancertype = Cancername
  total_gd_select_arrange = rbind(total_gd_select_arrange ,gd_select_arrange)

}

total_gd_select_arrange

stat_test <- total_gd_select_arrange %>%
  group_by(DRUG_NAME_new, cancertype) %>%
  t_test(Z_SCORE ~ cluster) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")

stat_test <- stat_test %>%
  add_xy_position(fun = "mean_sd", x = "DRUG_NAME_new", dodge = 0.8) 

total_gd_select_arrange$cancertype = factor(total_gd_select_arrange$cancertype , levels = c("STAD", "BLCA", "LUSC"))
total_gd_select_arrange = total_gd_select_arrange %>% arrange(cancertype) 

stat_test$DRUG_NAME_new = factor(stat_test$DRUG_NAME_new , levels = unique(total_gd_select_arrange$DRUG_NAME_new))
stat_test = stat_test %>% arrange(DRUG_NAME_new)
stat_test = stat_test %>% 
  group_by(cancertype) %>% 
  mutate(x = row_number()) %>% 
  ungroup() %>% 
  mutate(xmin = x - 0.2) %>% 
  mutate(xmax = x + 0.2)

stat_test_add_signif = stat_test %>% mutate(p_signif = case_when(p > 0.05 ~ "ns",
                                                                 p < 0.05 & p > 0.01 ~ "*",
                                                                 p < 0.01 & p > 0.001 ~ "**" ,
                                                                 p < 0.001 ~ "***", 
                                                                 .default = NA))

total_box = ggboxplot(total_gd_select_arrange ,
          x = "DRUG_NAME_new",
          y = "Z_SCORE",
          fill = "cluster",
          # palette = "npg",
          palette = c("#4dbbd5", "#e64b35"),
          facet.by = "cancertype",
          short.panel.labs = FALSE,
          panel.labs.background = list(fill = "steelblue", color = "steelblue")
          ) + 
  facet_grid(~ cancertype,scales = "free", space='free') +  
  stat_pvalue_manual(
    stat_test_add_signif,  
    label = "p_signif", 
    tip.length = 0.02,
    bracket.nudge.y = 0.6
  ) +
  ylim(-2.5, 2.5)+ 
  theme(axis.text.x = element_text(angle=45, hjust = 1),
        )+
  xlab("")  +
  guides(fill = "none")


setwd("~/nas/04.Results/drug/depmap/gdsc/")
ggsave(file= "total_repurposing_screening.svg", plot=total_box, width=8, height=4)

# lihc = readRDS(paste0("~/nas/04.Results/short_long/", CancerType,"_critical_features_short_long_with_drug.rds"))



