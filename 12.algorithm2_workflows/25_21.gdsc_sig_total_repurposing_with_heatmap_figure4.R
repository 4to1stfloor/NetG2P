library(readxl)
library(tidyverse)
library(svglite)
library(ggplot2)
library(ggsignif)
library(openxlsx)

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

# num_CancerType = "19.TCGA-LIHC"

total_repur_screening = data.frame()
for (num_CancerType in Cancerlist) {
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  Cancername = gsub('TCGA-' , '', CancerType)
  
  # call input
  gc_cellline = readRDS(paste0("~/nas/00.data/filtered_TCGA/", num_CancerType, "/",Cancername,"_cellline_dual_all_log.rds"))
  gc_TCGA = readRDS(paste0("~/nas/04.Results/short_long/", CancerType,"_critical_features_short_long.rds"))
  
  Cancerlist_edit = Cancerlist[Cancerlist != num_CancerType]
  # repur = "04.TCGA-CESC"
  wo_cli_drug = data.frame()
  
  for (repur in Cancerlist_edit) {
    
    repur_name = gsub('TCGA-' , '', gsub('[.]','',gsub('\\d','', repur)))
    
    tmp_cli_drug = read_xlsx(paste0(ref_path,"TCGA_clinical_drug/",repur_name, "_drug_info_update.xlsx")) # original
    
    tmp_cli_drug_filt = tmp_cli_drug %>% 
      mutate(cancertype = repur_name ) %>%
      select(cancertype , pharmaceutical_therapy_drug_name ) %>% 
      distinct(pharmaceutical_therapy_drug_name, .keep_all = TRUE)
    
    wo_cli_drug = rbind(wo_cli_drug , tmp_cli_drug_filt)
    
  }
  
  cli_drug_filt = wo_cli_drug %>% 
    filter(!pharmaceutical_therapy_drug_name %in% c("[Not Available]", "Unknown")) 
  
  cli_drug_filt_edit = left_join(cli_drug_filt, criteria_filt, by = c("pharmaceutical_therapy_drug_name" = "OldName")) %>%
    mutate(main_name_merge = coalesce(Correction, pharmaceutical_therapy_drug_name)) %>%
    select(cancertype , main_name_merge, Correction,pharmaceutical_therapy_drug_name, everything())
  
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
  
  gdsc_w_cluster_filt = gdsc_w_cluster %>% 
    filter(DRUG_NAME_new %in% 
             unique(c(unique(nih_each_edit$main_name_new), 
                      unique(cli_drug_filt_edit$main_name_merge), 
                      unique(anticancer_drug_filted_edit$main_name_new))))
  
  # drug_name = "Fulvestrant"
  
  for (drug_name in unique(gdsc_w_cluster_filt$DRUG_NAME_new )) {
    # n = n+1
    # print(drug_name)
    tmp_for_drug = gdsc_w_cluster_filt %>% filter(DRUG_NAME_new == drug_name)
    
    tmp_anti_long = tmp_for_drug %>% filter(cluster == "long")
    tmp_anti_short = tmp_for_drug %>% filter(cluster == "short")
    
    if (nrow(tmp_anti_long) < 2 | nrow(tmp_anti_short) < 2) {
      next
    } else {
      # print(anno_tmp)
      anno_tmp = t.test(tmp_anti_long$Z_SCORE, tmp_anti_short$Z_SCORE)$p.value
      nih_from = unique(trimws(unlist(str_split(nih_each_edit %>% 
                                                  filter(main_name_new == drug_name) %>% 
                                                  pull(Cancer_type), ","))))
      TCGA_from = unique(trimws(unlist(str_split(cli_drug_filt_edit %>% 
                                                   filter(main_name_merge == drug_name) %>% 
                                                   pull(cancertype), ","))))
      
      if (mean(tmp_anti_long$Z_SCORE) > mean(tmp_anti_short$Z_SCORE)) {
        
        direction = "short"
        
      } else {
        direction = "long"
      }
      # tmp_anti_other = anticancer_drug_filted_edit %>% 
      #   filter(main_name_new == drug_name) 
      # if (nrow(tmp_anti_other) != 0) {
      #   anti_from = trimws(colnames(tmp_anti_other[,which(tmp_anti_other == "Y")]))
      # }
      
      anti_from = anticancer_drug_filted_edit %>% 
        filter(main_name_new == drug_name) %>%
        slice(1) %>%
        select(5:13) %>% 
        mutate(across(where(~ any(. == "Y" & !is.na(.))), ~.)) %>%
        select(where(~ any(!is.na(.)))) %>% colnames() 
      
      tmp_repur_screen = data.frame(cancer_name = Cancername, 
                                    drug_name = drug_name, 
                                    pvalue = anno_tmp,
                                    direction = direction ,
                                    repur_from_nih = paste0(nih_from, collapse = "; "),
                                    repur_from_TCGA = paste0(TCGA_from, collapse = "; "),
                                    repur_from_anti = paste0(anti_from, collapse = "; "),
                                    target_gene_anti = ifelse(nrow(anticancer_drug_filted_edit %>% 
                                                                     filter(main_name_new == drug_name)) == 0,
                                                              "",anticancer_drug_filted_edit %>% 
                                                                filter(main_name_new == drug_name) %>%
                                                                slice(1) %>%
                                                                pull(Targets)),
                                    target_gene_gdsc = paste0(unique(tmp_for_drug$PUTATIVE_TARGET), collapse = "; ")
      )
      
    }
    total_repur_screening = rbind(total_repur_screening,tmp_repur_screen)
  }
  # lihc = readRDS(paste0("~/nas/04.Results/short_long/", CancerType,"_critical_features_short_long_with_drug.rds"))
}

# lihc_total_repur = total_repur_screening
# lihc_repur_cut = total_repur_screening %>% filter(direction == "short" & pvalue < 0.1)

# write.xlsx(lihc_total_repur , "~/nas/04.Results/drug/depmap/gdsc/lihc_repurposing_screening_w_target.xlsx")
# write.xlsx(lihc_repur_cut , "~/nas/04.Results/drug/depmap/gdsc/lihc_spe_repurposing_w_target.xlsx")

library(readxl)
library(tidyverse)
library(svglite)
library(ggplot2)
library(ggsignif)
library(openxlsx)
library(ggpubr)
library(rstatix)
library(ComplexHeatmap)
library(gridExtra)
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
total_repur = lihc_repur_cut
num_CancerType = "19.TCGA-LIHC"
# unique(total_repur$cancer_name)
total_gd_select_arrange = data.frame()
# Cancerlist = Cancerlist[c(5)]
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
  mutate(xmin = x - 0.2,
         xmax = x + 0.2) 

stat_test_add_signif = stat_test %>% mutate(p_signif = case_when(p > 0.05 ~ "ns",
                                                                 p <= 0.05 & p > 0.01 ~ "*",
                                                                 p <= 0.01 & p > 0.001 ~ "**" ,
                                                                 p <= 0.001 ~ "***", 
                                                                 .default = NA))
stat_test_add_signif$cancertype = factor(stat_test_add_signif$cancertype , levels = c("STAD","BLCA","LUSC"))
stat_test_add_signif = stat_test_add_signif %>% arrange(cancertype)

total_box = ggboxplot(total_gd_select_arrange ,
                      x = "DRUG_NAME_new",
                      y = "Z_SCORE",
                      fill = "cluster",
                      # palette = "npg",
                      palette = c("#4dbbd5", "#e64b35"),
                      facet.by = "cancertype",
                      short.panel.labs = FALSE,
                      panel.labs.background = list(fill = "steelblue", color = "steelblue"),
                      use.label = F,
                      use.labels = F
) + 
  # facet_grid(~ cancertype,scales = "free", space='free') +  
  stat_pvalue_manual(
    stat_test_add_signif,  
    label = "p", 
    tip.length = 0.02,
    bracket.nudge.y = 0.6
  ) +
  ylim(-2.5, 2.5)+ 
  theme(axis.text.x = element_blank())+
  # theme(axis.text.x = element_text(angle=45, hjust = 1))+
  xlab("")  +
  labs(subtitle = "")+
  rremove("x.ticks") +
  guides(fill = "none")

tmp_heat = total_repur %>% select(drug_name , repur_from_nih, repur_from_anti, repur_from_TCGA)

tmp_heat_df = tmp_heat %>%
  pivot_longer(cols = starts_with("repur_from_"), names_to = "variable", values_to = "value") %>%
  separate_rows(value, sep = "; ") %>%
  filter(value != "") %>% 
  select(drug_name, value) %>%
  mutate(value_present = ifelse(!is.na(value), 1, 0))

colnames(tmp_heat_df) = c("drug_name","cancertype", "present")

tmp_heat_filt_df = tmp_heat_df %>%
  mutate(origin = case_when(
    drug_name %in% total_gd_select_arrange$DRUG_NAME_new ~ 
      total_gd_select_arrange$cancertype[match(drug_name, total_gd_select_arrange$DRUG_NAME_new)],
    TRUE ~ NA_character_
  ))

tmp_heat_filt_df$drug_name = factor(tmp_heat_filt_df$drug_name , levels = unique(total_gd_select_arrange$DRUG_NAME_new))
tmp_heat_filt_df = tmp_heat_filt_df %>% arrange(drug_name)
tmp_heat_filt_df$origin = factor(tmp_heat_filt_df$origin , levels = unique(total_gd_select_arrange$cancertype))
tmp_heat_filt_df = tmp_heat_filt_df %>% arrange(origin)

cancer_colors = c(
  "LGG" = "#00A087FF", # LGG
  "OV" = "#3C5488FF", # OV
  "BRCA" = "#4DBBD5FF", # BRCA
  "KIDNEY" = "#A20056FF",
  "PAAD" = "#5F559BFF",
  "PRAD" = "#FFDC91FF",
  "SKCM" = "#808180FF",
  "STAD" = "#91D1C2FF", # STAD
  "LUSC" = "#B09C85FF", # LUSC
  "UVM" = "#008B45FF",
  "UCEC" = "#E64B35FF", # UCEC
  "LUAD" = "#F39B7FFF" # LUAD
)



p2 = tmp_heat_filt_df %>% 
  group_by(cancertype) %>%
  ggplot(., aes(x = drug_name, y = cancertype))+
  geom_tile(aes(fill = cancertype), color = "white") +
  facet_grid(~ origin,scales = "free", space='free') +
  theme_classic() + 
  scale_fill_manual(values = cancer_colors) +
  theme(
    legend.position="none",
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    text = element_text(family = "Helvetica", face = "bold")
  ) + xlab("")  +ylab("")

# theme(axis.text.x = element_text(angle=45, hjust = 1))+
# gA <- ggplotGrob(total_box)
# gB <- ggplotGrob(p2)
# grid::grid.newpage()
# grid::grid.draw(rbind(gA, gB))
# 
# grid.arrange(total_box, p2, nrow=2)
library(cowplot)
total_box_add_heat = plot_grid(total_box ,p2, align = "v",nrow = 2 , rel_heights = c(c(5/8, 3/8)))

setwd("~/nas/04.Results/drug/depmap/gdsc/")
ggsave(file= "lihc_repurposing_screening_with_heat.svg", plot=total_box_add_heat, width=8, height=4.8)

# lihc = readRDS(paste0("~/nas/04.Results/short_long/", CancerType,"_critical_features_short_long_with_drug.rds"))




