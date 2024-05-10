library(readxl)
library(tidyverse)
library(svglite)
library(ggplot2)
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

num_CancerType = "10.TCGA-BLCA"


sig_drugs = c()
annotation_col = data.frame()
for (num_CancerType in Cancerlist) {
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  Cancername = gsub('TCGA-' , '', CancerType)
  
  # call input
  gc_cellline = readRDS(paste0("~/nas/00.data/filtered_TCGA/", num_CancerType, "/",Cancername,"_cellline_dual_all_log.rds"))
  gc_TCGA = readRDS(paste0("~/nas/04.Results/short_long/", CancerType,"_critical_features_short_long.rds"))
  cli_drug = read_xlsx(paste0(ref_path,"TCGA_clinical_drug/",Cancername, "_drug_info_update.xlsx"))  # original

  # cli_drug = cli_drug %>% select(-X)
  cli_drug_filt = cli_drug %>% 
    filter(!pharmaceutical_therapy_drug_name %in% c("[Not Available]", "Unknown")) %>%
    select(-bcr_patient_uuid, -bcr_drug_uuid)
  
  cli_drug_filt_edit = left_join(cli_drug_filt, criteria_filt, by = c("pharmaceutical_therapy_drug_name" = "OldName")) %>%
    mutate(main_name_new = coalesce(Correction, pharmaceutical_therapy_drug_name)) %>%
    select(bcr_patient_barcode,main_name_new, Correction,pharmaceutical_therapy_drug_name, everything())
  
  filt_cancer_cell = meta_cell %>% 
    filter(DepMap.ID %in% rownames(gc_cellline)) %>%
    filter(DepMap == T & GDSC == T)
  
  setwd("~/nas/04.Results/drug/depmap/gdsc/")
  gdsc_each = gdsc %>% filter(COSMIC.ID %in% filt_cancer_cell$COSMIC.ID)
  
  gc_cellline_filt_df = readRDS(paste0(CancerType, "_DM_sl_cluster.rds"))

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
             unique(c(unique(nih_each_edit$main_name_new), unique(cli_drug_filt_edit$main_name_new), unique(anticancer_drug_filted_edit$main_name_new))))
  
  # drug_name = "Lapatinib"
  anno_drug = c()
  for (drug_name in unique(gdsc_w_cluster_filt$DRUG_NAME_new )) {
    # n = n+1
    # print(drug_name)
    tmp_for_drug = gdsc_w_cluster_filt %>% filter(DRUG_NAME_new == drug_name)
    
    tmp_anti_long = tmp_for_drug %>% filter(cluster == "long")
    tmp_anti_short = tmp_for_drug %>% filter(cluster == "short")
    
    if (nrow(tmp_anti_long) < 2 | nrow(tmp_anti_short) < 2) {
      next
    } else {
      anno_tmp = t.test(tmp_anti_long$Z_SCORE, tmp_anti_short$Z_SCORE)$p.value
      
      if (anno_tmp < 0.05 ) {
        sig_drugs = c(sig_drugs, drug_name)
        anno_drug = c(anno_drug, drug_name)
      } else {
        next
      }
      # remove(anno_tmp)
    }
  }
  
  if (length(anno_drug) == 0 ) {
    next
  } else {
    tmp_anno = data.frame(cancertype = Cancername, drugs = anno_drug)
    annotation_col = rbind(annotation_col, tmp_anno)
  }
  
  if (length(sig_drugs) == 0 )  {
    next
  } 
}

total_sig_df = data.frame()
total_sig_anno = data.frame()
# num_CancerType = "11.TCGA-STAD"

for (num_CancerType in Cancerlist) {
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  Cancername = gsub('TCGA-' , '', CancerType)
  
  # call input
  gc_cellline = readRDS(paste0("~/nas/00.data/filtered_TCGA/", num_CancerType, "/",Cancername,"_cellline_dual_all_log.rds"))
  gc_TCGA = readRDS(paste0("~/nas/04.Results/short_long/", CancerType,"_critical_features_short_long.rds"))
  cli_drug = read_xlsx(paste0(ref_path,"TCGA_clinical_drug/",Cancername, "_drug_info_update.xlsx"))  # original
  
  # cli_drug = cli_drug %>% select(-X)
  cli_drug_filt = cli_drug %>% 
    filter(!pharmaceutical_therapy_drug_name %in% c("[Not Available]", "Unknown")) %>%
    select(-bcr_patient_uuid, -bcr_drug_uuid)
  
  cli_drug_filt_edit = left_join(cli_drug_filt, criteria_filt, by = c("pharmaceutical_therapy_drug_name" = "OldName")) %>%
    mutate(main_name_new = coalesce(Correction, pharmaceutical_therapy_drug_name)) %>%
    select(bcr_patient_barcode,main_name_new, Correction,pharmaceutical_therapy_drug_name, everything())
  
  filt_cancer_cell = meta_cell %>% 
    filter(DepMap.ID %in% rownames(gc_cellline)) %>%
    filter(DepMap == T & GDSC == T)
  
  setwd("~/nas/04.Results/drug/depmap/gdsc/")
  gdsc_each = gdsc %>% filter(COSMIC.ID %in% filt_cancer_cell$COSMIC.ID)
  
  gc_cellline_filt_df = readRDS(paste0(CancerType, "_DM_sl_cluster.rds"))
  
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
             unique(c(unique(nih_each_edit$main_name_new), unique(cli_drug_filt_edit$main_name_new), unique(anticancer_drug_filted_edit$main_name_new))))
  
  # drug_name = "Lapatinib"
  gdsc_w_sig_cf_filt = gdsc_w_cluster_filt %>% filter(DRUG_NAME_new %in% sig_drugs)
  
  long_gdsc_cf_filt = gdsc_w_sig_cf_filt %>% filter(cluster == "long")
  short_gdsc_cf_filt = gdsc_w_sig_cf_filt %>% filter(cluster == "short")
  
  long_gdsc_cf_filt = long_gdsc_cf_filt %>% 
    select(DepMap.ID, DRUG_NAME_new,Z_SCORE) %>% 
    pivot_wider(names_from = DRUG_NAME_new, 
                values_from = Z_SCORE,
                values_fn = mean,
                values_fill = 0) %>%
    as.data.frame()
  
  rownames(long_gdsc_cf_filt) = long_gdsc_cf_filt$DepMap.ID
  long_gdsc_cf_filt$DepMap.ID = NULL
  
  short_gdsc_cf_filt = short_gdsc_cf_filt %>% 
    select(DepMap.ID, DRUG_NAME_new,Z_SCORE) %>% 
    pivot_wider(names_from = DRUG_NAME_new, 
                values_from = Z_SCORE,
                values_fn = mean,
                values_fill = 0) %>%
    as.data.frame()
  
  rownames(short_gdsc_cf_filt) = short_gdsc_cf_filt$DepMap.ID
  short_gdsc_cf_filt$DepMap.ID = NULL
  
  if (nrow(long_gdsc_cf_filt) < 2 | nrow(short_gdsc_cf_filt) < 2) {
    next
  } else {
    long_sum = long_gdsc_cf_filt %>% colMeans() %>% as.data.frame()
    short_sum = short_gdsc_cf_filt %>% colMeans() %>% as.data.frame()
    short_sum = short_sum[rownames(long_sum), , drop =F]
    delta_sum = long_sum - short_sum
    colnames(delta_sum) = c(paste0(Cancername , "_mean_delta"))
    
    posi_delta_sum = delta_sum[which(delta_sum > 0 ), , drop = F]
    neg_delta_sum = delta_sum[which(delta_sum < 0 ), , drop = F]
    
    if (nrow(posi_delta_sum) != 0) {
      posi_delta_sum$direction = ""
      posi_long_sum = long_sum[rownames(posi_delta_sum), , drop =F]
      posi_short_sum = short_sum[rownames(posi_delta_sum), , drop =F]
      
      if (nrow(posi_delta_sum[which(posi_long_sum < 0 & posi_short_sum < 0 ),]) != 0) {
        posi_delta_sum[rownames(posi_delta_sum[which(posi_long_sum < 0 & posi_short_sum < 0 ), , drop = F]),]$direction = "short_sensitive"
      } 
      
      if (nrow(posi_delta_sum[which(posi_long_sum < 0 & posi_short_sum > 0 ),]) != 0) {
        # posi_delta_sum[rownames(posi_delta_sum[which(posi_long_sum < 0 & posi_short_sum > 0 ), , drop = F]),]$direction = "long_sensitive"
      } 
      
      if (nrow(posi_delta_sum[which(posi_long_sum > 0 & posi_short_sum < 0 ),]) != 0) {
        posi_delta_sum[rownames(posi_delta_sum[which(posi_long_sum > 0 & posi_short_sum < 0 ), , drop = F]),]$direction = "short_sensitive"
      } 
      
      if (nrow(posi_delta_sum[which(posi_long_sum > 0 & posi_short_sum > 0 ),]) != 0) {
        posi_delta_sum[rownames(posi_delta_sum[which(posi_long_sum > 0 & posi_short_sum > 0 ), , drop = F]),]$direction = "short_sensitive"
      }
    }
   
    if (nrow(neg_delta_sum) != 0) {
      neg_delta_sum$direction = ""
      neg_long_sum = long_sum[rownames(neg_delta_sum), , drop =F]
      neg_short_sum = short_sum[rownames(neg_delta_sum), , drop =F]
      
      if (nrow(neg_delta_sum[which(neg_long_sum < 0 & neg_short_sum < 0 ),]) != 0) {
        neg_delta_sum[rownames(neg_delta_sum[which(neg_long_sum < 0 & neg_short_sum < 0 ), , drop = F]),]$direction = "long_sensitive"
      } 
      
      if (nrow(neg_delta_sum[which(neg_long_sum < 0 & neg_short_sum > 0 ),]) != 0) {
        neg_delta_sum[rownames(neg_delta_sum[which(neg_long_sum < 0 & neg_short_sum > 0 ), , drop = F]),]$direction = "long_sensitive"
      } 
      
      if (nrow(neg_delta_sum[which(neg_long_sum > 0 & neg_short_sum < 0 ),]) != 0) {
        # neg_delta_sum[rownames(neg_delta_sum[which(neg_long_sum > 0 & neg_short_sum < 0 ), , drop = F]),]$direction = "short"
      } 
      
      if (nrow(neg_delta_sum[which(neg_long_sum > 0 & neg_short_sum > 0 ),]) != 0) {
        
        neg_delta_sum[rownames(neg_delta_sum[which(neg_long_sum > 0 & neg_short_sum > 0 ), , drop = F]),]$direction = "long_sensitive"
      }
      
    }
    
    delta_sum = rbind(posi_delta_sum,neg_delta_sum)
    
    tmp_delta = as.data.frame(t(delta_sum %>% select(any_of(paste0(Cancername , "_mean_delta")))))
    # tmp_anno = delta_sum %>% select(direction) %>% t() %>% as.data.frame()
    tmp_anno = delta_sum %>% select(direction)
    tmp_anno$drugs = rownames(tmp_anno)
    rownames(tmp_anno) = NULL
    tmp_anno$cancertype = Cancername
    tmp_anno = tmp_anno %>% select( cancertype, drugs, direction)
    
    # tmp_delta = delta_sum %>% t() %>% as.data.frame()
    
    total_sig_df = bind_rows(total_sig_df, tmp_delta)
    total_sig_anno = rbind(total_sig_anno, tmp_anno)
    
  }
  
}

total_sig_rev = total_sig_df # delta 에서는 금지
# total_sig_rev = -total_sig_matrix # delta 에서는 금지
# total_rev_filt = total_rev + (-min(total_rev , na.rm = T))
total_sig_rev[is.na(total_sig_rev)] = 0

total_sig_rev_back = total_sig_rev

total_sig_rev_back[total_sig_rev_back > 0] = (total_sig_rev_back[total_sig_rev_back>0] - min(total_sig_rev_back[total_sig_rev_back>0])) / 
  (max(total_sig_rev_back[total_sig_rev_back>0]) - min(total_sig_rev_back[total_sig_rev_back>0])) *
  (1 -  min(total_sig_rev_back[total_sig_rev_back>0])) + min(total_sig_rev_back[total_sig_rev_back>0])

total_sig_rev_back[total_sig_rev_back < 0] = (total_sig_rev_back[total_sig_rev_back < 0] - min(total_sig_rev_back[total_sig_rev_back < 0])) / 
  (max(total_sig_rev_back[total_sig_rev_back < 0]) - min(total_sig_rev_back[total_sig_rev_back < 0])) *
  (max(total_sig_rev_back[total_sig_rev_back < 0]) - -1) + -1

annotation_col = annotation_col %>%
  left_join(total_sig_anno, by = c("cancertype", "drugs"))

other_df = data.frame(cancertype = "multiple_specific", 
                      drugs = annotation_col[which(duplicated(annotation_col$drugs)),]$drugs)

multi_drugs = other_df$drugs
multi_only = total_sig_rev %>% 
  select(any_of(multi_drugs))

multi_only = multi_only[grep(paste0(annotation_col %>% filter(drugs == other_df$drugs) %>% pull(cancertype), collapse = "|"),
                             rownames(multi_only)), , drop=F]

multi_only$direction = annotation_col %>% filter(drugs == other_df$drugs) %>% pull(direction)

# mg = "Fluorouracil"
# remove(mg)
other_sort_df = data.frame()
for (mg in colnames(multi_only %>% select(-direction))) {
  tmp_m = multi_only %>% select(any_of(mg), direction) 
  tmp_other = data.frame(cancertype = unlist(strsplit(rownames(tmp_m)[which(tmp_m$direction == "long_sensitive")], "_"))[1] ,
                         drugs = mg)
  
  tmp_other = cbind(tmp_other , total_sig_anno %>% 
                      filter(cancertype == tmp_other$cancertype) %>% 
                      filter(drugs == tmp_other$drugs) %>% 
                      select(direction))
  other_sort_df = rbind(other_sort_df,tmp_other) 
}

filtered_annotation_col = annotation_col %>%
  filter(!(drugs %in% other_sort_df$drugs))
final_annotation_col <- rbind(filtered_annotation_col, other_sort_df)

rownames(final_annotation_col) = final_annotation_col$drugs

final_annotation_col = final_annotation_col %>% select(-drugs)
final_annotation_col = final_annotation_col[colnames(total_sig_rev_back), , drop = FALSE]

## re-order
# cancername = "CESC"
ordered_drugs = c()
for (cancername in unique(sapply(strsplit(rownames(total_sig_rev_back), "_"), "[", 1))) {
  if (cancername %in% final_annotation_col$cancertype) {
    tmp_arrange = total_sig_rev_back %>% 
      select(final_annotation_col %>% 
               filter(cancertype == cancername) %>% 
               mutate(direction = factor(direction, levels = c("long_sensitive", "short_sensitive"))) %>%
               arrange(direction) %>% 
               rownames(.)) %>% 
      filter(grepl(cancername, rownames(.)))
    
    tmp_ordered = tmp_arrange %>% gather(key = "gene", value = "value") %>%
      arrange(desc(value)) %>% pull(gene)
    # tmp_ordered = colnames(tmp_arrange)
    
    ordered_drugs = c(ordered_drugs, tmp_ordered)
  }
  
}

if (sum(colnames(total_sig_rev_back) %in% ordered_drugs) == ncol(total_sig_rev_back) ) {
  total_sig_rev_back = total_sig_rev_back[,ordered_drugs]
} else {
  print("error")
}
final_annotation_col = final_annotation_col[colnames(total_sig_rev_back),]

ordered_drugs = c()
for (cancername in unique(sapply(strsplit(rownames(total_sig_rev_back), "_"), "[", 1))) {
  if (cancername %in% final_annotation_col$cancertype) {
    tmp_arrange = total_sig_rev_back %>% 
      select(final_annotation_col %>% 
               filter(cancertype == cancername) %>% 
               mutate(direction = factor(direction, levels = c("long_sensitive", "short_sensitive"))) %>%
               arrange(direction) %>% 
               rownames(.)) %>% 
      filter(grepl(cancername, rownames(.)))
    
    tmp_ordered = tmp_arrange %>% gather(key = "gene", value = "value") %>%
      arrange(desc(value)) %>% pull(gene)
    # tmp_ordered = colnames(tmp_arrange)
    
    ordered_drugs = c(ordered_drugs, tmp_ordered)
  }
}

if (sum(colnames(total_sig_rev_back) %in% ordered_drugs) == ncol(total_sig_rev_back) ) {
  total_sig_rev_back = total_sig_rev_back[,ordered_drugs]
} else {
  print("error")
}
final_annotation_col = final_annotation_col[colnames(total_sig_rev_back),]
# saveRDS(annotation_col, "~/nas/04.Results/drug/depmap/delta_cellline_heatmap_anno_col.rds")

# 
# cancer_colors = c(
#   "LGG" = "#00A087FF", # LGG
#   # "OV" = "#3C5488FF", # OV
#   # "BRCA" = "#4DBBD5FF", # BRCA
#   # "BLCA" = "#7E6148FF", # BLCA
#   # "LIHC" = "#8491B4FF", # LIHC
#   # "STAD" = "#91D1C2FF", # STAD
#   # "LUSC" = "#B09C85FF", # LUSC
#   "CESC" = "#DC0000FF", # CESC
#   # "UCEC" = "#E64B35FF", # UCEC
#   "LUAD" = "#F39B7FFF" # LUAD
# )

cancer_colors = c(

  "OV" = "#3C5488FF", # OV
  "BRCA" = "#4DBBD5FF", # BRCA
  "BLCA" = "#7E6148FF", # BLCA
  "LIHC" = "#8491B4FF", # LIHC
  "STAD" = "#91D1C2FF", # STAD
  "LUSC" = "#B09C85FF", # LUSC
  "UCEC" = "#E64B35FF" # UCEC

)

spe_colors = c(

  "long_sensitive"  = "#009E73", # effective both direction
  "short_sensitive" = "#FF3824" # effective both direction

)

row_colors_cluster <- c(
  # "LGG" = "#00A087FF", # LGG
  "OV" = "#3C5488FF", # OV
  "BRCA" = "#4DBBD5FF", # BRCA
  "BLCA" = "#7E6148FF", # BLCA
  "LIHC" = "#8491B4FF", # LIHC
  "STAD" = "#91D1C2FF", # STAD
  "LUSC" = "#B09C85FF", # LUSC
  "CESC" = "#DC0000FF", # CESC
  "UCEC" = "#E64B35FF", # UCEC
  "LUAD" = "#F39B7FFF"  # LUAD
)
# row_colors_sl <- c("short" = "red", "long" = "#009E73")

ann_colors_sl = list(direction = spe_colors, cancertype = cancer_colors ,cluster = row_colors_cluster)

final_total_sig_rev = total_sig_rev_back[grep(paste0(final_annotation_col$cancertype , collapse = "|"), rownames(total_sig_rev_back)), , drop = F]

annotation_row = data.frame(types = rownames(final_total_sig_rev))
annotation_row  = annotation_row %>% 
  mutate(cluster = sapply(strsplit(annotation_row$types, "_"), "[", 1)) 
rownames(annotation_row) = annotation_row$types

annotation_row = annotation_row %>% select(-types) %>% select(cluster)

svglite(filename = "gdsc_total_sig.svg" ,width = 10 , height = 5)

tmp_heatmap = ComplexHeatmap::pheatmap(final_total_sig_rev %>% as.matrix(),
                                       cluster_rows = F,
                                       cluster_cols = F,
                                       breaks = c(-1,max(final_total_sig_rev[final_total_sig_rev < 0]) , 0 , min(final_total_sig_rev[final_total_sig_rev > 0]), 1) ,
                                       column_split = factor(final_annotation_col$cancertype, levels = c(unique(final_annotation_col$cancertype))) ,
                                       row_split = factor(annotation_row$cluster, levels = c(unique(annotation_row$cluster))),
                                       annotation_colors = ann_colors_sl,
                                       annotation_col = final_annotation_col %>% select(direction),
                                       annotation_row = annotation_row,
                                       legend = T,
                                       # annotation_legend = T,
                                       # annotation_names_col = F,
                                       # annotation_names_row = F,
                                       show_colnames = T,
                                       show_rownames = F,
                                       cluster_column_slices = FALSE,
                                       color = c( colorRampPalette(c("#387eb8", "#F9FEFE", "#e21e26"))(length(unique(unlist(final_total_sig_rev))))),
                                       scale = "none",
                                       border_color = NA
) 

print(tmp_heatmap)

dev.off()



