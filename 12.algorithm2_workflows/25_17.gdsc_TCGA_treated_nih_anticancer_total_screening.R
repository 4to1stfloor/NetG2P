library(readxl)
library(tidyverse)
library(svglite)
library(ggplot2)
library(ggsignif)

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

# num_CancerType = "30.TCGA-BRCA"

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
             unique(c(unique(nih_each_edit$main_name_new), unique(cli_drug_filt_edit$main_name_new), unique(anticancer_drug_filted_edit$main_name_new))))
  
  fig_path = paste0(filepath,"/04.Results/drug/depmap/gdsc/", Cancername)
  if(!dir.exists(fig_path)){
    dir.create(fig_path)
    print(paste0("Created folder: ", fig_path))
  } else {
    print(paste0("Folder already exists: ", fig_path))
  }
  setwd(fig_path)
  
  total_gdsc_screening = data.frame()
  # drug_name = "Etoposide"
  for (drug_name in unique(gdsc_w_cluster_filt$DRUG_NAME_new )) {
    # n = n+1
    print(drug_name)
    tmp_for_drug = gdsc_w_cluster_filt %>% filter(DRUG_NAME_new == drug_name)
    
    tmp_anti_long = tmp_for_drug %>% filter(cluster == "long")
    tmp_anti_short = tmp_for_drug %>% filter(cluster == "short")
    
    if (nrow(tmp_anti_long) < 2 | nrow(tmp_anti_short) < 2) {
      next
    } else {
      anno_tmp = t.test(tmp_anti_long$Z_SCORE, tmp_anti_short$Z_SCORE)$p.value
      tmp_drug = ggplot(tmp_for_drug , aes( x = cluster , y = Z_SCORE, fill= cluster)) +
        # geom_violin(color ="black") +
        # geom_boxplot(width=0.1, color = "black" , fill="white")+
        geom_boxplot()+
        # scale_color_manual(values="black","black") +
        scale_fill_manual(values=c("#4DAF4A", "#E41A1C")) +

        geom_signif(
          annotation = paste0(formatC(anno_tmp, digits = 3)),
          map_signif_level = TRUE,
          comparisons = list(c("long", "short")),
          # y_position = 3.05,
          # xmin = 1,
          # xmax = 3,
          # tip_length = c(0.22, 0.02),
        ) +
        ggtitle(paste0(drug_name)) +
        # stat_compare_means(label.y = 10) +
        theme_minimal()

      if (anno_tmp < 0.05 ) {
        print(paste0(drug_name," : ",anno_tmp))
        fig_path = paste0(filepath,"/04.Results/drug/depmap/gdsc/", Cancername, "/sig/")
        if(!dir.exists(fig_path)){
          dir.create(fig_path)
          print(paste0("Created folder: ", fig_path))
        } else {
          print(paste0("Folder already exists: ", fig_path))
        }
        setwd(fig_path)
        
        ggsave(filename = paste0(CancerType,"_",drug_name, "_nih_anti_treated_screening_sig.svg"), tmp_drug)
        pathway_links = link_genes_filtered_df %>% filter(Genes %in% unique(tmp_for_drug$PUTATIVE_TARGET)) %>% pull(Pathway)
        pathway_single = single_genes %>% filter(Genes %in% unique(tmp_for_drug$PUTATIVE_TARGET)) %>% pull(Pathway)
        
        tmp_sig_gdsc = data.frame(cancer = Cancername , 
                                  drug_name = drug_name , 
                                  target_genes = unique(tmp_for_drug$PUTATIVE_TARGET), 
                                  critical_features = ifelse( length(unique(c(pathway_links,pathway_single))) == 0 , 
                                                              NA , 
                                                              unique(c(pathway_links,pathway_single))))
        total_gdsc_screening = rbind(total_gdsc_screening , tmp_sig_gdsc)
        
      } else {
        fig_path = paste0(filepath,"/04.Results/drug/depmap/gdsc/", Cancername, "/no_sig/")
        if(!dir.exists(fig_path)){
          dir.create(fig_path)
          print(paste0("Created folder: ", fig_path))
        } else {
          print(paste0("Folder already exists: ", fig_path))
        }
        setwd(fig_path)

        ggsave(filename = paste0(CancerType,"_",drug_name, "_nih_anti_treated_screening_nosig.svg"), tmp_drug)
      }
      remove(anno_tmp)
    }
  }
  
  fig_path = paste0(filepath,"/04.Results/drug/depmap/gdsc/", Cancername, "/sig/")
  if(!dir.exists(fig_path)){
    dir.create(fig_path)
    print(paste0("Created folder: ", fig_path))
  } else {
    print(paste0("Folder already exists: ", fig_path))
  }
  setwd(fig_path)
  
  write.csv(total_gdsc_screening, paste0("~/nas/04.Results/drug/depmap/gdsc/", Cancername, "/sig/", Cancername,"_gdsc_screening_sig.csv"))
  # lihc = readRDS(paste0("~/nas/04.Results/short_long/", CancerType,"_critical_features_short_long_with_drug.rds"))
}
