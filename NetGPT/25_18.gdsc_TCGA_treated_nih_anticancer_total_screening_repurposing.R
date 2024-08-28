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

# num_CancerType = "11.TCGA-STAD"

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
    print(drug_name)
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
        ggtitle(paste0(drug_name ,"\n",
                       "TCGA : " , paste0(TCGA_from, collapse = ","),"\n",
                       "nih : " , paste0(nih_from, collapse = ","),"\n",
                       "antifund : " , paste0(anti_from, collapse = ","))
                ) +
        # stat_compare_means(label.y = 10) +
        theme_minimal()
 
      if (anno_tmp < 0.05 ) {
        print(paste0(drug_name," : ",anno_tmp))
        if (mean(tmp_anti_long$Z_SCORE) > mean(tmp_anti_short$Z_SCORE) & mean(tmp_anti_short$Z_SCORE) < 0) {
          
          fig_path = paste0(filepath,"/04.Results/drug/depmap/gdsc/", Cancername, "/repur_sig/")
          if(!dir.exists(fig_path)){
            dir.create(fig_path)
            print(paste0("Created folder: ", fig_path))
          } else {
            print(paste0("Folder already exists: ", fig_path))
          }
          setwd(fig_path)
          
          ggsave(filename = paste0(CancerType,"_",drug_name, "_repur_screening_sig.svg"), tmp_drug)
          
          tmp_repur_screen = data.frame(cancer_name = Cancername, 
                                        drug_name = drug_name, 
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
          
          total_repur_screening = rbind(total_repur_screening, tmp_repur_screen)
        }
        
      } else {
        fig_path = paste0(filepath,"/04.Results/drug/depmap/gdsc/", Cancername, "/no_sig/")
        if(!dir.exists(fig_path)){
          dir.create(fig_path)
          print(paste0("Created folder: ", fig_path))
        } else {
          print(paste0("Folder already exists: ", fig_path))
        }
        setwd(fig_path)
        
        ggsave(filename = paste0(CancerType,"_",drug_name, "_repur_screening_nosig.svg"), tmp_drug)
      }
      remove(anno_tmp)
    }
    
  }
  # lihc = readRDS(paste0("~/nas/04.Results/short_long/", CancerType,"_critical_features_short_long_with_drug.rds"))
}

write.xlsx(total_repur_screening , "~/nas/04.Results/drug/depmap/gdsc/total_repurposing_screening_w_target.xlsx")

