
library(survival) 
library(survminer) 
library(readxl)
library(tidyverse)
library(svglite)

filepath = "/home/seokwon/nas/"
ref_path = paste0(filepath, "99.reference/")

Cancerlist = dir(paste0(filepath, "/00.data/filtered_TCGA/"))
Cancerlist = Cancerlist[c(-11,-12)]

anticancer_drug = read_xlsx("~/nas/99.reference/anticancer_fund_cancerdrugsdb.xlsx")

anticancer_drug_filted = anticancer_drug %>%
  separate_rows(Indications, sep = "\\;") %>%
  mutate(Indications = str_trim(Indications))

meta_for_TCGA = read.csv(paste0(ref_path, "/Depmap_meta_filt_to_TCGA.csv"))
gdsc = readRDS("/mnt/gluster_server/data/reference/GDSC/2024_01_11/GDSC_data_combined.rds")
meta_cell = readRDS("/mnt/gluster_server/data/reference/TLDR/meta_cells_primary.rds")
criteria = read.csv("~/nas/99.reference/drug_alias_single_SW.csv", encoding = "UTF-8")
criteria_filt = criteria %>%
  separate_rows(alias, sep = ",") %>%
  mutate(alias = trimws(alias))

# num_CancerType = "04.TCGA-CESC"
# setwd("~/nas/04.Results/drug/TCGA_drug/TCGA_drug/")

for (num_CancerType in Cancerlist) {
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  Cancername = gsub('TCGA-' , '', CancerType)
  
  # call input
  gc_cellline = readRDS(paste0("~/nas/00.data/filtered_TCGA/", num_CancerType, "/",Cancername,"_cellline_dual_all_log.rds"))
  gc_TCGA = readRDS(paste0("~/nas/04.Results/short_long/", CancerType,"_critical_features_short_long.rds"))
  cli_drug = read.csv(paste0(ref_path,"TCGA_clinical_drug/",Cancername, "_drug_info_update.csv"))
  unique(cli_drug$pharmaceutical_therapy_drug_name)
  # tcga_drug = readRDS( paste0("~/nas/04.Results/short_long/", CancerType,"_critical_features_short_long_with_drug.rds"))
  
  # cli_drug = cli_drug %>% select(-X)
  cli_drug_filt = cli_drug %>%
    slice(-1,-2) %>%
    filter(pharmaceutical_therapy_drug_name != "[Not Available]") %>%
    select(-bcr_patient_uuid, -bcr_drug_uuid)

  # cli_drug_filt_edit = cli_drug %>% 
  #   slice(-1,-2) %>% 
  #   filter(pharmaceutical_therapy_drug_name != "[Not Available]") %>%
  #   select(-bcr_patient_uuid, -bcr_drug_uuid)
  
  cli_drug_filt_edit = left_join(cli_drug_filt, criteria_filt, by = c("pharmaceutical_therapy_drug_name" = "alias")) %>%
    mutate(main_name_new = coalesce(main_name, pharmaceutical_therapy_drug_name)) %>%
    select(bcr_patient_barcode,main_name_new, main_name,pharmaceutical_therapy_drug_name, everything())

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
  
  gdsc_w_cluster = left_join(gdsc_w_cluster, criteria_filt, by = c("DRUG_NAME" = "alias")) %>%
    mutate(DRUG_NAME_new = coalesce(main_name, DRUG_NAME)) %>%
    select(DepMap.ID,cluster , DRUG_ID,DRUG_NAME_new , DRUG_NAME,main_name, everything())
  
  common_pat = intersect(rownames(gc_TCGA),cli_drug_filt_edit$bcr_patient_barcode)
  
  gc_TCGA_drug = gc_TCGA[common_pat,]
  
  for (cli_pat in unique(cli_drug_filt_edit$bcr_patient_barcode)) {
    tmp_pat = cli_drug_filt_edit %>% filter(bcr_patient_barcode == cli_pat)
    if (length(tmp_pat$main_name_new ) >1) {
      gc_TCGA_drug[cli_pat,"treated_drug"] = paste(tmp_pat$main_name_new  , collapse = ",")
    } else {
      gc_TCGA_drug[cli_pat,"treated_drug"] = tmp_pat$main_name_new 
    }
    
  }
  tcga_drug = gc_TCGA_drug
  ## TCGA_drug

  fig_path = paste0(filepath,"/04.Results/drug/TCGA_drug/", Cancername)
  if(!dir.exists(fig_path)){
    dir.create(fig_path)
    print(paste0("Created folder: ", fig_path))
  } else {
    print(paste0("Folder already exists: ", fig_path))
  }
  setwd(fig_path)
  
  tmp_tcga_drug = data.frame()
  tmp_gdsc_drug = data.frame()
  
  for (uni_drug in unique(cli_drug_filt_edit$main_name_new )) {
    num_pat = cli_drug_filt_edit %>% filter(main_name_new  == uni_drug) %>% nrow()
    num_cell = gdsc_w_cluster %>% filter(DRUG_NAME_new == uni_drug) %>% pull(DepMap.ID) %>% unique() %>% length()
    
    tmp_dd = data.frame(treated_drug = uni_drug , number_of_patient = num_pat)
    tmp_dc = data.frame(treated_drug = uni_drug , number_of_cellline = num_cell)
    
    tmp_tcga_drug = rbind(tmp_tcga_drug , tmp_dd)
    tmp_gdsc_drug = rbind(tmp_gdsc_drug , tmp_dc)
    
    tcga_drug_surv = tcga_drug %>%  
      mutate(drug_cluster = case_when(
        str_detect(treated_drug, uni_drug) ~ "treated",
        TRUE ~ "untreated")) 
    
    tcga_drug_surv_filt = tcga_drug_surv %>% 
      filter(drug_cluster == "treated")
    
    if (length(table(tcga_drug_surv_filt$cluster)) != 2) {
      next
    } else {
      
      fit = survfit(Surv(duration, status) ~ cluster, data = tcga_drug_surv_filt)
      
      if (grepl("/",uni_drug )) {
        uni_drug = paste(unlist(str_split(uni_drug , "/")), collapse = "_")
      } 
      
      svglite(filename =  paste0(CancerType,"_",uni_drug, "_sl_treated_TCGA_survival.svg"), 
              bg = "white", 
              pointsize = 12)
      
      tmp_drug_surv = ggsurvplot(fit, data = tcga_drug_surv_filt, risk.table = TRUE,
                                 palette = "npg", pval = TRUE, surv.median.line = "hv", xlab = "days")
      
      print(tmp_drug_surv)
      dev.off()
      }
  }
  
  tmp_tcga_drug_sort = tmp_tcga_drug %>% arrange(desc(number_of_patient))
  tmp_gdsc_drug_sort = tmp_gdsc_drug %>% arrange(desc(number_of_cellline))
  
  write.csv(tmp_tcga_drug_sort, paste0(Cancername, "_clinical_drug_table.csv"))
  write.csv(tmp_gdsc_drug_sort, paste0(Cancername, "_cellline_drug_table.csv"))
  
}
