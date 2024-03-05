
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

# num_CancerType = "10.TCGA-BLCA"

for (num_CancerType in Cancerlist) {
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  Cancername = gsub('TCGA-' , '', CancerType)
  
  # call input
  gc_cellline = readRDS(paste0("~/nas/00.data/filtered_TCGA/", num_CancerType, "/",Cancername,"_cellline_dual_all_log.rds"))
  gc_TCGA = readRDS(paste0("~/nas/04.Results/short_long/", CancerType,"_critical_features_short_long.rds"))
  cli_drug = read.csv(paste0(ref_path,"TCGA_clinical_drug/",Cancername, "_clinical_druginfo.csv"))
  
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
  
  cli_drug = cli_drug %>% select(-X)
  cli_drug_filt = cli_drug %>% 
    slice(-1,-2) %>% 
    filter(pharmaceutical_therapy_drug_name != "[Not Available]") %>%
    select(-bcr_patient_uuid, -bcr_drug_uuid)

  common_pat = intersect(rownames(gc_TCGA),cli_drug_filt$bcr_patient_barcode)
  
  gc_TCGA_drug = gc_TCGA[common_pat,]
  cli_drug_filt = cli_drug_filt %>% filter(bcr_patient_barcode %in% common_pat)

  gdsc_w_cluster_filt = gdsc_w_cluster %>% filter(DRUG_NAME %in% unique(cli_drug_filt$pharmaceutical_therapy_drug_name))

  ###  
  for (cli_pat in unique(cli_drug_filt$bcr_patient_barcode)) {
    tmp_pat = cli_drug_filt %>% filter(bcr_patient_barcode == cli_pat)
    if (length(tmp_pat$pharmaceutical_therapy_drug_name) >1) {
      gc_TCGA_drug[cli_pat,"treated_drug"] = paste(tmp_pat$pharmaceutical_therapy_drug_name , collapse = ",")
    } else {
      gc_TCGA_drug[cli_pat,"treated_drug"] = tmp_pat$pharmaceutical_therapy_drug_name
    }
    
  }
  
  sig_treated_target = data.frame()
  
  # for fic
  
  fig_path = paste0(filepath,"/04.Results/drug/depmap/gdsc/", Cancername)
  if(!dir.exists(fig_path)){
    dir.create(fig_path)
    print(paste0("Created folder: ", fig_path))
  } else {
    print(paste0("Folder already exists: ", fig_path))
  }
  setwd(fig_path)
  
  # drug_name = "Temsirolimus"
  for (drug_name in unique(gdsc_w_cluster_filt$DRUG_NAME)) {
    # n = n+1
    print(drug_name)
    tmp_for_drug = gdsc_w_cluster_filt %>% filter(DRUG_NAME == drug_name)
    
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
        ggsave(filename = paste0(CancerType,"_",drug_name, "_treated_TCGA_screening_sig.svg"), tmp_drug)
        
        tcga_drug = readRDS( paste0("~/nas/04.Results/short_long/", CancerType,"_critical_features_short_long_with_drug.rds"))
        
        tcga_drug_surv = tcga_drug %>%  
          mutate(drug_cluster = case_when(
            str_detect(treated_drug, drug_name) ~ "treated",
            TRUE ~ "untreated"))
        
        if (sum(tcga_drug_surv$drug_cluster == "treated") == 0) {
          next
        } else {
          
          fit = survfit(Surv(duration, status) ~ drug_cluster, data = tcga_drug_surv)
          
          svglite(filename =  paste0(CancerType,"_",drug_name, "_treated_TCGA_survival.svg"), 
                  bg = "white", 
                  pointsize = 12)
          
          tmp_drug_surv = ggsurvplot(fit, data = tcga_drug_surv, risk.table = TRUE,
                                     palette = "jco", pval = TRUE, surv.median.line = "hv", xlab = "days")
          
          print(tmp_drug_surv)
          dev.off()
    
        }
       
      } else {
        ggsave(filename = paste0(CancerType,"_",drug_name, "_treated_TCGA_screening_nosig.svg"), tmp_drug)
      }
     remove(anno_tmp)
    }
  }
  # lihc = readRDS(paste0("~/nas/04.Results/short_long/", CancerType,"_critical_features_short_long_with_drug.rds"))
}
