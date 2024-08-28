library(survminer) 
library(ggplot2) 
library(pheatmap)
library(survival) 
library(data.table) 
library(RColorBrewer) 
library(NbClust)
library(cluster) 
library(factoextra) 
library(dplyr)
library(stringr)
require(RSNNS)
require(clusterGeneration)

filepath = "/home/seokwon/nas/"
ref_path = paste0(filepath, "99.reference/")
type = "each"
bi_num_mode = "log"

Cancerlist = dir(paste0(filepath, "/00.data/filtered_TCGA/"))
# num_CancerType = "04.TCGA-CESC"
for (num_CancerType in Cancerlist) {
  
  # 1) prepare for NN
  # (1) for TCGA data upload
  
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  Cancername = gsub('TCGA-','', CancerType)
  if (CancerType == "Total-TCGA") {
    break
  }
  
  for (type in c("each","link")) {
    
    #   A. for each pathway
    if (type == "each") {
      phyper = readRDS(paste0(main.path_tc, "/phyper_cellline_enrichment_54_KEGG_GC_20220117.rds"))
    } else if (type == "link") {
      #   B. for pathwaylink
      phyper = readRDS(paste0(main.path_tc, "/phyper_cellline_enrichment_54_KEGG_net_GC_20220117.rds"))
    }
    
    phyper.filt = phyper
    
    #  transfer data 1 (>0.05) or 100 (<0.05)
    
    phyper.filt_t = t(phyper.filt)
    if (length(which(colSums(phyper.filt_t) == dim(phyper.filt_t)[1])) !=0) {
      phyper.filt_t_num= phyper.filt_t[,-which(colSums(phyper.filt_t) == dim(phyper.filt_t)[1])]
    }else {
      phyper.filt_t_num = phyper.filt_t
    }
    
    if (bi_num_mode == "log") {
      phyper.filt_t_num.log = -log(phyper.filt_t_num)
    } else if (bi_num_mode == "bi") {
      phyper.filt_t_num.log = ifelse(phyper.filt_t_num < 0.05, 1 , 0)
    } 
    
    # rownames(phyper.filt_t_num.log) = substr(rownames(phyper.filt_t_num.log), 1,12)
    
    df.phyper.filt_t_log_num = as.data.frame(phyper.filt_t_num.log)
    df.phyper.filt_t_log_num
    
    if (type == "link") {
      saveRDS(df.phyper.filt_t_log_num, file = paste0(main.path_tc,"/",Cancername,"_cellline_pathwaylink_all_log.rds"))
    } else if (type == "each") {
      saveRDS(df.phyper.filt_t_log_num, file = paste0(main.path_tc,"/",Cancername,"_cellline_pathwayeach_all_log.rds"))
    }
  }
  
  cellline_link = readRDS(paste0(main.path_tc,"/",Cancername,"_cellline_pathwaylink_all_log.rds"))
  cellline_each = readRDS(paste0(main.path_tc,"/",Cancername,"_cellline_pathwayeach_all_log.rds"))
  
  if (all.equal(rownames(cellline_link) , rownames(cellline_each))) {
    cellline_total = cbind(cellline_link,cellline_each)
  } else {
    print("error")
  }
  
  saveRDS(cellline_total, file = paste0(main.path_tc,"/",Cancername,"_cellline_dual_all_log.rds"))
}
