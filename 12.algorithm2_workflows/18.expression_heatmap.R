
## expression divide by common genes
library(stringr)
library(data.table)
library(ComplexHeatmap)
library(TCGAbiolinks)
library(SingleCellExperiment)
library(SummarizedExperiment)

filepath = "/home/seokwon/nas/"
ref_path = paste0(filepath, "99.reference/")

link_genes = readRDS(paste0(ref_path, "/KEGG_pathway_shared_genes.rds"))

folder_name = "h2o_bias_pval_dual_cut_50/"

# zscore normalization
tcga.calc.zscore = function(sce, target.genes){
  message("Calculating z-score with respect to all diploid cells. Version 2022.11.07")
  common.genes = intersect(rownames(sce), target.genes)
  if (length(common.genes) == 0) {
    stop("None of the target genes found in sce. Please check your nomenclature")
  } else if (all(common.genes %in% target.genes)) {
    print("Check!")
  }
  sce.sub = subset(sce, rownames(sce) %in% target.genes,)
  #i am not checking assay names.
  count.mat = assay(sce.sub, 1)
  cnv.mat = assay(sce.sub, 3)
  z.mat = matrix(data = NA, nrow = length(target.genes), ncol = ncol(count.mat))
  colnames(z.mat) = colnames(count.mat)
  rownames(z.mat) = target.genes
  for (i in 1:nrow(count.mat)) {
    idx.di = which(cnv.mat[i,] == 2)
    query.mean = mean(count.mat[i, idx.di], na.rm = T)
    query.sd = sd(count.mat[i, idx.di], na.rm = T)
    idx.symb = which(rownames(z.mat) == rownames(count.mat)[i])
    z.mat[i,] = (count.mat[i,] - query.mean)/query.sd
  }
  return(z.mat)
}

Cancerlist = dir(paste0(filepath, "/00.data/filtered_TCGA/"))
sce_path = "/mnt/gluster_server/data/raw/TCGA_data/00.data/"
Cancerlist = Cancerlist[-1:-11]
# for all
for (num_CancerType in Cancerlist) {
  
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  
  # call input
  duration_log_df = readRDS(paste0(main.path_tc,"/", CancerType, "_dual_add_duration_log.rds"))

  sce = readRDS(paste0(sce_path,num_CancerType,"/", CancerType,".sce_raw_snv_cnv.rds"))
  
  best_features = readRDS(paste0(main.path_tc,"/h2o_bias_pval_dual_cut_50/network/",CancerType,"_best_features_links.rds"))
  best_features$pathwaylinks = ifelse(best_features$from == best_features$to,best_features$from,  paste0(best_features$from, best_features$to) )
  
  # only has links
  total_link_genes <- unlist(strsplit(link_genes[which(link_genes$pathway_name %in% best_features$pathwaylinks),]$shared_genes, ","))

  # zscore normalizaiton 
  best_features_shared_zscore_exp_df = tcga.calc.zscore(sce = sce, target.genes = unique(total_link_genes))
  # (1) - filtering only 01A
  bf_sh_zs_exp_wo = as.data.frame(na.omit(best_features_shared_zscore_exp_df))
  
  if (CancerType %in% c("TCGA-COADREAD","TCGA-KIDNEY")) {
    cols_to_keep <- grep("^.*-01A-", names(bf_sh_zs_exp_wo), value = TRUE)
  } else {
    cols_to_keep <- grep("^.*-01A$", names(bf_sh_zs_exp_wo), value = TRUE)
  }

  bf_sh_zs_exp_wo_filt <- bf_sh_zs_exp_wo[, cols_to_keep]

  dup_names = grep("\\.1$",names(bf_sh_zs_exp_wo_filt),value = TRUE)
  if (CancerType %in% c("TCGA-COADREAD","TCGA-KIDNEY")) {
    dup_names_ori = substr(dup_names, 1, 28)
  }else {
    dup_names_ori = substr(dup_names, 1, 16)
  }
  
  
  if (length(dup_names) !=0) {
    for (dup_patients in 1:length(dup_names)) {
      first = dup_names[dup_patients]
      second = dup_names_ori[dup_patients]
      tmp_dup = as.data.frame(bf_sh_zs_exp_wo_filt[,which(names(bf_sh_zs_exp_wo_filt) == first)])
      tmp_dup2 = data.frame(lapply(tmp_dup, as.numeric))
      rownames(tmp_dup2) = rownames(bf_sh_zs_exp_wo_filt)
      colnames(tmp_dup2) = "first"
      
      tmp_dup3 = as.data.frame(bf_sh_zs_exp_wo_filt[,which(names(bf_sh_zs_exp_wo_filt) == second)])
      tmp_dup4 = data.frame(lapply(tmp_dup3, as.numeric))
      rownames(tmp_dup4) = rownames(bf_sh_zs_exp_wo_filt)
      colnames(tmp_dup4) = "second"
      averages <- as.data.frame(rowMeans(cbind(tmp_dup2, tmp_dup4)))
      colnames(averages) = second
      bf_sh_zs_exp_wo_filt[,first] =NULL
      bf_sh_zs_exp_wo_filt[,second] = NULL
      bf_sh_zs_exp_wo_filt = cbind(bf_sh_zs_exp_wo_filt,averages)
      
    }
  } 
  

  if (sum(duplicated(names(bf_sh_zs_exp_wo_filt))) == 0 ) {
    colnames(bf_sh_zs_exp_wo_filt) = substr(colnames(bf_sh_zs_exp_wo_filt), 1, 12)
  }
 
  # bf_sh_zs_exp_wo_filt <- data.frame(lapply(bf_sh_zs_exp_wo_filt_df, as.numeric))
  # bf_sh_zs_exp_wo_filt = bf_sh_zs_exp_wo_filt_df  
  # rownames(bf_sh_zs_exp_wo_filt) = rownames(bf_sh_zs_exp_wo_filt_df)
  
  
  # patients filt
  
  if (CancerType %in% c("TCGA-COADREAD","TCGA-KIDNEY")) {
    bf_zs_shared_exp_filt_df= bf_sh_zs_exp_wo_filt[,intersect(colnames(bf_sh_zs_exp_wo_filt) , rownames(duration_log_df))]
    
  } else {
    bf_zs_shared_exp_filt_df= bf_sh_zs_exp_wo_filt[,rownames(duration_log_df)]
  }
  
  
  # add vital_status
  bf_zs_shared_exp_filt_df = as.data.frame(t(bf_zs_shared_exp_filt_df))
  bf_zs_shared_exp_filt_df$vitalstatus  = NA
  
  for (patients in rownames(bf_zs_shared_exp_filt_df)) {
    bf_zs_shared_exp_filt_df[patients,]$vitalstatus = duration_log_df[patients,]$vitalstatus
  }
  
  # for complex_heatmap divide by a specific column annotation 
  col_colors = list(vitalstatus = c("black", "yellow"))
  col_ann = subset(bf_zs_shared_exp_filt_df[order(bf_zs_shared_exp_filt_df$vitalstatus),],select = vitalstatus)
  names(col_colors$vitalstatus) = unique(col_ann$vitalstatus)
  bf_zs_shared_exp_filt_df = bf_zs_shared_exp_filt_df[rownames(col_ann),]
  
  cheatmap_mat = as.matrix(t(bf_zs_shared_exp_filt_df[,-length(bf_zs_shared_exp_filt_df)]))
  # cheatmap_mat= as.matrix(lapply(cheatmap_mat, as.numeric))
  cheatmap_mat[cheatmap_mat > 4] = 4
  cheatmap_mat[cheatmap_mat < -4] = -4
  
  dir_path = paste0(filepath, "/00.data/filtered_TCGA/",num_CancerType,"/heatmap_with_best_features/")
  dir.create(dir_path)
  
  png(filename = paste0(dir_path,CancerType,"_with_cluster.png"),
      width = 1000, height = 1000, units = "px", pointsize = 12,
      bg = "white", res = NA, family = "")
  
  out = ComplexHeatmap::pheatmap(cheatmap_mat, 
                                 column_split = col_ann$vitalstatus,
                                 labels_col = "",
                                 show_rownames = T, 
                                 show_colnames = F, 
                                 annotation_col = col_ann,
                                 annotation_colors = col_colors,
                                 # clustering_method = "average",
                                 cluster_cols = T,
                                 cluster_rows = T)
  print(out)
  
  dev.off()
  
  png(filename = paste0(dir_path,CancerType,"_without_cluster.png"),
      width = 1000, height = 1000, units = "px", pointsize = 12,
      bg = "white", res = NA, family = "")
  
  out2 = pheatmap(cheatmap_mat, 
                  labels_col = "",
                  show_rownames = T, 
                  show_colnames = F, 
                  annotation_colors = col_colors,
                  annotation_col = col_ann,
                  # clustering_method = "average",
                  cluster_cols = T,
                  cluster_rows = T)
  print(out2)
  
  dev.off()
  
  remove(total_link_genes,
         best_features,
         best_features_shared_zscore_exp_df,
         bf_sh_zs_exp_wo_filt,
         bf_sh_zs_exp_wo,
         bf_zs_shared_exp_filt_df,col_colors,col_ann, cheatmap_mat, out,out2)
}



