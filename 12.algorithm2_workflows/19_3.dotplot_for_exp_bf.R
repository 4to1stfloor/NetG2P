library(stringr)
library(data.table)
library(ComplexHeatmap)
library(TCGAbiolinks)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(pheatmap)
library(readxl)
library(org.Hs.eg.db)
library(clusterProfiler)
library(dplyr)

filepath = "/home/seokwon/nas/"
ref_path = paste0(filepath, "99.reference/")

link_genes = readRDS(paste0(ref_path, "/KEGG_pathway_shared_genes.rds"))

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
Cancerlist = Cancerlist[-7]
   
# for all
for (num_CancerType in Cancerlist) {
  
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  
  # call input
  
  sce = readRDS(paste0(sce_path,num_CancerType,"/", CancerType,".sce_raw_snv_cnv.rds"))
  
  best_features = readRDS(paste0(main.path_tc,"/h2o_bias_pval_dual_cut_50/network/",CancerType,"_best_features_links.rds"))
  best_features$pathwaylinks = ifelse(best_features$from == best_features$to,best_features$from,  paste0(best_features$from, best_features$to) )
  annotate_best_features =  read.csv(paste0(filepath,"13.analysis/",CancerType,"_best_features.csv"))
  duration_log_df = readRDS(paste0(main.path_tc, "/", CancerType,"_dual_add_duration_log.rds"))
  
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
  
  best_features_df = duration_log_df[,annotate_best_features$variable]
  best_features_df$vitalstatus = duration_log_df$vitalstatus
  best_features_df$duration = duration_log_df$duration
  
  best_features_df$status = NA
  best_features_df$status[which(best_features_df$vitalstatus == "Dead")] = 1
  best_features_df$status[which(best_features_df$vitalstatus == "Alive")] = 0
  best_features_df$vitalstatus = NULL
  bf_zs_shared_exp_filt_df$vitalstatus = NULL
  pathwaylink_num = readxl::read_xlsx(paste0(main.path_tc,"/cluster_fig/",CancerType,"_result_survpval.xlsx"),sheet = "best")
  
  # 18 means best p-value when devide two cluster
  out = pheatmap::pheatmap((best_features_df[,1:as.numeric(pathwaylink_num$num)] > -log(0.05))*1 , cluster_cols = T,
                           cluster_rows = T,labels_cols = "", 
                           show_rownames = T)
  
  tmp_pheat_cut = as.data.frame (cutree(out$tree_row, 2) , out[["tree_row"]][["labels"]])
  colnames(tmp_pheat_cut) = "cluster"
  
  if (CancerType == "TCGA-KIDNEY") {
    bf_zs_shared_exp_filt_df = bf_zs_shared_exp_filt_df[rownames(best_features_df),]
    bf_zs_shared_exp_filt_df$cluster = tmp_pheat_cut$cluster
  } else if (all.equal(rownames(bf_zs_shared_exp_filt_df), rownames(best_features_df)) == TRUE) {
    bf_zs_shared_exp_filt_df$cluster = tmp_pheat_cut$cluster
  } else {
    bf_zs_shared_exp_filt_df = bf_zs_shared_exp_filt_df[rownames(best_features_df),]
    bf_zs_shared_exp_filt_df$cluster = tmp_pheat_cut$cluster
  }
  
  bad_group = bf_zs_shared_exp_filt_df[which(bf_zs_shared_exp_filt_df$cluster == 1),]
  bad_group = na.omit(bad_group)
  good_group = bf_zs_shared_exp_filt_df[which(bf_zs_shared_exp_filt_df$cluster == 2),]
  good_group = na.omit(good_group)
  
  bad_group[bad_group > 5] <- 5
  bad_group[bad_group < -5] <- -5
  
  good_group[good_group > 5] <- 5
  good_group[good_group < -5] <- -5

  bad_group$cluster = "bad" 
  good_group$cluster = "good" 
  
  total_group = rbind(good_group,bad_group)
  
  pvals <- as.data.frame(matrix(nrow = c(ncol(total_group)-1)))
  rownames(pvals) = colnames(total_group)[-length(total_group)]
  colnames(pvals) = "pval"
  for(i in 1:c(ncol(total_group)-1)) {
    gene_expr_bad <- bad_group[,i]
    gene_expr_good <- good_group[,i]
    ttest_result <- t.test(gene_expr_bad, gene_expr_good)
    pvals[i,"pval"] <- ttest_result$p.value
  }
  
  deg_good_bad = rownames(pvals)[which(pvals < 0.05)]
  deg_group = total_group[,c(deg_good_bad,"cluster")]
  
  bad_cluster1_gene = c()
  good_cluster2_gene = c()
  for (deg_genes in deg_good_bad) {
    if (mean(deg_group[which(deg_group$cluster == "bad"), deg_genes]) > mean(deg_group[which(deg_group$cluster == "good"), deg_genes])) {
      bad_cluster1_gene <- c(bad_cluster1_gene, deg_genes)
    } else {
      good_cluster2_gene <- c(good_cluster2_gene, deg_genes)
    }
  }
  
  bad_cluster1_gene_en = AnnotationDbi::select(org.Hs.eg.db, bad_cluster1_gene, 'ENTREZID', 'SYMBOL')[
    which(!is.na(AnnotationDbi::select(org.Hs.eg.db, bad_cluster1_gene, 'ENTREZID', 'SYMBOL')$ENTREZID)),]
  
  bad_cluster1_gene_en <- data.frame(bad_cluster1_gene_en, row.names = NULL)
  
  good_cluster2_gene_en = AnnotationDbi::select(org.Hs.eg.db, good_cluster2_gene, 'ENTREZID', 'SYMBOL')[
    which(!is.na(AnnotationDbi::select(org.Hs.eg.db, good_cluster2_gene, 'ENTREZID', 'SYMBOL')$ENTREZID)),]
  good_cluster2_gene_en <- data.frame(good_cluster2_gene_en, row.names = NULL)
  
  top_genes_group = list(bad_cluster = bad_cluster1_gene_en$ENTREZID,good_cluster = good_cluster2_gene_en$ENTREZID)
  ck <- compareCluster(geneCluster = top_genes_group, fun = "enrichKEGG")
  ck <- setReadable(ck, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
  
  fig_path = paste0(main.path_tc,"/",CancerType,"_analysis_exp_bf/")
  setwd(fig_path)
  
  png(filename = paste0(CancerType,"_dotplot.png"),
      width = 25, height = 25,  units = "cm" ,pointsize = 12,
      bg = "white", res = 1200, family = "")
  
  exp_dot = dotplot(ck)
  
  print(exp_dot)
  dev.off()
  

}  






