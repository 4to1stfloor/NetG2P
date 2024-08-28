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
library(DESeq2)
library(survival) 
library(survminer) 
library(DOSE)
library(UpSetR)
library(enrichplot)
library(dplyr)
library(readxl)
library(ggupset)

filepath = "/home/seokwon/nas/"
ref_path = paste0(filepath, "99.reference/")

link_genes = readRDS(paste0(ref_path, "/KEGG_pathway_shared_genes.rds"))
single_genes = readRDS(paste0(ref_path, "/Kegg_pathway_genes.rds"))

tcga.calc.zscore = function(sce, target.genes){
  message("Calculating z-score with respect to all diploid cells. Version 2023.05.03")
  common.genes = intersect(rownames(sce), target.genes)
  if (length(common.genes) == 0) {
    stop("None of the target genes found in sce. Please check your nomenclature")
  } else if (length(common.genes) != length(target.genes)) {
    message("Some of the genes from query does not exist in this cancer type. It will result in NAs")
    message("Missing genes are: ", paste(setdiff(target.genes, common.genes), collapse = ", "))
  }
  sce.sub = subset(sce, rownames(sce) %in% common.genes,)
  #i am not checking assay names.
  count.mat = assay(sce.sub, 1)
  cnv.mat = assay(sce.sub, 3)
  z.mat = matrix(data = NA, nrow = nrow(count.mat), ncol = ncol(count.mat))
  colnames(z.mat) = colnames(count.mat)
  rownames(z.mat) = rownames(count.mat)
  for (i in 1:nrow(count.mat)) {
    idx.di = which(cnv.mat[i,] == 2)
    query.mean = mean(count.mat[i, idx.di], na.rm = T)
    query.sd = sd(count.mat[i, idx.di], na.rm = T)
    z.mat[i,] = (count.mat[i,] - query.mean)/query.sd
  }
  return(z.mat)
}

Cancerlist = dir(paste0(filepath, "/00.data/filtered_TCGA/"))
sce_path = "/mnt/gluster_server/data/raw/TCGA_data/00.data/"

surv_total_results = read_xlsx("~/nas/04.Results/Total_results_survpval2.xlsx")

# for all
for (num_CancerType in Cancerlist) {
  
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  
  # call input
 
  sce = readRDS(paste0(sce_path,num_CancerType,"/", CancerType,".sce_raw_snv_cnv.rds"))
  short_long_features = read_xlsx(paste0(filepath, "04.Results/short_long/",CancerType, "_best_features_short_long.xlsx"))
  bf_short_long = readRDS(paste0(filepath, "04.Results/short_long/",CancerType,"_best_features_short_long.rds"))
  duration_log_df = readRDS(paste0(main.path_tc, "/", CancerType,"_dual_add_duration_log.rds"))
  
  # 일단 unique gene으로 해봄 
  total_link_genes <- unlist(strsplit(link_genes[which(link_genes$pathway_name %in% short_long_features$variable),]$shared_genes, ","))
  total_single_genes = single_genes[which(single_genes$Pathway %in% short_long_features$variable),]$Genes
  total_bf_genes = c(total_link_genes,total_single_genes)
  
  # zscore normalizaiton 
  best_features_shared_zscore_exp_df = tcga.calc.zscore(sce = sce, target.genes = unique(total_bf_genes))
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
  
  # patients filt
  if (CancerType %in% c("TCGA-COADREAD","TCGA-KIDNEY")) {
    bf_zs_shared_exp_filt_df= bf_sh_zs_exp_wo_filt[,intersect(colnames(bf_sh_zs_exp_wo_filt) , rownames(duration_log_df))]
    
  } else {
    bf_zs_shared_exp_filt_df= bf_sh_zs_exp_wo_filt[,rownames(duration_log_df)]
  }
  
  # add vital_status
  bf_zs_shared_exp_filt_df = as.data.frame(t(bf_zs_shared_exp_filt_df))
  
  bf_zs_shared_exp_filt_df$cluster  = NA
 
  bf_zs_shared_exp_filt_df[intersect(rownames(bf_zs_shared_exp_filt_df) , rownames(bf_short_long)),]$cluster = 
    bf_short_long[intersect(rownames(bf_zs_shared_exp_filt_df) , rownames(bf_short_long)),]$cluster

  bf_zs_shared_exp_filt_df = bf_zs_shared_exp_filt_df[which(!is.na(bf_zs_shared_exp_filt_df$cluster)),]
  
  short_gene = c()
  long_gene = c()
  
  # pvalue 일단 ㄴ (features 수준에서 pvalue컷 한거라)
  for (colnum in 1:(ncol(bf_zs_shared_exp_filt_df)-1)) {
    if (colMeans(bf_zs_shared_exp_filt_df[which(bf_zs_shared_exp_filt_df$cluster == "short"),][colnum]) >
        colMeans(bf_zs_shared_exp_filt_df[which(bf_zs_shared_exp_filt_df$cluster == "long"),][colnum])) {
      short_gene = c(short_gene , colnames(bf_zs_shared_exp_filt_df)[colnum])
    } else {
      long_gene =  c(long_gene,colnames(bf_zs_shared_exp_filt_df)[colnum])
    }
  }
  
  short_gene_en = AnnotationDbi::select(org.Hs.eg.db, short_gene, 'ENTREZID', 'SYMBOL')[
    which(!is.na(AnnotationDbi::select(org.Hs.eg.db, short_gene, 'ENTREZID', 'SYMBOL')$ENTREZID)),]
  short_gene_en <- data.frame(short_gene_en, row.names = NULL)
  
  long_cluster2_gene_en = AnnotationDbi::select(org.Hs.eg.db, long_gene, 'ENTREZID', 'SYMBOL')[
    which(!is.na(AnnotationDbi::select(org.Hs.eg.db, long_gene, 'ENTREZID', 'SYMBOL')$ENTREZID)),]
  long_cluster2_gene_en <- data.frame(long_cluster2_gene_en, row.names = NULL)

  top_genes_group = list(short_cluster = short_gene_en$ENTREZID,long_cluster = long_cluster2_gene_en$ENTREZID)
  
  ck <- compareCluster(geneCluster = top_genes_group, fun = "enrichKEGG")
  ck <- setReadable(ck, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
  
  # for top300 (from heatmap)
  fig_path = paste0(filepath,"04.Results/short_long/analysis_pathway_from_features/",CancerType,"_analysis_exp_bf/")
  if(!dir.exists(fig_path)){
    dir.create(fig_path)
    print(paste0("Created folder: ", fig_path))
  } else {
    print(paste0("Folder already exists: ", fig_path))
  }
  setwd(fig_path)
  
  png(filename = paste0(CancerType,"_cnetplot.png"),
      width = 25, height = 25,  units = "cm" ,pointsize = 12,
      bg = "white", res = 1200, family = "")
  
  exp_cnet = cnetplot(ck)
  
  print(exp_cnet)
  dev.off()
  
  png(filename = paste0(CancerType,"_dotplot.png"),
      width = 25, height = 25,  units = "cm" ,pointsize = 12,
      bg = "white", res = 1200, family = "")
  
  exp_dot = dotplot(ck) 
  
  print(exp_dot)
  dev.off()
  
  # for disease enrichment
  
  edo_short <- enrichDGN(top_genes_group$short_cluster)
  edo_long <- enrichDGN(top_genes_group$long_cluster)
  
  edox_short <- setReadable(edo_short, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
  edox_long <- setReadable(edo_long, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
  
  # barplot
  png(filename = paste0(CancerType,"_DGN_long_barplot.png"),
      width = 25, height = 25,  units = "cm" ,pointsize = 12,
      bg = "white", res = 1200, family = "")
  
  disease_good = barplot(edo_long, showCategory=25)  + ggtitle("top 25")
  
  print(disease_good)
  dev.off()
  
  png(filename = paste0(CancerType,"_DGN_short_barplot.png"),
      width = 25, height = 25,  units = "cm" ,pointsize = 12,
      bg = "white", res = 1200, family = "")
  
  disease_bad = barplot(edo_short, showCategory=25)  + ggtitle("top 25")
  
  print(disease_bad)
  dev.off()
  
  # dotplot
  png(filename = paste0(CancerType,"_DGN_long_dotplot.png"),
      width = 25, height = 25,  units = "cm" ,pointsize = 12,
      bg = "white", res = 1200, family = "")
  
  disease_dot_good = dotplot(edo_long, showCategory=30) + ggtitle("dotplot for good")
  
  print(disease_dot_good)
  dev.off()
  
  png(filename = paste0(CancerType,"_DGN_short_dotplot.png"),
      width = 25, height = 25,  units = "cm" ,pointsize = 12,
      bg = "white", res = 1200, family = "")
  
  disease_dot_bad =  dotplot(edo_short, showCategory=30) + ggtitle("dotplot for bad")
  
  print(disease_dot_bad)
  dev.off()
  
  # tree plot
  edox_long <- pairwise_termsim(edo_long)
  png(filename = paste0(CancerType,"_DGN_long_treeplot.png"),
      width = 25, height = 25,  units = "cm" ,pointsize = 12,
      bg = "white", res = 1200, family = "")
  
  disease_tree_good =  treeplot(edox_long) + ggtitle("treeplot for good")
  
  print(disease_tree_good)
  dev.off()
  
  edox_short <- pairwise_termsim(edo_short)
  
  png(filename = paste0(CancerType,"_DGN_short_treeplot.png"),
      width = 25, height = 25,  units = "cm" ,pointsize = 12,
      bg = "white", res = 1200, family = "")
  
  disease_tree_bad =  treeplot(edox_short) + ggtitle("treeplot for bad")
  
  print(disease_tree_bad)
  dev.off()
  
  # emapplot
  
  png(filename = paste0(CancerType,"_DGN_long_emapplot.png"),
      width = 25, height = 25,  units = "cm" ,pointsize = 12,
      bg = "white", res = 1200, family = "")
  
  disease_emap_good = emapplot(edox_long) + ggtitle("emapplot for long")
  
  print(disease_emap_good)
  dev.off()
  
  png(filename = paste0(CancerType,"_DGN_short_emapplot.png"),
      width = 25, height = 25,  units = "cm" ,pointsize = 12,
      bg = "white", res = 1200, family = "")
  
  disease_emap_bad = emapplot(edox_short) + ggtitle("emapplot for short")
  
  print(disease_emap_bad)
  dev.off()
  
  # upsetplot
  png(filename = paste0(CancerType,"_DGN_long_upsetplot.png"),
      width = 25, height = 25,  units = "cm" ,pointsize = 12,
      bg = "white", res = 1200, family = "")
  
  disease_upset_good = upsetplot(edox_long) + ggtitle("upsetplot for long")
  
  print(disease_upset_good)
  dev.off()
  
  png(filename = paste0(CancerType,"_DGN_short_upsetplot.png"),
      width = 25, height = 25,  units = "cm" ,pointsize = 12,
      bg = "white", res = 1200, family = "")
  
  disease_upset_bad = upsetplot(edox_short) + ggtitle("usetplot for short")
  
  print(disease_upset_bad)
  dev.off()
  
  remove(ck,top_genes_group,bf_zs_shared_exp_filt_df,bf_sh_zs_exp_wo_filt)
}  


