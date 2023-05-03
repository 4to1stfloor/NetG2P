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

filepath = "/home/seokwon/nas/"
ref_path = paste0(filepath, "99.reference/")

link_genes = readRDS(paste0(ref_path, "/KEGG_pathway_shared_genes.rds"))
single_genes = readRDS(paste0(ref_path, "/Kegg_pathway_genes.rds"))

folder_name = "h2o_bias_pval_dual_cut_50/"

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
Cancerlist = Cancerlist[-7]
surv_total_results = read.csv("~/nas/04.Results/Total_results_survpval.csv")

# for all
for (num_CancerType in Cancerlist) {
  
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  
  # call input
  
  sce = readRDS(paste0(sce_path,num_CancerType,"/", CancerType,".sce_raw_snv_cnv.rds"))
  best_features = readRDS(paste0(main.path_tc,"/h2o_bias_pval_dual_cut_50/network/",CancerType,"_best_features_links.rds"))
  best_features = best_features[1:surv_total_results[which(surv_total_results$CancerType == CancerType),]$num_of_features,]
  
  best_features$pathwaylinks = ifelse(best_features$from == best_features$to,best_features$from,  paste0(best_features$from, best_features$to) )
  
  annotate_best_features =  read.csv(paste0(filepath,"13.analysis/",CancerType,"_best_features.csv"))
  duration_log_df = readRDS(paste0(main.path_tc, "/", CancerType,"_dual_add_duration_log.rds"))
  
  # 편의상 total_link_genes 라고 함 -> 사실 link + single임 single은 pathway안에 있는 유전자 모두 link 는 shared_gene만 
  total_link_genes <- unlist(strsplit(link_genes[which(link_genes$pathway_name %in% best_features$pathwaylinks),]$shared_genes, ","))
  total_single_genes = single_genes[which(single_genes$Pathway %in% best_features$pathwaylinks),]$Genes
  total_link_genes = c(total_link_genes,total_single_genes)

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
  
  # get best pval on survplot (It could be changed and also the method of )
  pathwaylink_num = surv_total_results[which(surv_total_results$CancerType == CancerType),]$num_of_features
  
  # 18 means best p-value when devide two cluster
  out = pheatmap::pheatmap((best_features_df[,1:as.numeric(pathwaylink_num)] > -log(0.05))*1 ,
                           cluster_cols = T,
                           cluster_rows = T,
                           labels_cols = "", 
                           show_rownames = T,
                           silent = T)
  
  tmp_pheat_cut = as.data.frame (cutree(out$tree_row, 2) , out[["tree_row"]][["labels"]])
  colnames(tmp_pheat_cut) = "cluster"
  
  # rownames(bf_zs_shared_exp_filt_df) %in% rownames(best_features_df)
  if (CancerType == "TCGA-KIDNEY") {
    bf_zs_shared_exp_filt_df = bf_zs_shared_exp_filt_df[rownames(best_features_df),]
    bf_zs_shared_exp_filt_df$cluster = tmp_pheat_cut$cluster
    
    } else if (all.equal(rownames(bf_zs_shared_exp_filt_df), rownames(best_features_df)) == TRUE) {
      bf_zs_shared_exp_filt_df$cluster = tmp_pheat_cut$cluster
      }  else {
        bf_zs_shared_exp_filt_df = bf_zs_shared_exp_filt_df[rownames(best_features_df),]
        bf_zs_shared_exp_filt_df$cluster = tmp_pheat_cut$cluster
        }

  bf_zs_shared_exp_filt_df$duration = duration_log_df$duration
  bf_zs_shared_exp_filt_df$vital_status = duration_log_df$vitalstatus
  bf_zs_shared_exp_filt_df$status = NA
  bf_zs_shared_exp_filt_df$status = ifelse(bf_zs_shared_exp_filt_df$vital_status == "Alive", 0 , 1)
  bf_zs_shared_exp_filt_df$vital_status = NULL
  bf_zs_shared_exp_filt_df = na.omit(bf_zs_shared_exp_filt_df)
  
  bf_zs_shared_exp_filt_df = bf_zs_shared_exp_filt_df[which((bf_zs_shared_exp_filt_df$duration >= 0)),]
  
  bf_zs_shared_exp_filt_df2 = bf_zs_shared_exp_filt_df[,which(!colnames(bf_zs_shared_exp_filt_df) %in% c("cluster","duration","status"))]
  bf_zs_shared_exp_filt_df2[bf_zs_shared_exp_filt_df2 >5] = 5
  bf_zs_shared_exp_filt_df2[bf_zs_shared_exp_filt_df2 <-5] = -5
  
  bf_zs_shared_exp_filt_max_df = cbind( bf_zs_shared_exp_filt_df2, bf_zs_shared_exp_filt_df[,which(colnames(bf_zs_shared_exp_filt_df) %in% c("cluster","duration","status"))])
  
  
  fit = survfit(Surv(duration, status) ~ cluster, data = bf_zs_shared_exp_filt_max_df)
  # ggsurvplot(fit, data = bf_zs_shared_exp_filt_df, risk.table = TRUE,
  #            palette = "jco", pval = TRUE, surv.median.line = "hv", xlab = "days")
  
  # hypothesis : cluster 1 = better prognosis
  if (mean(fit$surv[1:fit[['strata']][['cluster=1']]]) > mean(fit$surv[fit[['strata']][['cluster=1']] + 1: fit[['strata']][['cluster=2']]])) {
    bad_group = bf_zs_shared_exp_filt_max_df[which(bf_zs_shared_exp_filt_max_df$cluster == 2),]
    good_group = bf_zs_shared_exp_filt_max_df[which(bf_zs_shared_exp_filt_max_df$cluster == 1),]
    
  } else if (mean(fit$surv[1:fit[['strata']][['cluster=1']]]) < mean(fit$surv[fit[['strata']][['cluster=1']] + 1: fit[['strata']][['cluster=2']]])) {
    bad_group = bf_zs_shared_exp_filt_max_df[which(bf_zs_shared_exp_filt_max_df$cluster == 1),]
    good_group = bf_zs_shared_exp_filt_max_df[which(bf_zs_shared_exp_filt_max_df$cluster == 2),]
  } else {
    print("I don't know")
  }

  bad_group = na.omit(bad_group)
  good_group = na.omit(good_group)
  
  bad_out = pheatmap::pheatmap(bad_group %>% select(-cluster,-duration,-status), 
                               cluster_cols = T,
                               cluster_rows = T, 
                               labels_cols = "",
                               show_rownames = T,
                               silent = T)
  
  good_out = pheatmap::pheatmap(good_group %>% select(-cluster,-duration,-status), 
                                cluster_cols = T,
                                cluster_rows = T, 
                                labels_cols = "",
                                show_rownames = T,
                                silent = T)
 
  # cut genes 300 / total genes about 1100
  bad_cluster1_gene = head(colnames(bad_group[,bad_out$tree_col[["order"]]]), n= 300)
  good_cluster2_gene = head(colnames(good_group[,good_out$tree_col[["order"]]]), n= 300)
  
  bad_cluster1_gene_en = AnnotationDbi::select(org.Hs.eg.db, bad_cluster1_gene, 'ENTREZID', 'SYMBOL')[
    which(!is.na(AnnotationDbi::select(org.Hs.eg.db, bad_cluster1_gene, 'ENTREZID', 'SYMBOL')$ENTREZID)),]
  bad_cluster1_gene_en <- data.frame(bad_cluster1_gene_en, row.names = NULL)
  
  good_cluster2_gene_en = AnnotationDbi::select(org.Hs.eg.db, good_cluster2_gene, 'ENTREZID', 'SYMBOL')[
    which(!is.na(AnnotationDbi::select(org.Hs.eg.db, good_cluster2_gene, 'ENTREZID', 'SYMBOL')$ENTREZID)),]
  good_cluster2_gene_en <- data.frame(good_cluster2_gene_en, row.names = NULL)
  
  top_genes_group = list(bad_cluster = bad_cluster1_gene_en$ENTREZID,good_cluster = good_cluster2_gene_en$ENTREZID)
  
  ck <- compareCluster(geneCluster = top_genes_group, fun = "enrichKEGG")
  ck <- setReadable(ck, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
  
  # for top300 (from heatmap)
  fig_path = paste0(main.path_tc,"/",CancerType,"_analysis_exp_bf/")
  setwd(fig_path)
  
  png(filename = paste0(CancerType,"_cnetplot.png"),
      width = 25, height = 25,  units = "cm" ,pointsize = 12,
      bg = "white", res = 1200, family = "")
  
  exp_cnet = cnetplot(ck) + ggtitle("Exp_for_Top300")
  
  print(exp_cnet)
  dev.off()
  
  png(filename = paste0(CancerType,"_dotplot.png"),
      width = 25, height = 25,  units = "cm" ,pointsize = 12,
      bg = "white", res = 1200, family = "")
  
  exp_dot = dotplot(ck) + ggtitle("Exp_for_Top300")
  
  print(exp_dot)
  dev.off()
  
  # for disease enrichment
  
  edo_bad <- enrichDGN(top_genes_group$bad_cluster)
  edo_good <- enrichDGN(top_genes_group$good_cluster)
  
  edox_bad <- setReadable(edo_bad, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
  edox_good <- setReadable(edo_good, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
  
  # barplot
  png(filename = paste0(CancerType,"_DGN_good_barplot.png"),
      width = 25, height = 25,  units = "cm" ,pointsize = 12,
      bg = "white", res = 1200, family = "")
  
  disease_good = barplot(edo_good, showCategory=25)  + ggtitle("top 25")
  
  print(disease_good)
  dev.off()
  
  png(filename = paste0(CancerType,"_DGN_bad_barplot.png"),
      width = 25, height = 25,  units = "cm" ,pointsize = 12,
      bg = "white", res = 1200, family = "")
  
  disease_bad = barplot(edo_bad, showCategory=25)  + ggtitle("top 25")
  
  print(disease_bad)
  dev.off()
  
  # dotplot
  png(filename = paste0(CancerType,"_DGN_good_dotplot.png"),
      width = 25, height = 25,  units = "cm" ,pointsize = 12,
      bg = "white", res = 1200, family = "")
  
  disease_dot_good = dotplot(edo_good, showCategory=30) + ggtitle("dotplot for good")
  
  print(disease_dot_good)
  dev.off()
  
  png(filename = paste0(CancerType,"_DGN_bad_dotplot.png"),
      width = 25, height = 25,  units = "cm" ,pointsize = 12,
      bg = "white", res = 1200, family = "")
  
  disease_dot_bad =  dotplot(edo_bad, showCategory=30) + ggtitle("dotplot for bad")
  
  print(disease_dot_bad)
  dev.off()
  
  # tree plot
  edox_good <- pairwise_termsim(edo_good)
  png(filename = paste0(CancerType,"_DGN_good_treeplot.png"),
      width = 25, height = 25,  units = "cm" ,pointsize = 12,
      bg = "white", res = 1200, family = "")
  
  disease_tree_good =  treeplot(edox_good) + ggtitle("treeplot for good")
  
  print(disease_tree_good)
  dev.off()
  
  edox_bad <- pairwise_termsim(edo_bad)
  
  png(filename = paste0(CancerType,"_DGN_bad_treeplot.png"),
      width = 25, height = 25,  units = "cm" ,pointsize = 12,
      bg = "white", res = 1200, family = "")
  
  disease_tree_bad =  treeplot(edox_bad) + ggtitle("treeplot for bad")
  
  print(disease_tree_bad)
  dev.off()
  
  # emapplot
  
  png(filename = paste0(CancerType,"_DGN_good_emapplot.png"),
      width = 25, height = 25,  units = "cm" ,pointsize = 12,
      bg = "white", res = 1200, family = "")
  
  disease_emap_good = emapplot(edox_good) + ggtitle("emapplot for good")
  
  print(disease_emap_good)
  dev.off()
  
  png(filename = paste0(CancerType,"_DGN_bad_emapplot.png"),
      width = 25, height = 25,  units = "cm" ,pointsize = 12,
      bg = "white", res = 1200, family = "")
  
  disease_emap_bad = emapplot(edox_bad) + ggtitle("emapplot for bad")
  
  print(disease_emap_bad)
  dev.off()
  
  # upsetplot
  png(filename = paste0(CancerType,"_DGN_good_upsetplot.png"),
      width = 25, height = 25,  units = "cm" ,pointsize = 12,
      bg = "white", res = 1200, family = "")
  
  disease_upset_good = upsetplot(edox_good) + ggtitle("upsetplot for good")
  
  print(disease_upset_good)
  dev.off()
  
  png(filename = paste0(CancerType,"_DGN_bad_upsetplot.png"),
      width = 25, height = 25,  units = "cm" ,pointsize = 12,
      bg = "white", res = 1200, family = "")
  
  disease_upset_bad = upsetplot(edox_bad) + ggtitle("usetplot for bad")
  
  print(disease_upset_bad)
  dev.off()
  
  # for ttest method (bad vs good)
  bad_group$cluster = "bad" 
  good_group$cluster = "good" 
  
  total_group = rbind(good_group,bad_group)
  
  pvals <- as.data.frame(matrix(nrow = c(ncol(total_group)-3)))
  rownames(pvals) = colnames(total_group)[1:(ncol(total_group)-3)]
  colnames(pvals) = "pval"

  # t_test for divide genes by good or bad 
  for(i in 1:c(ncol(total_group)-3)) {
    gene_expr_bad <- bad_group[,i]
    gene_expr_good <- good_group[,i]
  
    ttest_result <- tryCatch({
      t.test(gene_expr_bad, gene_expr_good)
      
    }, error = function(e) {
       NA
    })
    if (sum(is.na(ttest_result)) != 0) {
      pvals[i,"pval"] <- NA
    }else {
      pvals[i,"pval"] <- ttest_result$p.value
    }
    
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
  
  hreason = c(0.1,0.2,0.3,0.4,0.5)
  
  for (not_spe in hreason) {
    if (length(bad_cluster1_gene) == 0 || length(good_cluster2_gene) == 0) {
      
      deg_good_bad = rownames(pvals)[which(pvals < 0.1)]
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
    } else {
      break
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
  
  # for exp enrichment
  fig_path = paste0(main.path_tc,"/",CancerType,"_analysis_exp_bf/")
  setwd(fig_path)
  
  png(filename = paste0(CancerType,"_ttest_cnetplot.png"),
      width = 25, height = 25,  units = "cm" ,pointsize = 12,
      bg = "white", res = 1200, family = "")
  
  exp_cnet = cnetplot(ck) + ggtitle("Exp_for_t_test")
  
  print(exp_cnet)
  dev.off()
  
  png(filename = paste0(CancerType,"_ttest_dotplot.png"),
      width = 25, height = 25,  units = "cm" ,pointsize = 12,
      bg = "white", res = 1200, family = "")
  
  exp_dot = dotplot(ck)
  
  print(exp_dot)
  dev.off()
  
  # for disease enrichment
  
  edo_bad <- enrichDGN(top_genes_group$bad_cluster)
  edo_good <- enrichDGN(top_genes_group$good_cluster)
  
  edox_bad <- setReadable(edo_bad, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
  edox_good <- setReadable(edo_good, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
  
  # barplot
  png(filename = paste0(CancerType,"_ttest_DGN_good_barplot.png"),
      width = 25, height = 25,  units = "cm" ,pointsize = 12,
      bg = "white", res = 1200, family = "")
  
  disease_good = barplot(edo_good, showCategory=25)  + ggtitle("top 25")
  
  print(disease_good)
  dev.off()
  
  png(filename = paste0(CancerType,"_ttest_DGN_bad_barplot.png"),
      width = 25, height = 25,  units = "cm" ,pointsize = 12,
      bg = "white", res = 1200, family = "")
  
  disease_bad = barplot(edo_bad, showCategory=25)  + ggtitle("top 25")
  
  print(disease_bad)
  dev.off()

  # dotplot
  png(filename = paste0(CancerType,"_ttest_DGN_good_dotplot.png"),
      width = 25, height = 25,  units = "cm" ,pointsize = 12,
      bg = "white", res = 1200, family = "")
  
  disease_dot_good = dotplot(edo_good, showCategory=30) + ggtitle("dotplot for good")
  
  print(disease_dot_good)
  dev.off()
  
  png(filename = paste0(CancerType,"_ttest_DGN_bad_dotplot.png"),
      width = 25, height = 25,  units = "cm" ,pointsize = 12,
      bg = "white", res = 1200, family = "")
  
  disease_dot_bad =  dotplot(edo_bad, showCategory=30) + ggtitle("dotplot for bad")
  
  print(disease_dot_bad)
  dev.off()
  
  # tree plot
  edox_good <- pairwise_termsim(edo_good)
  png(filename = paste0(CancerType,"_ttest_DGN_good_treeplot.png"),
      width = 25, height = 25,  units = "cm" ,pointsize = 12,
      bg = "white", res = 1200, family = "")
  
  disease_tree_good =  treeplot(edox_good) + ggtitle("treeplot for good")
  
  print(disease_tree_good)
  dev.off()
  
  edox_bad <- pairwise_termsim(edo_bad)
  
  png(filename = paste0(CancerType,"_ttest_DGN_bad_treeplot.png"),
      width = 25, height = 25,  units = "cm" ,pointsize = 12,
      bg = "white", res = 1200, family = "")
  
  disease_tree_bad =  treeplot(edox_bad) + ggtitle("treeplot for bad")
  
  print(disease_tree_bad)
  dev.off()
  
  # emapplot
  
  png(filename = paste0(CancerType,"_ttest_DGN_good_emapplot.png"),
      width = 25, height = 25,  units = "cm" ,pointsize = 12,
      bg = "white", res = 1200, family = "")
  
  disease_emap_good = emapplot(edox_good) + ggtitle("emapplot for good")
  
  print(disease_emap_good)
  dev.off()
  
  png(filename = paste0(CancerType,"_ttest_DGN_bad_emapplot.png"),
      width = 25, height = 25,  units = "cm" ,pointsize = 12,
      bg = "white", res = 1200, family = "")
  
  disease_emap_bad = emapplot(edox_bad) + ggtitle("emapplot for bad")
  
  print(disease_emap_bad)
  dev.off()
  
  # upsetplot
  png(filename = paste0(CancerType,"_ttest_DGN_good_upsetplot.png"),
      width = 25, height = 25,  units = "cm" ,pointsize = 12,
      bg = "white", res = 1200, family = "")
  
  disease_upset_good = upsetplot(edox_good) + ggtitle("upsetplot for good")
  
  print(disease_upset_good)
  dev.off()
  
  png(filename = paste0(CancerType,"_ttest_DGN_bad_upsetplot.png"),
      width = 25, height = 25,  units = "cm" ,pointsize = 12,
      bg = "white", res = 1200, family = "")
  
  disease_upset_bad = upsetplot(edox_bad) + ggtitle("usetplot for bad")
  
  print(disease_upset_bad)
  dev.off()
  
  remove(exp_cnet,fig_path,ck,top_genes_group,good_cluster2_gene_en,bad_cluster1_gene_en,good_cluster2_gene,bad_cluster1_gene,bad_out
         ,good_out,bad_group,good_group,out,tmp_pheat_cut,bf_zs_shared_exp_filt_df,best_features_df,bf_sh_zs_exp_wo_filt)
}  


