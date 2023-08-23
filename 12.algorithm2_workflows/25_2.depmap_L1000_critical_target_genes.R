library(tidyverse)
library(ggplot2)
library(ggridges)
library(reshape2)
library(readxl)
library(org.Hs.eg.db)
library(cmapR)
library(GSVA)
library(ggsignif)
library(ggrepel)

filepath = "/home/seokwon/nas/"
ref_path = paste0(filepath, "99.reference/")

link_genes = readRDS(paste0(ref_path, "/KEGG_pathway_shared_genes.rds"))
single_genes = readRDS(paste0(ref_path, "/Kegg_pathway_genes.rds"))
link_genes$n_genes = NULL
link_genes_filtered = link_genes[which(link_genes$shared_genes != ""),]

link_genes_filtered <- link_genes_filtered %>%
  mutate(shared_genes = strsplit(shared_genes, ",")) %>%
  unnest(cols = shared_genes)
link_genes_filtered_df = as.data.frame(link_genes_filtered)
colnames(link_genes_filtered_df) = colnames(single_genes)


cancer_specific_critical_features = read_csv("~/nas/04.Results/critical_features/unique_short_long_critical_features_for_SW.csv")
meta_for_TCGA = read.csv(paste0(ref_path, "/Depmap_meta_filt_to_TCGA.csv"))
meta_for_cellline = read.csv("~/nas/99.reference/meta_L1000_matched_TCGA.csv")
LINCS_meta <- read.table("~/nas/99.reference/GSE92742_Broad_LINCS_inst_info.txt", sep="\t", header = T)

# for fic
fig_path = "~/nas/04.Results/drug/L1000/"
if(!dir.exists(fig_path)){
  dir.create(fig_path)
  print(paste0("Created folder: ", fig_path))
} else {
  print(paste0("Folder already exists: ", fig_path))
}
setwd(fig_path)

## Read Dataset 
meta.tbl <- read.csv(file = paste0(ref_path,'DepMap/Model.csv'), 
                     sep=',', header=T, fill=T)

ge.tbl <- read.csv(file = paste0(ref_path, 'DepMap/CRISPRGeneEffect.csv'), 
                   sep=',', header=T, fill=T, row.names=1)
# rnai
# ge.tbl <- read.csv(file = '/mnt/gluster_server/data/raw/DepMap/23Q2/', 
#                    sep=',', header=T, fill=T, row.names=1)
# num_CancerType = "18.TCGA-LUAD"
Cancerlist = dir(paste0(filepath, "/00.data/filtered_TCGA/"))

for (num_CancerType in Cancerlist) {
  
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  ModelType = gsub('TCGA-','', CancerType)
  
  # call input
  
  if (sum(cancer_specific_critical_features$which_cancer == CancerType) == 0) {
    next
  } else {
    
    cancer_critical_features = cancer_specific_critical_features[which(cancer_specific_critical_features$which_cancer == CancerType),]
    cancer_critical_features = cancer_critical_features[which(cancer_critical_features[,CancerType] != "common"),]
    total_link_genes = link_genes_filtered_df[which(link_genes_filtered_df$Pathway %in% cancer_critical_features$features),]$Genes
    total_single_genes = single_genes[which(single_genes$Pathway %in% cancer_critical_features$features),]$Genes
    interested_gene_set = c(total_link_genes,total_single_genes)
  }
  
  ## Interested gene set
  depmap_common_genes <- colnames(ge.tbl)[which(colnames(ge.tbl) %in% interested_gene_set)]

  ## Select cells you want
  cell_oi <- meta.tbl %>% 
    #subset(GrowthPattern != 'Organoid' & DepmapModelType == 'COAD') %>%
    #subset(GrowthPattern != 'Organoid' & DepmapModelType == 'COAD' & PrimaryOrMetastasis == 'Primary') %>%
    # subset(GrowthPattern != 'Organoid' & DepmapModelType == 'BLCA' & PrimaryOrMetastasis == 'Primary' & SampleCollectionSite == 'Colon') %>%
    subset(GrowthPattern != 'Organoid' & 
             PrimaryOrMetastasis == 'Primary' &
             OncotreeLineage == unique(meta_for_TCGA[which(meta_for_TCGA$TCGA == CancerType),]$oncotree) &
             OncotreeCode %in% unique(meta_for_TCGA[which(meta_for_TCGA$TCGA == CancerType),]$oncotreecode)) %>%
    dplyr::select(ModelID) %>% unlist(use.names=F)
  
  ## Subset dataset
  ge.sub <- ge.tbl[rownames(ge.tbl) %in% cell_oi, c('cell_id', depmap_common_genes)] 
  ge.interest <- melt(ge.sub, id.vars = 'cell_id')
  
  depmap_sig_genes = as.vector(unique(ge.interest[which(ge.interest$value < -1),]$variable))
  
  # 
  fig_path = "~/nas/04.Results/drug/depmap/"
  setwd(fig_path)
  
  png(filename = paste0(CancerType,"_depmap_specific_critical_features.png"),
      width = 35, height = 35,  units = "cm" ,pointsize = 12,
      bg = "white", res = 1200, family = "")
  
  ## Kenel Density Curve
  total_out = ggplot(ge.interest, aes(x = value, y = variable, fill = after_stat(x))) +
    geom_density_ridges_gradient(scale = 2, gradient_lwd = .5, color = "black",
                                 jittered_points = TRUE,
                                 position = position_points_jitter(width = 0.05, height = 0),
                                 point_shape = '|', point_size = 3, point_alpha = 1, alpha = 0.7) +
    scale_fill_viridis_c(option = "plasma", name = "Gene Effect") +
    geom_vline(xintercept = -1, color='darkred')+
    geom_vline(xintercept = c(1,0,-2), color = 'darkgrey') +
    geom_vline(xintercept = c(-2.5, -1.5, -0.5, 0.5), color = 'darkgrey', linetype = 'dashed') +
    
    ylab('Gene') +
    xlab('Gene Effect \n(Chronos score, 23Q2)') +
    theme_bw() 
  
  print(total_out)
  
  dev.off()
  
  #################################################################################################
  #################################################################################################
  #### 2. L1000 Plot for target genes                                                       #######
  #################################################################################################
  #################################################################################################
  
  # filtered cell line by CancerType
  cell_id = meta_for_cellline[meta_for_cellline$cancer_type == CancerType,]$cell_line
  if (length(cell_id) == 0) {
    next
  }
  l1000_meta <- LINCS_meta[LINCS_meta$cell_id %in% cell_id,] # LOVO SW480 SW620 HCT116  # change a cancer cell line 
  # 
  DMSO_meta <- l1000_meta %>% dplyr::filter( pert_iname == "DMSO" ) # CTRL
  # # USP7_samples_meta <- l1000_meta %>% dplyr::filter( pert_type == "trt_sh" &  pert_iname == "USP7") # trt_sh : Knock down, pert_iname : gene_name
  candi_samples_meta <- l1000_meta %>% dplyr::filter( pert_type == "trt_sh" )
  
  if (nrow(candi_samples_meta) == 0) {
    next
  }
  # candi_samples_meta <- l1000_meta %>% dplyr::filter(pert_iname %in% interested_gene_set)
  candi_merged_meta <- rbind(DMSO_meta, candi_samples_meta)

  L1000_exp <- cmapR::parse_gctx("/mnt/gluster_server/dkshin/data/LINCS/GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx", cid = candi_merged_meta$inst_id)
  # L1000_zscore_rep <- parse_gctx("~/nas/LINCS/cmapZscores.gctx", cid = l1000_meta$inst_id)

  l1000_rownames <- L1000_exp@rid

  l1000_gene_symbol <- AnnotationDbi::select(org.Hs.eg.db, keys=l1000_rownames, columns="SYMBOL", keytype="ENTREZID")
  l1000_gene_symbol <- l1000_gene_symbol$SYMBOL
  if (l1000_rownames %>% length != l1000_gene_symbol %>% length) {
    print("Something is wrong")
    break
  }

  merged_exp <- L1000_exp@mat
  rownames(merged_exp) <- l1000_gene_symbol

  signatures = list( interested_gene_set)
  
  gs_res <- gsva(merged_exp, signatures, parallel.sz = 50L)
  
  if (all.equal(rownames(t(gs_res)),candi_merged_meta$inst_id) ) {
    gs_res_t = as.data.frame(t(gs_res))
    gs_res_t$gene_id = candi_merged_meta$pert_iname
  }

  per.gene.value <- data.frame(gs_value = gs_res_t$V1, 
                               pert_gene = gs_res_t$gene_id)
  
  per_gene_value_edit = per.gene.value[which(per.gene.value$pert_gene %in% c("DMSO",interested_gene_set)),]
  
  for (ck_gene in depmap_sig_genes) {
    if (nrow(per_gene_value_edit[which(per_gene_value_edit$pert_gene == ck_gene),]) < 2) {
      depmap_sig_genes = depmap_sig_genes[!depmap_sig_genes %in% ck_gene]
      
    }
    
  }
  
  depmap_list <- lapply(depmap_sig_genes, function(gene) c("DMSO", gene))
  
  fig_path = "~/nas/04.Results/drug/L1000/"
  setwd(fig_path)
  # A really basic boxplot.
  png(filename = paste0(CancerType,"_l1000_specific_critical_genes.png"),
      width = 35, height = 35,  units = "cm" ,pointsize = 12,
      bg = "white", res = 1200, family = "")
  
  l1000_out = ggplot(per_gene_value_edit, aes(x=pert_gene, y=gs_value)) + 
    geom_boxplot(fill="orange", alpha=0.5) + 
    geom_signif(comparisons = depmap_list, 
                map_signif_level = TRUE, step_increase=0.1) +
    xlab("Pert_genes")
  
  print(l1000_out)
  
  dev.off()
  
  
  #################################################################################################
  #################################################################################################
  #### 3. (L1000 Signature score vs. DepMap effect score) Plot                              #######
  #################################################################################################
  #################################################################################################

  final.genes = intersect(unique(per_gene_value_edit$pert_gene), depmap_common_genes)
  
  final.genes_edit = final.genes
  for (test_gene in final.genes) {
    if (nrow(per_gene_value_edit[which(per_gene_value_edit$pert_gene == test_gene),]) == 1) {
      final.genes_edit = final.genes_edit[! final.genes_edit %in% test_gene] 
    }
  }
  
  gene.filter.plot <- as.data.frame(matrix(nrow = length(final.genes_edit), ncol = 5))
  colnames(gene.filter.plot) <- c("Genes", "mean_depmap", "mean_eff_depmap", "pval_l1000", "gs_l1000")
  
  gene.filter.plot$Genes <- final.genes_edit
  # Mean(depmap score)
  gene.filter.plot$mean_depmap <- sapply(1:length(final.genes_edit), 
                                         function(x) -mean(dplyr::filter(ge.interest, variable == final.genes_edit[x])$value) ) 
  
  # Mean(depmap score <-1)
  temp.depmap <-sapply(1:length(final.genes_edit), 
                       function(x) -mean(dplyr::filter(dplyr::filter(ge.interest, variable == final.genes_edit[x]), value <= -1)$value) ) 
  temp.depmap[is.na(temp.depmap)] <- 0
  gene.filter.plot$mean_eff_depmap <- temp.depmap
  
  # L1000 signature score
  gene.filter.plot$gs_l1000 <- sapply(1:length(final.genes_edit),
                                      function(x) mean(dplyr::filter(per_gene_value_edit, pert_gene == final.genes_edit[x])$gs_value)) 
  # L1000 p-value
  
  gene.filter.plot$pval_l1000 <- sapply(1:length(final.genes_edit),
                                        function(x) -log10(t.test(dplyr::filter(per_gene_value_edit, pert_gene == final.genes_edit[x])$gs_value, 
                                                                  y=dplyr::filter(per_gene_value_edit, pert_gene == "DMSO")$gs_value, alternatives="greater")$p.value)) 
  gene.filter.plot$sig <- ifelse(gene.filter.plot$pval_l1000 > -log10(0.05) & gene.filter.plot$mean_eff_depmap !=0, "Dual_sig", "Not Sig")
  
  ## Plots
  
  # ggplot(gene.filter.plot, aes(mean_depmap, pval_l1000, color = sig)) +
  #   scale_color_manual(values = c("red","grey")) +
  #   geom_text_repel(
  #     data = subset(gene.filter.plot, pval_l1000 > -log10(0.05)),
  #     aes(label = Genes),
  #     size = 5,
  #     box.padding = unit(0.35, "lines"),
  #     point.padding = unit(0.3, "lines")) +
  #   geom_point()+
  #   theme_bw(base_size = 12)
  
  fig_path = "~/nas/04.Results/drug/"
  setwd(fig_path)
  # A really basic boxplot.
  png(filename = paste0(CancerType,"_depmap_l1000_specific_critical_genes.png"),
      width = 35, height = 35,  units = "cm" ,pointsize = 12,
      bg = "white", res = 1200, family = "")
  
  dual_out = ggplot(gene.filter.plot, aes(mean_depmap, pval_l1000, color = sig)) +
    scale_color_manual(values = c("red","grey")) +
    geom_text_repel(
      # data = subset(gene.filter.plot, pval_l1000 > -log10(0.05)),
      aes(label = Genes),
      size = 5,
      box.padding = unit(0.35, "lines"),
      point.padding = unit(0.3, "lines")) +
    geom_point()+
    theme_bw(base_size = 12)+
    ggtitle(paste0(CancerType, "Depmap & L1000 significant target genes"))
  
  
  print(dual_out)
  
  dev.off()
  
  

  
}
