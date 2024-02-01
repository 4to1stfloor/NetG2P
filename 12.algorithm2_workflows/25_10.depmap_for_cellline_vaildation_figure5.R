library(ggbiplot)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggfortify)
library(pheatmap)
library(RColorBrewer)
library(Rtsne)
# library(tsne)
# library(umap)
library(tidyr)
library(dplyr)
library(reshape2)
library(readxl)
library(ggridges)
library(ggdist)
library(ggpubr)

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

## Read Dataset 
meta.tbl <- read.csv(file = paste0(ref_path,'DepMap/Model.csv'), 
                     sep=',', header=T, fill=T)

ge.tbl <- read.csv(file = paste0(ref_path, 'DepMap/CRISPRGeneEffect.csv'), 
                   sep=',', header=T, fill=T, row.names=1)

meta_for_TCGA = read.csv(paste0(ref_path, "/Depmap_meta_filt_to_TCGA.csv"))

Cancerlist = dir(paste0(filepath, "/00.data/filtered_TCGA/"))
setwd("~/nas/04.Results/drug/depmap/")
####
Cancerlist = Cancerlist[c(-11,-12)]
# num_CancerType = "32.TCGA-UCEC"

for (num_CancerType in Cancerlist) {
  
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  Cancername = gsub('TCGA-' , '', CancerType)
  
  # call input
  gc_cellline = readRDS(paste0("~/nas/00.data/filtered_TCGA/", num_CancerType, "/",Cancername,"_cellline_dual_all_log.rds"))
  gc_TCGA = readRDS(paste0("~/nas/04.Results/short_long/", CancerType,"_critical_features_short_long.rds"))
  
  short_cluster = gc_TCGA %>% filter(cluster == "short") %>% dplyr::select(-vitalstatus, -duration,-status, -cluster)
  long_cluster = gc_TCGA %>% filter(cluster == "long") %>% dplyr::select(-vitalstatus, -duration,-status, -cluster)
  wo_info = gc_TCGA %>% dplyr::select(-vitalstatus, -duration, -status)
  
  gc_cellline$cluster = "cell"
  gc_cellline_filt_df = gc_cellline[,colnames(wo_info)]

  ### Dimension reduction (PCA)
  
  set.seed(1)

  combine_df = rbind(wo_info, gc_cellline_filt_df)
  
  pca_res_combine <- prcomp(combine_df %>% dplyr::select(-cluster), retx=TRUE, center=TRUE, scale=TRUE)
  library(broom)
  variance_exp <- pca_res_combine %>%  
    tidy("pcs") %>% 
    pull(percent)
  
  pca_plot = pca_res_combine %>%
    ggplot(aes(x = PC1, y = PC2, color= combine_df$cluster))+
    geom_point()+
    labs(x = paste0("PC1: ",round(variance_exp[1]*100), "%"),
         y = paste0("PC2: ",round(variance_exp[2]*100), "%")) +
    theme_classic()+  
    scale_color_manual(values=c("#6633CC", "#4DAF4A","#E41A1C")) 
  
  ggsave(file=paste0(CancerType , "_PCA_plot.svg"), plot=pca_plot, width=10, height=10)
  
  ## distance calculate
  
  pca_res = as.data.frame(pca_res_combine$x)
  
  if (all.equal(rownames(pca_res) , c(rownames(wo_info), rownames(gc_cellline_filt_df)))) {
    pca_res$cluster = c(wo_info$cluster, gc_cellline_filt_df$cluster)
    pca_long = pca_res %>% filter(cluster == "long")
    pca_short = pca_res %>% filter(cluster == "short")
    pca_cell = pca_res %>% filter(cluster == "cell")
    
  }
  
  long_xy = c(mean(pca_long[,1]), mean(pca_long[,2]))
  short_xy = c(mean(pca_short[,1]), mean(pca_short[,2]))
  
  tmp_cluster =c()
  for (num_cell in 1:nrow(pca_cell)) {
    tmp_long = sqrt((long_xy[1] - pca_cell[num_cell,1])^2 + (long_xy[2] - pca_cell[num_cell,2])^2)
    tmp_short = sqrt((short_xy[1] - pca_cell[num_cell,1])^2 + (short_xy[2] - pca_cell[num_cell,2])^2)
    if (tmp_long > tmp_short) {
      tmp_cluster = c(tmp_cluster, "short")
    } else {
      tmp_cluster = c(tmp_cluster, "long")
    }
  }
  
  gc_cellline_filt_df$cluster = tmp_cluster
  
  ### heatmap 
  cancer_bf = read.csv(paste0("~/nas/04.Results/short_long/ttest_common/",CancerType,"_critical_features_short_long_common.csv"))
  
  total_cell_for_fig = gc_cellline_filt_df
  total_cell_meta = data.frame(cluster = total_cell_for_fig[,c("cluster")])
  total_cell_for_tmp = total_cell_for_fig[,which(!colnames(total_cell_for_fig) %in% c("cluster"))]
  
  # # Convert the matrix to a numeric matrix
  
  tmp_numeric <- matrix(as.numeric(unlist(total_cell_for_tmp)), nrow = nrow(total_cell_for_tmp))
  vec <- as.vector(tmp_numeric)
  # hist(vec)
  total_cell_for_tmp[total_cell_for_tmp > quantile(vec, probs = 0.75) + 1.5*IQR(vec)] <- round(quantile(vec, probs = 0.75) + 1.5*IQR(vec),1)
  
  total_cell_for_fig_final = cbind(total_cell_for_tmp,total_cell_meta)
  
  annotation_col <- data.frame(patients_group= total_cell_for_fig_final$cluster)
  rownames(annotation_col) <- rownames(total_cell_for_fig_final)
  
  num_features = c("short" = "#E41A1C", "long" = "#4DAF4A" , "common" = "#377EB8")
  short_long_colors <- c("short" = "red", "long" = "#009E73")
  
  annotation_row = data.frame(types = cancer_bf$classification)
  rownames(annotation_row) <- cancer_bf$variable
  
  ann_colors_sl = list(patients_group = short_long_colors, types = num_features )
  Colors = brewer.pal(9, "YlOrRd")
  
  png(filename = paste0(CancerType , "_cellline_heatmap_plot.png"),
      width = 10, height = 10,  units = "cm" ,pointsize = 12,
      bg = "white", res = 900, family = "")
  
  cell_heatmap = ComplexHeatmap::pheatmap(as.matrix(t(total_cell_for_fig_final %>% dplyr::select(-cluster))),
                                          column_split = factor(annotation_col$patients_group, levels = c("short","long")),
                                          annotation_col = annotation_col,
                                          annotation_row = annotation_row,
                                          annotation_colors = ann_colors_sl,
                                          cluster_rows = T,
                                          cluster_cols = T,
                                          legend = T,
                                          annotation_legend = T,
                                          show_colnames = T,
                                          cluster_column_slices = FALSE,
                                          color = Colors,
                                          show_rownames = F
  )
  print(cell_heatmap)
  
  dev.off()

  #### validation 
  cellline_long = gc_cellline_filt_df %>% filter(cluster == "long")
  cellline_short = gc_cellline_filt_df %>% filter(cluster == "short")
  
  depmap_common_link_genes <- unique(link_genes_filtered_df[which(link_genes_filtered_df$Pathway %in% colnames(gc_cellline_filt_df %>% dplyr::select(-cluster))),]$Genes)
  depmap_common_each_genes <- unique(link_genes_filtered_df[which(single_genes$Pathway %in% colnames(gc_cellline_filt_df %>% dplyr::select(-cluster))),]$Genes)
  depmap_common_genes = unique(c(depmap_common_link_genes,depmap_common_each_genes))
  
  filtered_df = data.frame()
  for (critical_genes in depmap_common_genes) {
    
    ge_sub_short = ge.tbl %>% 
      filter(rownames(ge.tbl) %in% rownames(cellline_short)) %>%
      dplyr::select(cell_id,any_of(critical_genes))
    
    ge_sub_long = ge.tbl %>% 
      filter(rownames(ge.tbl) %in% rownames(cellline_long)) %>%
      dplyr::select(cell_id,any_of(critical_genes))
    
    ge_short_interest <- melt(ge_sub_short, id.vars = 'cell_id')
    ge_short_interest = ge_short_interest[!is.na(ge_short_interest$value),]
    ge_long_interest <- melt(ge_sub_long, id.vars = 'cell_id')
    ge_long_interest = ge_long_interest[!is.na(ge_long_interest$value),]
    
    if (!is.character(ge_long_interest) & !is.logical(ge_short_interest)) {
      if (nrow(ge_long_interest) <=1 | length(ge_long_interest) <=1 ) {
        next
      } else if (nrow(ge_short_interest) <=1 | length(ge_short_interest) <=1) {
        next
      } else {
        t_test_res =  t.test(ge_short_interest$value, ge_long_interest$value)
        
        if (t_test_res$p.value <0.05) {
          tmp_test_df = data.frame(genes = critical_genes,
                                   delta_long_to_short =  mean(ge_short_interest$value) - mean(ge_long_interest$value),
                                   sig = "significant")
          
          filtered_df = rbind(filtered_df, tmp_test_df)
          
        } else {
          tmp_test_df = data.frame(genes = critical_genes,
                                   delta_long_to_short = mean(ge_short_interest$value) - mean(ge_long_interest$value),
                                   sig = "not")
          
          filtered_df = rbind(filtered_df, tmp_test_df)
        }
      }
    } else {
      next
    }
    
  }
  
  if (nrow(filtered_df) == 0) {
    next
  } else {
    filtered_ordered_df = filtered_df %>% arrange(delta_long_to_short)
  }
  
  dotplot_line = ggplot(filtered_ordered_df, aes(x = factor(genes , levels = filtered_ordered_df$genes), y = delta_long_to_short )) + 
    geom_point(aes(color = sig), 
               alpha = 1.0, size = 1) +
    geom_hline(yintercept = -0.15, colour = "red") + 
    geom_hline(yintercept = 0.15, colour = "red") +
    scale_color_manual(values=c("#999999", "#E69F00")) +
    theme_classic() +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  
  ggsave(file=paste0(CancerType , "_cellline_dotplotline_plot.svg"), plot=dotplot_line, width=10, height=10)
  
  difference_df = filtered_ordered_df[which(abs(filtered_ordered_df$delta_long_to_short) > 0.15),]
  sig_dif_df = difference_df %>% filter(sig == "significant")
  
  if (nrow(sig_dif_df) == 0) {
    next
  } else {
    
    sig_sub_short = ge.tbl %>% 
      filter(rownames(ge.tbl) %in% rownames(cellline_short)) %>%
      dplyr::select(cell_id,any_of(sig_dif_df$genes))
    sig_sub_long = ge.tbl %>% 
      filter(rownames(ge.tbl) %in% rownames(cellline_long)) %>%
      dplyr::select(cell_id,any_of(sig_dif_df$genes))
    
    sig_short_interest <- melt(sig_sub_short, id.vars = 'cell_id')
    sig_short_interest = sig_short_interest[!is.na(sig_short_interest$value),]
    sig_long_interest <- melt(sig_sub_long, id.vars = 'cell_id')
    sig_long_interest = sig_long_interest[!is.na(sig_long_interest$value),]
    
    if (length(sig_long_interest) == 0 | length(sig_long_interest) == 0) {
      
    }
    gene.filter_short_plot <- as.data.frame(matrix(nrow = length(sig_dif_df$genes), ncol = 2))
    gene.filter_long_plot <- as.data.frame(matrix(nrow = length(sig_dif_df$genes), ncol = 2))
    
    colnames(gene.filter_short_plot) <- c("Genes", "mean_depmap")
    colnames(gene.filter_long_plot) <- c("Genes", "mean_depmap")
    
    gene.filter_short_plot$Genes <- sig_dif_df$genes
    gene.filter_long_plot$Genes <- sig_dif_df$genes
    # Mean(depmap score)
    gene.filter_short_plot$mean_depmap <- sapply(1:length(sig_dif_df$genes),
                                                 function(x) -mean(dplyr::filter(sig_short_interest, variable == sig_dif_df$genes[x])$value) )
    
    gene.filter_long_plot$mean_depmap <- sapply(1:length(sig_dif_df$genes),
                                                function(x) -mean(dplyr::filter(sig_long_interest, variable == sig_dif_df$genes[x])$value) )
    
    gene_filtered_short_edit = gene.filter_short_plot %>% filter(mean_depmap > 0.5)
    gene_filtered_short_low = gene.filter_short_plot %>% filter(mean_depmap < 0.5)
    
    gene_filtered_long_edit = gene.filter_long_plot %>% filter(mean_depmap > 0.5)
    gene_filtered_long_low = gene.filter_long_plot %>% filter(mean_depmap < 0.5)
    
    ### both > 0.5 
    sig_1_genes = intersect(gene_filtered_short_edit$Genes , gene_filtered_long_edit$Genes)
    
    if (length(sig_1_genes) == 0 ) {
      next
    } else {
      
      gene_filtered_short_edit = gene_filtered_short_edit %>% filter(Genes %in% sig_1_genes)
      gene_filtered_long_edit = gene_filtered_long_edit %>% filter(Genes %in% sig_1_genes)
      
      sig_short_edit_specific = sig_short_interest[which(sig_short_interest$variable %in% gene_filtered_short_edit$Genes),]
      sig_long_edit_specific = sig_long_interest[which(sig_long_interest$variable %in% gene_filtered_long_edit$Genes),]
      
      sig_short_edit_specific$group = "short"
      sig_long_edit_specific$group = "long"
      
      sig_total_edit_specific = rbind(sig_short_edit_specific,sig_long_edit_specific)
      
      kernal_plot = sig_total_edit_specific %>%
        ggdensity(x = "value",
                  # add = "mean",
                  rug = T,
                  color = "group",
                  fill = "group",
                  palette = c("#4DAF4A", "#E41A1C")) +
        scale_x_continuous(limits = c(-3.5, 0.5)) +
        facet_grid(rows = "variable", scales = "free", space = "free")
      
      ggsave(file=paste0(CancerType , "_kernel_both_plot.svg"), plot=kernal_plot, width=10, height=10)
      
    }
    ### short > 0.5 & long < 0.5 
    sig_2_genes = intersect(gene_filtered_short_edit$Genes , gene_filtered_long_low$Genes)
    
    if (length(sig_2_genes) == 0 ) {
      next
    } else {
      gene_filtered_short_edit = gene_filtered_short_edit %>% filter(Genes %in% sig_2_genes)
      gene_filtered_long_low = gene_filtered_long_low %>% filter(Genes %in% sig_2_genes)
      
      sig_short_edit_specific = sig_short_interest[which(sig_short_interest$variable %in% gene_filtered_short_edit$Genes),]
      sig_long_edit_specific = sig_long_interest[which(sig_long_interest$variable %in% gene_filtered_long_low$Genes),]
      
      sig_short_edit_specific$group = "short"
      sig_long_edit_specific$group = "long"
      
      sig_total_edit_specific = rbind(sig_short_edit_specific,sig_long_edit_specific)
      
      kernal_plot = sig_total_edit_specific %>%
        ggdensity(x = "value",
                  # add = "mean",
                  rug = T,
                  color = "group",
                  fill = "group",
                  palette = c("#4DAF4A", "#E41A1C")) +
        scale_x_continuous(limits = c(-3.5, 0.5)) +
        facet_grid(rows = "variable", scales = "free", space = "free")
      
      ggsave(file=paste0(CancerType , "_kernel_short_plot.svg"), plot=kernal_plot, width=10, height=10)
    }
    ### short < 0.5 & long > 0.5 
    sig_3_genes = intersect(gene_filtered_short_low$Genes , gene_filtered_long_edit$Genes)
    
    if (length(sig_3_genes) == 0 ) {
      next
    } else {
      gene_filtered_short_low = gene_filtered_short_low %>% filter(Genes %in% sig_3_genes)
      gene_filtered_long_edit = gene_filtered_long_edit %>% filter(Genes %in% sig_3_genes)
      
      sig_short_edit_specific = sig_short_interest[which(sig_short_interest$variable %in% gene_filtered_short_low$Genes),]
      sig_long_edit_specific = sig_long_interest[which(sig_long_interest$variable %in% gene_filtered_long_edit$Genes),]
      
      sig_short_edit_specific$group = "short"
      sig_long_edit_specific$group = "long"
      
      sig_total_edit_specific = rbind(sig_short_edit_specific,sig_long_edit_specific)
      
      kernal_plot = sig_total_edit_specific %>%
        ggdensity(x = "value",
                  # add = "mean",
                  rug = T,
                  color = "group",
                  fill = "group",
                  palette = c("#4DAF4A", "#E41A1C")) +
        scale_x_continuous(limits = c(-3.5, 0.5)) +
        facet_grid(rows = "variable", scales = "free", space = "free")
      
      ggsave(file=paste0(CancerType , "_kernel_long_plot.svg"), plot=kernal_plot, width=10, height=10)
    }
    ### short < 0.5 & long < 0.5
    sig_4_genes = intersect(gene_filtered_short_low$Genes , gene_filtered_long_low$Genes)
    
    if (length(sig_4_genes) == 0 ) {
      next
    } else {
      gene_filtered_short_low = gene_filtered_short_low %>% filter(Genes %in% sig_4_genes)
      gene_filtered_long_low = gene_filtered_long_low %>% filter(Genes %in% sig_4_genes)
      
      sig_short_edit_specific = sig_short_interest[which(sig_short_interest$variable %in% gene_filtered_short_low$Genes),]
      sig_long_edit_specific = sig_long_interest[which(sig_long_interest$variable %in% gene_filtered_long_low$Genes),]
      
      sig_short_edit_specific$group = "short"
      sig_long_edit_specific$group = "long"
      
      sig_total_edit_specific = rbind(sig_short_edit_specific,sig_long_edit_specific)
      
      kernal_plot = sig_total_edit_specific %>%
        ggdensity(x = "value",
                  # add = "mean",
                  rug = T,
                  color = "group",
                  fill = "group",
                  palette = c("#4DAF4A", "#E41A1C")) +
        scale_x_continuous(limits = c(-3.5, 0.5)) +
        facet_grid(rows = "variable", scales = "free", space = "free")
      
      ggsave(file=paste0(CancerType , "_kernel_not_plot.svg"), plot=kernal_plot, width=10, height=10)
    }
    
  }
  
}


