
ucec_cellline = readRDS("~/nas/00.data/filtered_TCGA/32.TCGA-UCEC/UCEC_cellline_dual_all_log.rds")
ucec_TCGA = readRDS("~/nas/04.Results/short_long/TCGA-UCEC_critical_features_short_long.rds")

short_cluster = ucec_TCGA %>% filter(cluster == "short") %>% select(-vitalstatus, -duration,-status, -cluster)
long_cluster = ucec_TCGA %>% filter(cluster == "long") %>% select(-vitalstatus, -duration,-status, -cluster)

short_cluster %>%
  rownames_to_column(var = "features") %>%
  gather(key = "patient", value = "value", -features)

library(ggbiplot)

test_df = ucec_TCGA %>% select(-vitalstatus, -duration, -status)
ucec_cellline$cluster = "cell"
ucec_cellline_filt_df = ucec_cellline[,colnames(test_df)]

combine_df = rbind(test_df, ucec_cellline_filt_df)

pca_data = prcomp(combine_df %>% select(-cluster),
                  center = T,
                  scale. = T)

autoplot(pca_data, 
         data = combine_df, 
         colour = 'cluster'
         # frame = TRUE, 
         # frame.type = 'norm'
)

pca_res = as.data.frame(pca_data$x)

rownames(pca_res) == c(rownames(test_df), rownames(ucec_cellline_filt_df))
pca_res$cluster = c(test_df$cluster, ucec_cellline_filt_df$cluster)

pca_long = pca_res %>% filter(cluster == "long")
pca_short = pca_res %>% filter(cluster == "short")
pca_cell = pca_res %>% filter(cluster == "cell")

long_xy = c(mean(pca_long[,1]), mean(pca_long[,2]))
short_xy = c(mean(pca_short[,1]), mean(pca_short[,2]))

tmp_test =c()
for (num_cell in 1:nrow(pca_cell)) {
  tmp_long = sqrt((long_xy[1] - pca_cell[num_cell,1])^2 + (long_xy[2] - pca_cell[num_cell,2])^2)
  tmp_short = sqrt((short_xy[1] - pca_cell[num_cell,1])^2 + (short_xy[2] - pca_cell[num_cell,2])^2)
  if (tmp_long > tmp_short) {
    tmp_test = c(tmp_test, "short")
  } else {
    tmp_test = c(tmp_test, "long")
  }
}

ucec_cellline_filt_df$cluster = tmp_test

pheatmap(ucec_cellline_filt_df %>% select(-cluster))


cancer_bf = read.csv(paste0("~/nas/04.Results/short_long/ttest_common/",CancerType,"_critical_features_short_long_common.csv"))

annotation_col <- data.frame(patients_group= ucec_cellline_filt_df$cluster)
rownames(annotation_col) <- rownames(ucec_cellline_filt_df)

num_features = c("short" = "#E41A1C", "long" = "#4DAF4A" , common = "#377EB8")
short_long_colors <- c("short" = "red", "long" = "#009E73")

annotation_row = data.frame(types = cancer_bf$classification)
rownames(annotation_row) <- cancer_bf$variable

ann_colors_sl = list(patients_group = short_long_colors, types = num_features )
Colors = brewer.pal(9, "YlOrRd")
ComplexHeatmap::pheatmap(as.matrix(t(ucec_cellline_filt_df %>% select(-cluster))),
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




library(tidyr)
library(dplyr)
library(reshape2)
library(readxl)
library(ggridges)
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

ucec_cellline_long = ucec_cellline_filt_df %>% filter(cluster == "long")
ucec_cellline_short = ucec_cellline_filt_df %>% filter(cluster == "short")

filtered_df = data.frame()
for (critical_features in colnames(ucec_cellline_filt_df %>% select(-cluster))) {
  
  depmap_common_link_genes <- unique(link_genes_filtered_df[which(link_genes_filtered_df$Pathway %in% critical_features),]$Genes)
  depmap_common_each_genes <- unique(link_genes_filtered_df[which(single_genes$Pathway %in% critical_features),]$Genes)
  
  depmap_common_genes = c(depmap_common_link_genes,depmap_common_each_genes)
  
  # cell_oi <- meta.tbl %>% 
  #   #subset(GrowthPattern != 'Organoid' & DepmapModelType == 'COAD') %>%
  #   #subset(GrowthPattern != 'Organoid' & DepmapModelType == 'COAD' & PrimaryOrMetastasis == 'Primary') %>%
  #   # subset(GrowthPattern != 'Organoid' & DepmapModelType == 'BLCA' & PrimaryOrMetastasis == 'Primary' & SampleCollectionSite == 'Colon') %>%
  #   subset(GrowthPattern != 'Organoid' & 
  #            PrimaryOrMetastasis == 'Primary' &
  #            OncotreeLineage %in% unique(meta_for_TCGA[which(meta_for_TCGA$TCGA %in% c("TCGA-OV", "TCGA-BLCA", "TCGA-UCEC")),]$oncotree) &
  #            OncotreeCode %in% unique(meta_for_TCGA[which(meta_for_TCGA$TCGA %in% c("TCGA-OV", "TCGA-BLCA", "TCGA-UCEC")),]$oncotreecode)) %>%
  #   dplyr::select(ModelID) %>% unlist(use.names=F)
  
  ## Subset dataset

  ge_sub_short = ge.tbl %>% 
    filter(rownames(ge.tbl) %in% rownames(ucec_cellline_short)) %>%
    select(cell_id,any_of(depmap_common_genes))
  ge_sub_long = ge.tbl %>% 
    filter(rownames(ge.tbl) %in% rownames(ucec_cellline_long)) %>%
    select(cell_id,any_of(depmap_common_genes))
  
  if (ncol(ge_sub_short) <=1 | ncol(ge_sub_long) <=1) {
    next
  } else {
    ge_short_interest <- melt(ge_sub_short, id.vars = 'cell_id')
    ge_short_interest = ge_short_interest[!is.na(ge_short_interest$value),]
    ge_long_interest <- melt(ge_sub_long, id.vars = 'cell_id')
    ge_long_interest = ge_long_interest[!is.na(ge_long_interest$value),]
    
 
    for (common_genes in depmap_common_genes) {
      if (sum(ge_short_interest$variable == common_genes) != 0 & sum(ge_long_interest$variable == common_genes) != 0 ) {
        tmp_short = ge_short_interest %>% filter(variable == common_genes)
        tmp_long = ge_long_interest %>% filter(variable == common_genes)
        t_test_res =  t.test(tmp_short$value, tmp_long$value)
        
        if (t_test_res$p.value <0.05) {
          tmp_test_df = data.frame(features = critical_features, genes = common_genes)
          filtered_df = rbind(filtered_df, tmp_test_df)
          # print(critical_features)
          # print(common_genes)
          
        }
      } else {
        next
      }
      
    }
  }
  
  
  
}


LAMA5
LAMB4
THBS3
ITGA11
HGF
tmp_genes = filtered_df %>% filter(features == names(which(table(filtered_df$features) == max(table(filtered_df$features)))))

depmap_common_genes = tmp_genes$genes

depmap_common_genes = filtered_df$genes

ge_sub_short = ge.tbl %>% 
  filter(rownames(ge.tbl) %in% rownames(ucec_cellline_short)) %>%
  select(cell_id,any_of(depmap_common_genes))
ge_sub_long = ge.tbl %>% 
  filter(rownames(ge.tbl) %in% rownames(ucec_cellline_long)) %>%
  select(cell_id,any_of(depmap_common_genes))

ge_short_interest <- melt(ge_sub_short, id.vars = 'cell_id')
ge_short_interest = ge_short_interest[!is.na(ge_short_interest$value),]
ge_long_interest <- melt(ge_sub_long, id.vars = 'cell_id')
ge_long_interest = ge_long_interest[!is.na(ge_long_interest$value),]


gene.filter_short_plot <- as.data.frame(matrix(nrow = length(depmap_common_genes), ncol = 2))
gene.filter_long_plot <- as.data.frame(matrix(nrow = length(depmap_common_genes), ncol = 2))

colnames(gene.filter_short_plot) <- c("Genes", "mean_depmap")
colnames(gene.filter_long_plot) <- c("Genes", "mean_depmap")

gene.filter_short_plot$Genes <- depmap_common_genes
gene.filter_long_plot$Genes <- depmap_common_genes
# Mean(depmap score)
gene.filter_short_plot$mean_depmap <- sapply(1:length(depmap_common_genes),
                                             function(x) -mean(dplyr::filter(ge_short_interest, variable == depmap_common_genes[x])$value) )

gene.filter_long_plot$mean_depmap <- sapply(1:length(depmap_common_genes),
                                            function(x) -mean(dplyr::filter(ge_long_interest, variable == depmap_common_genes[x])$value) )

gene_filtered_short_edit = gene.filter_short_plot %>% filter(mean_depmap > 0.5)
gene_filtered_long_edit = gene.filter_long_plot %>% filter(mean_depmap > 0.5)

ge_short_edit_specific = ge_short_interest[which(ge_short_interest$variable %in% gene_filtered_short_edit$Genes),]
ge_long_edit_specific = ge_long_interest[which(ge_long_interest$variable %in% gene_filtered_long_edit$Genes),]

## Kenel Density Curve

ge_short_edit_specific$group = "short"
ge_long_edit_specific$group = "long"

ge_total_edit_specific = rbind(ge_short_edit_specific,ge_long_edit_specific)

library(ggdist)

table(filtered_df$features)

ge_total_edit_specific %>%
  filter(variable != "RFC1") %>%
  ggdensity(x = "value",
            # add = "mean",
            rug = T,
            color = "group",
            fill = "group",
            palette = c("#4DAF4A", "#E41A1C")) +
  scale_x_continuous(limits = c(-3.5, 0.5)) +
  facet_grid(rows = "variable", scales = "free", space = "free")


# unique(ge_total_edit_specific$variable)
# ge_total_edit_specific

# ggplot(ge_short_edit_specific %>% filter(variable != "RFC1"), aes(x = value, y = variable, fill = after_stat(x))) +
#   geom_density_ridges_gradient(scale = 2, gradient_lwd = .5, color = "black",
#                                jittered_points = TRUE,
#                                position = position_points_jitter(width = 0.05, height = 0),
#                                point_shape = '|', point_size = 3, point_alpha = 1, alpha = 0.7) +
#   scale_fill_viridis_c(option = "plasma", name = "Gene Effect") +
#   # geom_vline(xintercept = -1, color='darkred')+
#   geom_vline(xintercept = c(1,0,-2), color = 'darkgrey') +
#   geom_vline(xintercept = c(-2.5, -1.5, -0.5, 0.5), color = 'darkgrey', linetype = 'dashed') +
#   scale_x_continuous(limits = c(-3.2, 1))+
#   ylab('Gene') +
#   xlab('Gene Effect \n(Chronos score, 23Q2)') +
#   theme_bw()
# 
# ggplot(ge_long_edit_specific, aes(x = value, y = variable, fill = after_stat(x))) +
#   geom_density_ridges_gradient(scale = 2, gradient_lwd = .5, color = "black",
#                                jittered_points = TRUE,
#                                position = position_points_jitter(width = 0.05, height = 0),
#                                point_shape = '|', point_size = 3, point_alpha = 1, alpha = 0.7) +
#   scale_fill_viridis_c(option = "plasma", name = "Gene Effect") +
#   # geom_vline(xintercept = -1, color='darkred')+
#   geom_vline(xintercept = c(1,0,-2), color = 'darkgrey') +
#   geom_vline(xintercept = c(-2.5, -1.5, -0.5, 0.5), color = 'darkgrey', linetype = 'dashed') +
#   scale_x_continuous(limits = c(-3.2, 1))+
#   ylab('Gene') +
#   xlab('Gene Effect \n(Chronos score, 23Q2)') +
#   theme_bw()
# 
# 
# ggplot(ge_long_edit_specific, aes(x = value, y = variable, fill = after_stat(x))) +
#   geom_density_ridges_gradient(scale = 2, gradient_lwd = .5, color = "black",
#                                jittered_points = TRUE,
#                                position = position_points_jitter(width = 0.05, height = 0),
#                                point_shape = '|', point_size = 3, point_alpha = 1, alpha = 0.7) +
#   scale_fill_viridis_c(option = "plasma", name = "Gene Effect") +
#   # geom_vline(xintercept = -1, color='darkred')+
#   geom_vline(xintercept = c(1,0,-2), color = 'darkgrey') +
#   geom_vline(xintercept = c(-2.5, -1.5, -0.5, 0.5), color = 'darkgrey', linetype = 'dashed') +
#   
#   ylab('Gene') +
#   xlab('Gene Effect \n(Chronos score, 23Q2)') +
#   theme_bw()

scale_edge_color_manual(values=c(
  "grey", # LGG
  "#72BC6C", # OV
  "#C0392B", # BRCA
  "#72BC6C", # BLCA
  "#C0392B", # LIHC
  "grey", # STAD
  "#C0392B", # LUSC
  "grey", # CESC
  "#72BC6C", # UCEC
  "#C0392B"  # LUAD
))