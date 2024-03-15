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
num_CancerType =  "04.TCGA-CESC"
total_matrix = data.frame()
for (num_CancerType in Cancerlist) {
  
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  Cancername = gsub('TCGA-' , '', CancerType)
  
  # call input

  gc_TCGA = readRDS(paste0("~/nas/04.Results/short_long/", CancerType,"_critical_features_short_long.rds"))
  gc_cellline_filt_df = readRDS(paste0("~/nas/04.Results/drug/depmap/gdsc/",CancerType, "_DM_sl_cluster.rds"))

  depmap_common_link_genes <- unique(link_genes_filtered_df[which(link_genes_filtered_df$Pathway %in% colnames(gc_cellline_filt_df %>% dplyr::select(-cluster))),]$Genes)
  depmap_common_each_genes <- unique(link_genes_filtered_df[which(single_genes$Pathway %in% colnames(gc_cellline_filt_df %>% dplyr::select(-cluster))),]$Genes)
  depmap_common_genes = unique(c(depmap_common_link_genes,depmap_common_each_genes))

  depmap_score_critical_features = ge.tbl %>% 
    filter(rownames(ge.tbl) %in% rownames(gc_cellline_filt_df)) %>%
    dplyr::select(any_of(depmap_common_genes))
  
  common_depmap = gc_cellline_filt_df[rownames(depmap_score_critical_features),]
  # at_depmap
  if (all.equal(rownames(depmap_score_critical_features), rownames(common_depmap))) {
    depmap_score_critical_features$cluster = common_depmap$cluster
  }
  
  dep_score_cf_filt = depmap_score_critical_features %>% select(-names(which(colSums(is.na(depmap_score_critical_features)) != 0)))
  
  long_dep_score_cf_filt = dep_score_cf_filt %>% filter(cluster == "long")
  short_dep_score_cf_filt = dep_score_cf_filt %>% filter(cluster == "short")

  long_sum = long_dep_score_cf_filt %>% select(-cluster) %>% colMeans() %>% as.data.frame()
  short_sum = short_dep_score_cf_filt %>% select(-cluster) %>% colMeans() %>% as.data.frame()
  
  mid_sum = cbind(long_sum, short_sum)
  colnames(mid_sum) = c(paste0(Cancername , "_long_mean"), paste0(Cancername , "_short_mean"))
  tmp = as.data.frame(t(mid_sum))

  total_matrix = bind_rows(total_matrix, tmp)

}
total_rev = -total_matrix
# total_rev_filt = total_rev + (-min(total_rev , na.rm = T))
total_rev[is.na(total_rev)] = 0

pheatmap(total_rev,
         cluster_rows = F ,
         color = colorRampPalette(c("blue","white" ,"red"))(100))


total_sig_matrix = data.frame()
for (num_CancerType in Cancerlist) {
  
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  Cancername = gsub('TCGA-' , '', CancerType)
  
  # call input
  
  gc_TCGA = readRDS(paste0("~/nas/04.Results/short_long/", CancerType,"_critical_features_short_long.rds"))
  gc_cellline_filt_df = readRDS(paste0("~/nas/04.Results/drug/depmap/gdsc/",CancerType, "_DM_sl_cluster.rds"))
  
  depmap_common_link_genes <- unique(link_genes_filtered_df[which(link_genes_filtered_df$Pathway %in% colnames(gc_cellline_filt_df %>% dplyr::select(-cluster))),]$Genes)
  depmap_common_each_genes <- unique(link_genes_filtered_df[which(single_genes$Pathway %in% colnames(gc_cellline_filt_df %>% dplyr::select(-cluster))),]$Genes)
  depmap_common_genes = unique(c(depmap_common_link_genes,depmap_common_each_genes))
  
  depmap_score_critical_features = ge.tbl %>% 
    filter(rownames(ge.tbl) %in% rownames(gc_cellline_filt_df)) %>%
    dplyr::select(any_of(depmap_common_genes))
  
  common_depmap = gc_cellline_filt_df[rownames(depmap_score_critical_features),]
  # at_depmap
  if (all.equal(rownames(depmap_score_critical_features), rownames(common_depmap))) {
    depmap_score_critical_features$cluster = common_depmap$cluster
  }
  
  dep_score_cf_filt = depmap_score_critical_features %>% select(-names(which(colSums(is.na(depmap_score_critical_features)) != 0)))
  # depmap_gene = "CFLAR"
  
  sig_genes = c()
  for (depmap_gene in dep_score_cf_filt %>% select(-cluster) %>% colnames()) {
    tmp_dep = dep_score_cf_filt %>% select(any_of(depmap_gene) , cluster)
    
    tmp_dep_short = tmp_dep %>% filter(cluster == "short") 
    tmp_dep_long = tmp_dep %>% filter(cluster == "long")

    if (nrow(tmp_dep_short) < 2 | nrow(tmp_dep_long) < 2) {
      next
    }
    anno_tmp = t.test(tmp_dep_short %>% 
                        select(any_of(depmap_gene)) %>% 
                        pull(), tmp_dep_long %>% 
                        select(any_of(depmap_gene)) %>% 
                        pull())$p.value
    
    if (anno_tmp < 0.05) {
      sig_genes = c(sig_genes, depmap_gene)
    } else {
      next
    }
    
  }
  
  if (length(sig_genes) == 0 )  {
    next
  }
  dep_score_sig_cf_filt = dep_score_cf_filt[,c(sig_genes,"cluster")]
  
  long_dep_score_cf_filt = dep_score_sig_cf_filt %>% filter(cluster == "long")
  short_dep_score_cf_filt = dep_score_sig_cf_filt %>% filter(cluster == "short")
  
  long_sum = long_dep_score_cf_filt %>% select(-cluster) %>% colMeans() %>% as.data.frame()
  short_sum = short_dep_score_cf_filt %>% select(-cluster) %>% colMeans() %>% as.data.frame()
  
  mid_sum = cbind(long_sum, short_sum)
  colnames(mid_sum) = c(paste0(Cancername , "_long_mean"), paste0(Cancername , "_short_mean"))
  tmp = as.data.frame(t(mid_sum))
  
  total_sig_matrix = bind_rows(total_sig_matrix, tmp)
  
}

total_sig_matrix

total_sig_rev = -total_sig_matrix
# total_rev_filt = total_rev + (-min(total_rev , na.rm = T))
total_sig_rev[is.na(total_sig_rev)] = 0

# # Convert the matrix to a numeric matrix
# total_sig_rev_back

tmp_numeric <- matrix(as.numeric(unlist(total_sig_rev)), nrow = nrow(total_sig_rev))
vec <- as.vector(tmp_numeric)

hist(vec[vec != 0])
boxplot(vec[vec != 0])

total_sig_rev[total_sig_rev > quantile(vec[vec != 0], probs = 0.75) + 1.5*IQR(vec[vec != 0])] <- 
  round(quantile(vec[vec != 0], probs = 0.75) + 1.5*IQR(vec[vec != 0]),1)

tmp_tt <- matrix(as.numeric(unlist(total_sig_rev)), nrow = nrow(total_sig_rev))
vec_tt <- as.vector(tmp_tt)

hist(vec_tt)
boxplot(vec_tt[vec_tt!=0])

# # Convert the matrix to a numeric matrix
total_sig_rev_back = total_sig_rev

total_sig_rev_back[total_sig_rev_back > 0] = (total_sig_rev_back[total_sig_rev_back>0] - min(total_sig_rev_back[total_sig_rev_back>0])) / 
  (max(total_sig_rev_back[total_sig_rev_back>0]) - min(total_sig_rev_back[total_sig_rev_back>0])) *
  (1 -  min(total_sig_rev_back[total_sig_rev_back>0])) + min(total_sig_rev_back[total_sig_rev_back>0])

total_sig_rev_back[total_sig_rev_back < 0] = (total_sig_rev_back[total_sig_rev_back < 0] - min(total_sig_rev_back[total_sig_rev_back < 0])) / 
  (max(total_sig_rev_back[total_sig_rev_back < 0]) - min(total_sig_rev_back[total_sig_rev_back < 0])) *
  (max(total_sig_rev_back[total_sig_rev_back < 0]) - -1) + -1


annotation_col = data.frame()

for (cancername in unique(sapply(strsplit(rownames(total_sig_rev_back), "_"), "[", 1))) {
  
  nonzero_cols = colnames(total_sig_rev_back)[
    apply(total_sig_rev_back, 2, function(x) all(x[grep(paste0("^",cancername), rownames(total_sig_rev_back))] != 0) & 
            all(x[-grep(paste0("^",cancername), rownames(total_sig_rev_back))] == 0))
  ]
  
  tmp_annocol = data.frame(cancertype = cancername , genes = nonzero_cols)
  annotation_col = rbind(annotation_col, tmp_annocol)
}

other_df = data.frame(cancertype = "multiple_specific", genes = colnames(total_sig_rev_back)[!colnames(total_sig_rev_back) %in% annotation_col$genes])

multi_genes = other_df$genes
multi_only = total_sig_rev_back %>% 
  select(any_of(multi_genes))
# mg = "CSF1"
other_sort_df = data.frame()
for (mg in colnames(multi_only)) {
  tmp_m = multi_only %>% select(any_of(mg)) %>% filter(. !=0)
  tmp_other = data.frame(cancertype = unlist(strsplit(rownames(tmp_m)[which(tmp_m == max(tmp_m))], "_"))[1] , genes = mg)
  other_sort_df = rbind(other_sort_df,tmp_other) 
}

annotation_col = rbind(annotation_col, other_sort_df)
rownames(annotation_col) <- annotation_col$genes
annotation_col = annotation_col %>% select(-genes)
annotation_col = annotation_col[colnames(total_sig_rev_back), , drop = FALSE]

annotation_col = annotation_col %>% mutate(spe_cluster = case_when(rownames(.) %in% multi_genes ~ "multi_specific",
                                                  .default = "cancer_specific"))

annotation_row = data.frame(types = rownames(total_sig_rev_back))
annotation_row  = annotation_row %>% 
  mutate(cluster = sapply(strsplit(annotation_row$types, "_"), "[", 1)) %>%
  mutate(sl = sapply(strsplit(annotation_row$types, "_"), "[", 2))
rownames(annotation_row) = annotation_row$types

annotation_row = annotation_row %>% select(-types) %>% select(sl, cluster)

# c(
#   "LGG" = "#00A087FF", # LGG
#   "OV" = "#3C5488FF", # OV
#   "BRCA" = "#4DBBD5FF", # BRCA
#   "BLCA" = "#7E6148FF", # BLCA
#   "LIHC" = "#8491B4FF", # LIHC
#   "STAD" = "#91D1C2FF", # STAD
#   "LUSC" = "#B09C85FF", # LUSC
#   "CESC" = "#DC0000FF", # CESC
#   "UCEC" = "#E64B35FF", # UCEC
#   "LUAD" = "#F39B7FFF"  # LUAD
# )
cancer_colors = c(
  "LGG" = "#00A087FF", # LGG
  "OV" = "#3C5488FF", # OV
  "BRCA" = "#4DBBD5FF", # BRCA
  "BLCA" = "#7E6148FF", # BLCA
  "LIHC" = "#8491B4FF", # LIHC
  "STAD" = "#91D1C2FF", # STAD
  "LUSC" = "#B09C85FF", # LUSC
  "CESC" = "#DC0000FF", # CESC
  "UCEC" = "#E64B35FF", # UCEC
  "LUAD" = "#F39B7FFF" # LUAD
)

spe_colors = c(
  "cancer_specific" = "#e62a89", # cancerspe
  "multi_specific" = "#189d77" # muti
)
row_colors_cluster <- c(
  "LGG" = "#00A087FF", # LGG
  "OV" = "#3C5488FF", # OV
  "BRCA" = "#4DBBD5FF", # BRCA
  "BLCA" = "#7E6148FF", # BLCA
  "LIHC" = "#8491B4FF", # LIHC
  "STAD" = "#91D1C2FF", # STAD
  "LUSC" = "#B09C85FF", # LUSC
  "CESC" = "#DC0000FF", # CESC
  "UCEC" = "#E64B35FF", # UCEC
  "LUAD" = "#F39B7FFF"  # LUAD
)
row_colors_sl <- c("short" = "red", "long" = "#009E73")

ann_colors_sl = list(spe_cluster = spe_colors, cancertype = cancer_colors ,cluster = row_colors_cluster ,  sl =  row_colors_sl)

ComplexHeatmap::pheatmap(total_sig_rev_back %>% as.matrix(),
                         cluster_rows = F,
                         cluster_cols = F,
                         # breaks = c(-Inf,max(total_sig_rev[total_sig_rev < 0]) , 0 , min(total_sig_rev[total_sig_rev > 0]), Inf) ,
                         column_split = factor(annotation_col$cancertype, levels = c(unique(annotation_col$cancertype))) ,
                         row_split = factor(annotation_row$cluster, levels = c(unique(annotation_row$cluster))),
                         annotation_colors = ann_colors_sl,
                         annotation_col = annotation_col %>% select(spe_cluster),
                         annotation_row = annotation_row,
                         legend = T,
                         annotation_legend = T,
                         annotation_names_col = F,
                         annotation_names_row = F,
                         show_colnames = T,
                         show_rownames = F,
                         cluster_column_slices = FALSE,
                         color = c( colorRampPalette(c("#387eb8", "#F9FEFE", "#e21e26"))(100)),
                         scale = "none",
                         border_color = NA
                         )

dev.off()
#FCFDFE
#F9FEFE
#FDFDFE