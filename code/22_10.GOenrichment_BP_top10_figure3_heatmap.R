#making heatmaps
#Short/long analysis

library(tidyverse)
library(readxl)
library(org.Hs.eg.db)
library(dplyr)
library(clusterProfiler)
library(tidyverse)
library(ComplexHeatmap)
library(svglite)

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

Cancerlist = dir(paste0(filepath, "/00.data/filtered_TCGA/"))

Cancerlist = Cancerlist[c(-11,-12)]
total_features = list()

# for all
for (num_CancerType in Cancerlist) {
  
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  Cancernum = gsub('[.-]','',gsub('[a-zA-Z]','', num_CancerType))
  # call input
  
  # cancer_bf = read.csv(paste0(filepath,"04.Results/bestfeatures/",CancerType, "_critical_features.csv"))
  cancer_bf_sl = read.csv(paste0(filepath,"04.Results/short_long/ttest_common/",CancerType, "_critical_features_short_long_common.csv"))
  
  # cancer_bf_cut = cancer_bf[1:surv_total_results[which(surv_total_results$CancerType == CancerType),]$num_of_features,]
  Cancername = gsub('TCGA-','', CancerType)
  
  total_features[[Cancername]] = cancer_bf_sl$variable 
  
}

# call input

total_genes_list = list()
for (Cancername in names(total_features)) {
  
  tmp_gene = c()
  
  sg = single_genes[which(single_genes$Pathway %in% total_features[[Cancername]]),]$Genes
  lg = link_genes_filtered_df[which(link_genes_filtered_df$Pathway %in% total_features[[Cancername]]),]$Genes
  tmp_gene = c(sg, lg)
  tmp_gene = unique(tmp_gene)
  
  total_genes_list[[Cancername]] = tmp_gene
  
  
}
library(org.Hs.eg.db)
library(DOSE)

top_genes_group = list()

for (Cancername in names(total_features)) {
  
  cancer_gene_en = AnnotationDbi::select(org.Hs.eg.db, total_genes_list[[Cancername]], 'ENTREZID', 'SYMBOL')[
    which(!is.na(AnnotationDbi::select(org.Hs.eg.db, total_genes_list[[Cancername]], 'ENTREZID', 'SYMBOL')$ENTREZID)),]
  
  cancer_gene_en <- data.frame(cancer_gene_en, row.names = NULL)
  
  top_genes_group[[paste0(Cancername, "_cluster")]] = cancer_gene_en$ENTREZID
  
}

# short_enrichGO = enrichGO(short_gene_en$ENTREZID ,
#                           OrgDb = org.Hs.eg.db ,
#                           keyType = "ENTREZID" , 
#                           ont = "BP" ,
#                           pvalueCutoff  = 0.05,
#                           qvalueCutoff = 0.01)
# 
# short_enrichGO <- setReadable(short_enrichGO, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
?compareCluster
ck <- compareCluster(geneCluster = top_genes_group, 
                     fun = "enrichGO",
                     ont = "BP",
                     OrgDb='org.Hs.eg.db',
                     pvalueCutoff  = 0.05,
                     qvalueCutoff = 0.0001 )

ck <- setReadable(ck, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

total_enrich = ck@compareClusterResult

total_enrich_df = total_enrich%>%
  mutate(GeneRatio_numeric = sapply(str_split(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))) 

min_max_scale <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

crop_enrich = total_enrich_df %>%
  group_by(Cluster) %>%
  arrange(desc(GeneRatio_numeric)) %>%
  slice_head(n = 10)

crop_enrich = total_enrich_df %>% filter(Description %in% unique(crop_enrich$Description))
# crop_enrich = test %>% filter(Description %in% unique(crop_enrich$Description))


crop_enrich$log_qvalue = -log10(crop_enrich$qvalue)
crop_enrich = crop_enrich %>% 
  group_by(Cluster) %>%
  mutate(nor_qvalue = min_max_scale(log_qvalue))
crop_enrich_matrix = reshape2::dcast(crop_enrich, Cluster ~ Description, value.var = "nor_qvalue", fill = 0)
rownames(crop_enrich_matrix) <- gsub("_cluster","",crop_enrich_matrix$Cluster)
crop_enrich_matrix <- crop_enrich_matrix[,-1]

managed_colname = c()
# cln = "oxidoreduction-driven active transmembrane transporter activity"
for (cln in colnames(crop_enrich_matrix)) {
  tmp_cln = unlist(str_split(cln , " "))
  if (length(tmp_cln) %% 2 == 0 & length(tmp_cln) >2) {
    tmp_edit = paste0(paste(tmp_cln[1:(length(tmp_cln) / 2)] , collapse = " "), "\n" , paste(tmp_cln[(length(tmp_cln) / 2 + 1): length(tmp_cln)], collapse = " " ))
    managed_colname = c(managed_colname,tmp_edit)
    
  } else if (length(tmp_cln) %% 2 != 0 & length(tmp_cln) >3) {
    tmp_edit = paste0(paste(tmp_cln[1:ceiling(length(tmp_cln) / 2)] , collapse = " "), "\n" , paste(tmp_cln[(ceiling(length(tmp_cln) / 2) + 1): length(tmp_cln)], collapse = " " ))
    managed_colname = c(managed_colname,tmp_edit)
    
  } else {
    managed_colname = c(managed_colname,cln)
  }
}


setwd("~/nas/04.Results/GOenrichment_test/BP/")

svglite(filename = "figure3_GOenrichment_BP_top10.svg", 
        bg = "white", 
        pointsize = 12,
        width = 18,
        height = 10)
go_heat = pheatmap(crop_enrich_matrix, 
         cluster_rows = T,
         cluster_cols = T,
         show_row_dend = T,
         fontfamily = "arial",
         fontface = "bold",
         labels_col = str_to_title(managed_colname),
         col = circlize::colorRamp2(c( 0, 1), c( "white", "red")),
         row_names_side = "left", 
         cellheight = 10,
         # cellwidth = 10,
         # col_names_max_width = unit(12, "cm"),
         clustering_method = "complete")
print(go_heat)
dev.off()

# "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid"
