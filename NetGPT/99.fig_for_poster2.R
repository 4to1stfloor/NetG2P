


library(tidyverse)
library(readxl)
library(org.Hs.eg.db)
library(dplyr)
library(clusterProfiler)

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

surv_total_results = read_xlsx("~/nas/04.Results/Total_results_survpval2.xlsx")
num_CancerType = "19.TCGA-LIHC"

main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))

# call input
cancer_slc = read.csv(paste0("~/nas/04.Results/short_long/ttest_common/",CancerType,"_critical_features_short_long_common.csv"))

short_features = cancer_slc[which(cancer_slc$classification == "short"),]$variable 
long_features = cancer_slc[which(cancer_slc$classification == "long"),]$variable 
common_features = cancer_slc[which(cancer_slc$classification == "common"),]$variable 

short_gene = c()
long_gene = c()
common_gene = c()

# short
for (sf in short_features) {
  count <- str_count(sf, "P")
  if (count == 1) {
    tmp_short_gene = single_genes[which(single_genes$Pathway == sf),]$Genes
  } else {
    tmp_short_gene = link_genes_filtered_df[which(link_genes_filtered_df$Pathway == sf),]$Genes
  }
  short_gene = c(short_gene, tmp_short_gene)
  
}

# long
for (lf in long_features) {
  count <- str_count(lf, "P")
  if (count == 1) {
    tmp_long_gene = single_genes[which(single_genes$Pathway == lf),]$Genes
  } else {
    tmp_long_gene = link_genes_filtered_df[which(link_genes_filtered_df$Pathway == lf),]$Genes
  }
  long_gene = c(long_gene, tmp_long_gene)
  
}

# long
for (cmf in common_features) {
  count <- str_count(lf, "P")
  if (count == 1) {
    tmp_common_gene = single_genes[which(single_genes$Pathway == cmf),]$Genes
  } else {
    tmp_common_gene = link_genes_filtered_df[which(link_genes_filtered_df$Pathway == cmf),]$Genes
  }
  common_gene = c(common_gene, tmp_common_gene)
  
}


short_gene_en = AnnotationDbi::select(org.Hs.eg.db, short_gene, 'ENTREZID', 'SYMBOL')[
  which(!is.na(AnnotationDbi::select(org.Hs.eg.db, short_gene, 'ENTREZID', 'SYMBOL')$ENTREZID)),]
short_gene_en <- data.frame(short_gene_en, row.names = NULL)

long_gene_en = AnnotationDbi::select(org.Hs.eg.db, long_gene, 'ENTREZID', 'SYMBOL')[
  which(!is.na(AnnotationDbi::select(org.Hs.eg.db, long_gene, 'ENTREZID', 'SYMBOL')$ENTREZID)),]
long_gene_en <- data.frame(long_gene_en, row.names = NULL)

common_gene_en = AnnotationDbi::select(org.Hs.eg.db, common_gene, 'ENTREZID', 'SYMBOL')[
  which(!is.na(AnnotationDbi::select(org.Hs.eg.db, common_gene, 'ENTREZID', 'SYMBOL')$ENTREZID)),]
common_gene_en <- data.frame(common_gene_en, row.names = NULL)


short_enrichGO = enrichGO(short_gene_en$ENTREZID ,
                          OrgDb = org.Hs.eg.db ,
                          keyType = "ENTREZID" , 
                          ont = "BP" ,
                          pvalueCutoff  = 0.05,
                          qvalueCutoff = 0.01)

short_enrichGO <- setReadable(short_enrichGO, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

dotplot(short_enrichGO) 

long_enrichGO = enrichGO(long_gene_en$ENTREZID ,
                         OrgDb = org.Hs.eg.db ,
                         keyType = "ENTREZID" , 
                         ont = "BP" ,
                         pvalueCutoff  = 0.05,
                         qvalueCutoff = 0.01 )

long_enrichGO <- setReadable(long_enrichGO, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

dotplot(long_enrichGO) 

common_enrichGO = enrichGO(common_gene_en$ENTREZID ,
                         OrgDb = org.Hs.eg.db ,
                         keyType = "ENTREZID" , 
                         ont = "BP" ,
                         pvalueCutoff  = 0.05,
                         qvalueCutoff = 0.01 )

common_enrichGO <- setReadable(common_enrichGO, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

dotplot(common_enrichGO) 


top_genes_group = list(short_cluster = short_gene_en$ENTREZID, long_cluster = long_gene_en$ENTREZID)

ck <- compareCluster(geneCluster = top_genes_group, fun = "enrichKEGG")
ck <- setReadable(ck, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

png(filename = paste0(CancerType,"_dotplot.png"),
    width = 25, height = 25,  units = "cm" ,pointsize = 12,
    bg = "white", res = 1200, family = "")

exp_dot = dotplot(ck) + ggtitle("Exp_for_Top300")

print(exp_dot)
dev.off()
