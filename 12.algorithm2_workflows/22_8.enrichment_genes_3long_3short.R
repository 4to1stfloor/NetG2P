#making heatmaps
#Short/long analysis

library(openxlsx)
crit.sl.features = read.xlsx("~/nas/04.Results/critical_features/most_common_with_short_long_critical_features.xlsx")

remove.cancer = c("TCGA-COADREAD", "TCGA-LUSC")
crit.sl.features = crit.sl.features[,!colnames(crit.sl.features) %in% remove.cancer]


#making short long common balance
crit.sl.features$short.rat = NA
crit.sl.features$long.rat = NA
crit.sl.features$common.rat = NA

for (j in 1:nrow(crit.sl.features)) {
  query.cancer.type = crit.sl.features$which_cancer[j]
  query.cancer = sort(unlist(strsplit(query.cancer.type, ", ")))
  query.cancer.clean= query.cancer[!query.cancer %in% remove.cancer]
  
  crit.sl.features$count[j] = length(query.cancer.clean)
  crit.sl.features$which_cancer[j] = paste(query.cancer.clean, collapse = ", ")
  
  crit.sl.features$short.rat[j] = length(grep("short", crit.sl.features[j,])) / crit.sl.features$count[j]
  crit.sl.features$long.rat[j] = length(grep("long", crit.sl.features[j,])) / crit.sl.features$count[j]
  crit.sl.features$common.rat[j] = length(grep("common", crit.sl.features[j,])) / crit.sl.features$count[j]
}

crit.sl.features = crit.sl.features[order(crit.sl.features$count, decreasing = T),]
common.features = crit.sl.features[crit.sl.features$count > 1,]

common.features$`TCGA-LGG`


#make a heatmap
################################################################################
common.feat.mat = common.features[,grep("TCGA", colnames(common.features))]
rownames(common.feat.mat) = common.features$features

common.feat.mat[common.feat.mat == "long"] = -1
common.feat.mat[common.feat.mat == "common"] = 0
common.feat.mat[common.feat.mat == "short"] = 1

common.feat.mat = as.matrix(common.feat.mat)
common.n.mat = apply(common.feat.mat, 2, as.numeric)
rownames(common.n.mat) = rownames(common.feat.mat)

common.n.raw = common.n.mat
#masking NAs
#dont do it too often
common.n.mat[is.na(common.n.mat)] = 0
#removing with more than 8 zeros
zero.test = rowSums(common.n.mat != 0)

common.nn.mat = common.n.mat[zero.test > 2,]

library(pheatmap)
pheatmap(common.n.mat, cluster_rows = T, cluster_cols = T)
pheatmap(common.nn.mat, cluster_rows = T, cluster_cols = T,clustering_method = "mcquitty")


#
loss.count = rowSums(common.nn.mat[,colnames(common.nn.mat) %in% c("TCGA-LIHC", "TCGA-KIDNEY")])
loss.feat = names(loss.count)[abs(loss.count) > 1]

gain.count = rowSums(common.nn.mat[,colnames(common.nn.mat) %in% c("TCGA-STAD", "TCGA-OV", "TCGA-BLCA", "TCGA-UCEC")])
gain.feat = names(gain.count)[abs(gain.count) > 1]

common.filt = common.nn.mat[rownames(common.nn.mat) %in% unique(c(loss.feat, gain.feat)),]
test = pheatmap(common.filt, cluster_rows = T, cluster_cols = T)


cancer.order = c("TCGA-LIHC", "TCGA-KIDNEY", "TCGA-BRCA", "TCGA-LUAD", "TCGA-LGG", "TCGA-CESC", "TCGA-STAD", "TCGA-OV", "TCGA-BLCA", "TCGA-UCEC")
#extracting gene order
gene.order = test$tree_row$labels[test$tree_row$order]


common.n.mat = apply(common.feat.mat, 2, as.numeric)
rownames(common.n.mat) = rownames(common.feat.mat)

common.custom.mat = common.n.mat[match(gene.order, rownames(common.n.mat)), match(cancer.order, colnames(common.n.mat))]

common.cus.pheat = pheatmap(common.custom.mat, cluster_rows = F, cluster_cols = F)
################################################################################




#making 3long graph
feat.sel = rowSums(common.n.mat[, colnames(common.n.mat) %in% c("TCGA-OV", "TCGA-BLCA", "TCGA-UCEC")])
# feat.sel = rowSums(common.n.mat[, colnames(common.n.mat) %in% c("TCGA-OV", "TCGA-BLCA", "TCGA-UCEC")])

feat.3l = names(feat.sel)[which(feat.sel == -3)]

common_mat = common.n.mat[, colnames(common.n.mat) %in% c("TCGA-LIHC", "TCGA-KIDNEY", "TCGA-LUAD")]
common_mat[common_mat == 0] = 5
common_mat[is.na(common_mat)] = 0

feat.see = rowSums(common_mat)
feat.see_filt = feat.see[!feat.see %in% c(-1,1,0)]

feat.see_filt[feat.see_filt == 5] = 0
feat.see_filt[feat.see_filt == 10] = 0
feat.see_filt[feat.see_filt == 4] = -1


feat.3s = names(feat.see_filt)[which(feat.see_filt == 2)]
feat.3l

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


main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))

# call input

# short - > KIDNEY , LUAD, LIHC
# long - > OV, BLCA, UCEC

short_features = feat.3s
long_features = feat.3l

short_gene = c()
long_gene = c()

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

short_gene_en = AnnotationDbi::select(org.Hs.eg.db, short_gene, 'ENTREZID', 'SYMBOL')[
  which(!is.na(AnnotationDbi::select(org.Hs.eg.db, short_gene, 'ENTREZID', 'SYMBOL')$ENTREZID)),]
short_gene_en <- data.frame(short_gene_en, row.names = NULL)

long_gene_en = AnnotationDbi::select(org.Hs.eg.db, long_gene, 'ENTREZID', 'SYMBOL')[
  which(!is.na(AnnotationDbi::select(org.Hs.eg.db, long_gene, 'ENTREZID', 'SYMBOL')$ENTREZID)),]
long_gene_en <- data.frame(long_gene_en, row.names = NULL)

short_enrichGO = enrichGO(short_gene_en$ENTREZID ,
                          OrgDb = org.Hs.eg.db ,
                          keyType = "ENTREZID" , 
                          ont = "BP" ,
                          pvalueCutoff  = 0.05,
                          qvalueCutoff = 0.01)

short_enrichGO <- setReadable(short_enrichGO, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

dotplot(short_enrichGO) 

long_enrichGO = enrichGO(long_gene_en$ENTREZID,
                         OrgDb = org.Hs.eg.db ,
                         keyType = "ENTREZID" , 
                         ont = "BP" ,
                         pvalueCutoff  = 0.05,
                         qvalueCutoff = 0.0001 )

long_enrichGO <- setReadable(long_enrichGO, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

dotplot(long_enrichGO) + dotplot(short_enrichGO) 

library(DOSE)
data(geneList)

gseDO(as.numeric(long_gene_en$ENTREZID))

top_genes_group = list(short_cluster = short_gene_en$ENTREZID, long_cluster = long_gene_en$ENTREZID)
?compareCluster
ck <- compareCluster(geneCluster = top_genes_group, 
                     fun = "enrichGO",
                     OrgDb='org.Hs.eg.db')

ck <- setReadable(ck, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

png(filename = paste0(CancerType,"_dotplot.png"),
    width = 25, height = 25,  units = "cm" ,pointsize = 12,
    bg = "white", res = 1200, family = "")

exp_dot = dotplot(ck)

cnetplot(ck)

print(exp_dot)
dev.off()
