#making heatmaps
#Short/long analysis

library(tidyverse)
library(readxl)
library(org.Hs.eg.db)
library(dplyr)
library(clusterProfiler)
library(tidyverse)

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

ck <- compareCluster(geneCluster = top_genes_group, 
                     fun = "enrichGO",
                     OrgDb='org.Hs.eg.db',
                     pvalueCutoff  = 0.05,
                     qvalueCutoff = 0.0001 )

ck <- setReadable(ck, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

# png(filename = paste0(CancerType,"_dotplot.png"),
#     width = 25, height = 25,  units = "cm" ,pointsize = 12,
#     bg = "white", res = 1200, family = "")

total_dot = dotplot(ck)
setwd("~/nas/04.Results/GOenrichment_test/BP/")
ggsave(file="total_cancer_total_cf.svg", plot=total_dot, width=20, height=10)

# print(exp_dot)
# dev.off()


### divide features

total_divide_features = list()

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
  
  total_divide_features[[paste0(Cancername, "_long")]] = cancer_bf_sl[which(cancer_bf_sl$classification == "long"),]$variable 
  total_divide_features[[paste0(Cancername, "_short")]] = cancer_bf_sl[which(cancer_bf_sl$classification == "short"),]$variable 
  total_divide_features[[paste0(Cancername, "_common")]] = cancer_bf_sl[which(cancer_bf_sl$classification == "common"),]$variable 
  
}

total_divide_genes_list = list()
for (cancer_features in names(total_divide_features)) {
  
  tmp_gene = c()
  
  sg = single_genes[which(single_genes$Pathway %in% total_divide_features[[cancer_features]]),]$Genes
  lg = link_genes_filtered_df[which(link_genes_filtered_df$Pathway %in% total_divide_features[[cancer_features]]),]$Genes
  tmp_gene = c(sg, lg)
  tmp_gene = unique(tmp_gene)
  
  total_divide_genes_list[[cancer_features]] = tmp_gene
  
  
}
library(org.Hs.eg.db)
library(DOSE)

divide_genes_group = list()

for (cancer_genes in names(total_divide_features)) {
  
  cancer_divide_gene_en = AnnotationDbi::select(org.Hs.eg.db, total_divide_genes_list[[cancer_genes]], 'ENTREZID', 'SYMBOL')[
    which(!is.na(AnnotationDbi::select(org.Hs.eg.db, total_divide_genes_list[[cancer_genes]], 'ENTREZID', 'SYMBOL')$ENTREZID)),]
  
  cancer_divide_gene_en <- data.frame(cancer_divide_gene_en, row.names = NULL)
  
  divide_genes_group[[paste0(cancer_genes, "_cluster")]] = cancer_divide_gene_en$ENTREZID
  
}

ck_divide <- compareCluster(geneCluster = divide_genes_group, 
                     fun = "enrichGO",
                     OrgDb='org.Hs.eg.db',
                     ont = "BP",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff = 0.0001 )

ck_divide <- setReadable(ck_divide, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

# png(filename = paste0(CancerType,"_dotplot.png"),
#     width = 25, height = 25,  units = "cm" ,pointsize = 12,
#     bg = "white", res = 1200, family = "")

total_divide_dot = dotplot(ck_divide)
ggsave(file="total_cancer_total_cf_divide.svg", plot=total_divide_dot, width=20, height=10)

?compareCluster
scl_res = ck_divide@compareClusterResult

test_df = data.frame()
aaa = "CESC_long_cluster"
for (aaa in names(divide_genes_group)) {
  if (length(divide_genes_group[[aaa]]) != 0) {
    tmp_df = data.frame(str_split(aaa, "_")[[1]][1], str_split(aaa, "_")[[1]][2], divide_genes_group[[aaa]])
    colnames(tmp_df) = c("cancertype" , "classification", "ENTREZID")
    test_df = rbind(test_df, tmp_df)
  } else {
    next
  }
  
}

group_res = compareCluster(
  ENTREZID~classification+cancertype,
  data = test_df, 
  fun = "enrichGO",
  OrgDb='org.Hs.eg.db',
  ont = "BP",
  pvalueCutoff  = 0.05,
  qvalueCutoff = 0.0001 )

ttt = dotplot(group_res)

ggsave(file="test_ttt.svg", plot=ttt, width=50, height=60,limitsize = FALSE)

test_df = group_res@compareClusterResult

metabolism = c("metabolic ", " ion ","cadmium ion", "glycolytic ", "tricarboxylic acid",
               "glucose ", "catabolic ", "ATP ", 
               "proton transmembrane ", "acetyl-CoA", "electron transport","biosynthetic process","hexose","pentose-phosphate",
               "NADH regeneration","urea cycle","hydroxylation","NADH regeneration")

immune = c("cytokine", "leukocyte", "immune", "T cell", "lymphocyte", "Fc receptor", "cell activation","mononuclear cell",
           " lipopolysaccharide", " defense response", " antigen", 
           "phagocytosis", "response to virus","hemopoiesis", "chemotaxis", "cell killing", "autophagy", "MHC","repair","immunity", "response to biotic stimulus",
           "interferon","viral","inflammatory","bacterial origin","repair","complement activation")

signaling = c("signaling", "kinase", "cascade", "autophosphorylation","phosphorylation", "dephosphorylation","cellular response", "protein binding",
              "secretion","signal","signal transduction", "transcription factor","GTPase","peptide transport",
              "SMAD protein", "protein transport", "intracellular transport","zymogen","chemical synaptic transmission","immunoglobulin")

cell_cycle_develop = c("cell cycle","nuclear division", "cell division", "proliferation", "cell growth", "gland development", "gilogenesis", "cell development"
                       ,"renal system development", "protein localization", "angiogenesis","DNA replication", "differentiation", "epithelial tube",
                       "cell projection", "cellular component","pericardium dyevelopment",
                       "anterior/posterior pattern specification","ureteric bud development"," development","developmental growth","cell fate",
                       "telomere")

migration = c("migration", "actin", "lamellipodium", "epithelial to mesenchymal", "fiber assembly","mesenchyme development")  

apoptosis = c("cell death", "apoptosis","neuron death")


test_filted_df = test_df %>% mutate(total_imp = case_when(Description %in% 
                                                            test_df$Description[grep(paste(metabolism, collapse = "|"), test_df$Description)] ~ "Metabolism",
                                                          Description %in% 
                                                            test_df$Description[grep(paste(immune, collapse = "|"), test_df$Description)] ~ "Immune",
                                                          Description %in%
                                                            test_df$Description[grep(paste(signaling, collapse = "|"), test_df$Description)] ~ "Signaling",
                                                          Description %in%
                                                            test_df$Description[grep(paste(cell_cycle_develop, collapse = "|"), test_df$Description)] ~ "Cell cycle",
                                                          Description %in%
                                                            test_df$Description[grep(paste(migration, collapse = "|"), test_df$Description)] ~ "Migration",
                                                          Description %in%
                                                            test_df$Description[grep(paste(apoptosis, collapse = "|"), test_df$Description)] ~ "Apoptosis",
                                                          .default = "")) 
  
color_category = c("#7E6148FF","#8491B4FF","#00A087FF","#4DBBD5FF","#E64B35FF","#3C5488FF")


test_res = group_res
test_res@compareClusterResult = test_filted_df
clusterProfiler::dotplot(test_res, split="cancertype", showCategory = 5)
class(test_res)

head(group_res@compareClusterResult)



test_go = enrichGO(cancer_divide_gene_en$ENTREZID,
         OrgDb = org.Hs.eg.db ,
         keyType = "ENTREZID" , 
         ont = "BP" ,
         pvalueCutoff  = 0.05,
         qvalueCutoff = 0.0001 )

test_go_df = test_go@result


total_slc_features = list()

# for all
tmp_short = c()
tmp_long = c()
tmp_common = c()

for (num_CancerType in Cancerlist) {
  
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  Cancernum = gsub('[.-]','',gsub('[a-zA-Z]','', num_CancerType))
  # call input
  
  # cancer_bf = read.csv(paste0(filepath,"04.Results/bestfeatures/",CancerType, "_critical_features.csv"))
  cancer_bf_sl = read.csv(paste0(filepath,"04.Results/short_long/ttest_common/",CancerType, "_critical_features_short_long_common.csv"))
  
  # cancer_bf_cut = cancer_bf[1:surv_total_results[which(surv_total_results$CancerType == CancerType),]$num_of_features,]
  Cancername = gsub('TCGA-','', CancerType)
  
  tmp_long = c(tmp_short, cancer_bf_sl[which(cancer_bf_sl$classification == "long"),]$variable )
  tmp_short = c(tmp_short, cancer_bf_sl[which(cancer_bf_sl$classification == "short"),]$variable )
  tmp_common = c(tmp_short, cancer_bf_sl[which(cancer_bf_sl$classification == "common"),]$variable )
  
}

total_slc_features[["long"]] = unique(tmp_long)
total_slc_features[["short"]] = unique(tmp_short)
total_slc_features[["common"]] = unique(tmp_common)

total_slc_genes_list = list()
for (slc in names(total_slc_features)) {
  
  tmp_gene = c()
  
  sg = single_genes[which(single_genes$Pathway %in% total_slc_features[[slc]]),]$Genes
  lg = link_genes_filtered_df[which(link_genes_filtered_df$Pathway %in% total_slc_features[[slc]]),]$Genes
  tmp_gene = c(sg, lg)
  tmp_gene = unique(tmp_gene)
  
  total_slc_genes_list[[slc]] = tmp_gene
  
}

total_slc_genes_list

library(org.Hs.eg.db)
library(DOSE)

slc_genes_group = list()

for (slc_genes in names(total_slc_features)) {
  
  cancer_slc_gene_en = AnnotationDbi::select(org.Hs.eg.db, total_slc_genes_list[[slc_genes]], 'ENTREZID', 'SYMBOL')[
    which(!is.na(AnnotationDbi::select(org.Hs.eg.db, total_slc_genes_list[[slc_genes]], 'ENTREZID', 'SYMBOL')$ENTREZID)),]
  
  cancer_slc_gene_en <- data.frame(cancer_slc_gene_en, row.names = NULL)
  
  slc_genes_group[[paste0(slc_genes, "_cluster")]] = cancer_slc_gene_en$ENTREZID
  
}

ck_slc <- compareCluster(geneCluster = slc_genes_group, 
                            fun = "enrichGO",
                            OrgDb='org.Hs.eg.db',
                            pvalueCutoff  = 0.05,
                            qvalueCutoff = 0.0001 )

ck_slc <- setReadable(ck_slc, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

# png(filename = paste0(CancerType,"_dotplot.png"),
#     width = 25, height = 25,  units = "cm" ,pointsize = 12,
#     bg = "white", res = 1200, family = "")

total_slc_dot = dotplot(ck_slc)
ggsave(file="total_cancer_total_cf_slc.svg", plot=total_slc_dot, width=10, height=10)

