library(dplyr)
library(tidyr)
library(ComplexHeatmap)
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

Cancerlist = dir(paste0(filepath, "/00.data/filtered_TCGA/"))
cancerhallmark_geneset = read.csv("~/nas/99.reference/Zhang2020_CHG_filt.csv")

convert_to_term_gene <- function(df) {
  term_gene_list <- do.call(rbind, lapply(names(df), function(term) {
    data.frame(term = term, gene = df[[term]], stringsAsFactors = FALSE)
  }))
  return(term_gene_list)
}

term_gene_df <- convert_to_term_gene(cancerhallmark_geneset)
term_gene_filtered_df = term_gene_df %>% filter(gene != "")

####
Cancerlist = Cancerlist[c(-11,-12)]
# num_CancerType = "04.TCGA-CESC"

total_counts = sum(sapply(cancerhallmark_geneset, function(x) sum(x != "")))
total_hallmark_hyper = data.frame()
total_enrich_df = data.frame()
for (num_CancerType in Cancerlist) {
  
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  Cancername = gsub('TCGA-' , '', CancerType)
  
  critical_features = read.csv(paste0("~/nas/04.Results/bestfeatures/",CancerType,"_critical_features.csv"))
  single_critical_genes = single_genes %>% 
    filter(Pathway %in% critical_features$variable) %>%
    pull(Genes) %>%
    unique()
  
  link_critical_genes = link_genes_filtered_df %>% 
    filter(Pathway %in% critical_features$variable) %>%
    pull(Genes) %>%
    unique()
  
  total_critical_genes = unique(c(single_critical_genes , link_critical_genes))
  # hallmark = "Activating_Invasion_and_Metastasis"
  
  tmp = enricher(total_critical_genes, 
                 pvalueCutoff = 0.05,
                 TERM2GENE = term_gene_filtered_df, 
                 minGSSize=1, 
                 maxGSSize = 1500)
  
  tmp_df = as.data.frame(tmp@result)
  tmp_df_select = tmp_df %>% select(ID,pvalue)
  rownames(tmp_df_select) = NULL
  tmp_df_select = tmp_df_select %>% 
    mutate(cancername = Cancername)
  
  total_enrich_df = rbind(total_enrich_df, tmp_df_select)
  # 
  # rownames(hallmark_hyper) = hallmark_hyper$hallmark
  # hallmark_hyper_filt = hallmark_hyper %>% select(pvalue)
  # colnames(hallmark_hyper_filt) = Cancername
  # hallmark_hyper_filt[hallmark_hyper_filt > 0.05] = 1
  # hallmark_hyper_filt = -log(hallmark_hyper_filt)
  # 
  # if (sum(hallmark_hyper_filt == 0) != 0) {
  #   hallmark_hyper_filt = (hallmark_hyper_filt - min(hallmark_hyper_filt)) / (max(hallmark_hyper_filt) - min(hallmark_hyper_filt))
  # } else {
  #   hallmark_hyper_filt = (hallmark_hyper_filt - 0.1) / (max(hallmark_hyper_filt) - 0.1)
  # }
  # 
  # if (Cancername == "CESC") {
  #   total_hallmark_hyper = hallmark_hyper_filt
  # } else {
  #   total_hallmark_hyper = cbind(total_hallmark_hyper , hallmark_hyper_filt)
  # }
  
}
library(reshape)
total_enrich_df_reshape = total_enrich_df %>% 
  pivot_wider(names_from = cancername , values_from = pvalue) %>%
  as.data.frame()
rownames(total_enrich_df_reshape) = total_enrich_df_reshape$ID
total_enrich_df_reshape$ID = NULL

total_enrich_df_reshape[total_enrich_df_reshape > 0.05] = 1
total_enrich_df_reshape_log = -log10(total_enrich_df_reshape)

# 
# total_hallmark_hyper_filt = total_hallmark_hyper[names(sort(rowSums(total_hallmark_hyper), decreasing = T)),]
# total_hallmark_hyper_filt = total_hallmark_hyper_filt %>% select(LGG,BRCA,STAD,LUAD,LIHC,LUSC,OV,BLCA,UCEC,CESC)
# 
# library(viridis)
# library(svglite)
# 
# total_hallmark_hyper_filt[total_hallmark_hyper_filt > 0.5] = 1
# total_hallmark_hyper_filt[total_hallmark_hyper_filt <= 0.5] = 0

total_enrich_df_reshape_log_bi = total_enrich_df_reshape_log
total_distri = total_enrich_df_reshape_log_bi %>% 
  select(OV, UCEC,BLCA)
total_distri[total_distri <= 29.7] = 0
total_distri[total_distri > 29.7] = 1

total_module = total_enrich_df_reshape_log_bi %>% 
  select(STAD,LIHC,BRCA,CESC,LGG)

total_module[total_module <= 20] = 0
total_module[total_module > 20] = 1

total_other = total_enrich_df_reshape_log_bi %>% 
  select(LUAD,LUSC)
# summary(total_other)
total_other[total_other <= 50] = 0
total_other[total_other > 50] = 1

total_enrich_rebuild = cbind(total_distri,total_module,total_other)

setwd("~/nas/04.Results/cancerhallmark/")

svglite::svglite(filename = "suppl_fig_cancerhallmark_hyper_bi.svg")

fig_hall = pheatmap(total_enrich_rebuild %>% as.matrix() %>% t() ,
         # color = c( colorRampPalette(c("#F9FEFE", "black"))(100)),
         # color = inferno(50),
         scale = "none",
         cluster_cols = T,
         cluster_rows = T,
         border_color = NA,
         fontfamily = "Helvetica",
         fontfamily_row = "Helvetica",
         fontfamily_col = "Helvetica",
         fontface = "bold",
         fontface_row = "bold",
         fontface_col = "bold",
         fontsize_row = 10,  
         fontsize_col = 10,
         clustering_method = "average",
         # cutree_rows = 3,
         # cutree_cols = 4
)

print(fig_hall)
dev.off()
