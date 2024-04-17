library(dplyr)
library(tidyr)
library(ComplexHeatmap)
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

####
Cancerlist = Cancerlist[c(-11,-12)]
# num_CancerType = "04.TCGA-CESC"

total_counts = sum(sapply(cancerhallmark_geneset, function(x) sum(x != "")))
total_hallmark_hyper = data.frame()
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

  hallmark_hyper = data.frame()
  for (hallmark in colnames(cancerhallmark_geneset)) {
    
    tmp_hallmark = cancerhallmark_geneset %>% 
      select(all_of(hallmark)) 
    
    tmp_hallmark_filt = tmp_hallmark[tmp_hallmark != "" , ,drop = F]
    
    x = sum(tmp_hallmark_filt %>% pull(all_of(hallmark)) %in% total_critical_genes)
    m = nrow(tmp_hallmark_filt)
    n = total_counts - m
    k = length(total_critical_genes)

    # pvalue = phyper(x-1,m,n,k , lower.tail = F)
    pvalue_test = dhyper(x,m,n,k )
    
    # tmp_hyper = data.frame(hallmark = colnames(tmp_hallmark_filt), 
    #                        pvalue = pvalue)
    tmp_hyper = data.frame(hallmark = colnames(tmp_hallmark_filt), 
                           pvalue = pvalue_test)
    
    hallmark_hyper = rbind(hallmark_hyper, tmp_hyper)
  }
  
  rownames(hallmark_hyper) = hallmark_hyper$hallmark
  hallmark_hyper_filt = hallmark_hyper %>% select(pvalue)
  colnames(hallmark_hyper_filt) = Cancername
  hallmark_hyper_filt[hallmark_hyper_filt > 0.05] = 1
  hallmark_hyper_filt = -log(hallmark_hyper_filt)
  
  if (sum(hallmark_hyper_filt == 0) != 0) {
    hallmark_hyper_filt = (hallmark_hyper_filt - min(hallmark_hyper_filt)) / (max(hallmark_hyper_filt) - min(hallmark_hyper_filt))
  } else {
    hallmark_hyper_filt = (hallmark_hyper_filt - 0.1) / (max(hallmark_hyper_filt) - 0.1)
  }
  
  if (Cancername == "CESC") {
    total_hallmark_hyper = hallmark_hyper_filt
  } else {
    total_hallmark_hyper = cbind(total_hallmark_hyper , hallmark_hyper_filt)
  }

}

total_hallmark_hyper_filt = total_hallmark_hyper[names(sort(rowSums(total_hallmark_hyper), decreasing = T)),]
library(viridis)
library(svglite)

setwd("~/nas/04.Results/cancerhallmark/")

svglite::svglite(filename = "figure3D_cancerhallmark_hyper_maxnor.svg")

fig_hall = pheatmap(total_hallmark_hyper_filt %>% as.matrix() %>% t() ,
         # color = c( colorRampPalette(c("#F9FEFE", "black"))(100)),
         color = inferno(50),
         scale = "none",
         cluster_cols = F,
         border_color = NA,
         fontfamily = "Helvetica",
         fontfamily_row = "Helvetica",
         fontfamily_col = "Helvetica",
         fontface = "bold",
         fontface_row = "bold",
         fontface_col = "bold",
         fontsize_row = 10,  
         fontsize_col = 10
)


print(fig_hall)
dev.off()
