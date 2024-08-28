library(dplyr)
library(tidyr)
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
# num_CancerType = "04.TCGA-CESC"

for (num_CancerType in Cancerlist) {
  
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  Cancername = gsub('TCGA-' , '', CancerType)
  if (Cancername == "LGG") {
    next
  }
  depmap_sig = read.csv(paste0(filepath, "04.Results/drug/depmap/",Cancername, "_depmap_sig_total.csv"))
  cf = read.csv(paste0("~/nas/04.Results/critical_features/eachcancer/",CancerType,"_critical_features.csv"))
  
  link_genes_filtered_df %>% filter(Pathway %in% cf$variable) %>% filter(Genes %in% depmap_sig$Genes)
  single_genes %>% filter(Pathway %in% cf$variable) %>% filter(Genes %in% depmap_sig$Genes)

  depmap_sig = depmap_sig %>% select(-X)
  colnames(depmap_sig) = c("Genes", "delta_long_to_short", "sig", "specific")
  
  depmap_sig_with_link = left_join(depmap_sig, link_genes_filtered_df , by = "Genes")
  
  colnames(depmap_sig_with_link) = c("Genes", "delta_long_to_short", "sig", "specific" , "pathwaylink")
  depmap_sig_with_cf = left_join(depmap_sig_with_link, single_genes , by = "Genes", relationship = "many-to-many")
  
  depmap_sig_with_cf_filt = depmap_sig_with_cf %>% filter(pathwaylink %in% cf$variable | Pathway %in% cf$variable)
  depmap_sig_tmp = depmap_sig_with_cf_filt %>% mutate(pathwaylink = case_when(pathwaylink %in% cf$variable ~ pathwaylink,
                                                                                       .default = NA),
                                                               Pathway = case_when(Pathway %in% cf$variable ~ Pathway,
                                                                                   .default = NA))

  if (all.equal(unique(depmap_sig$Genes) , unique(depmap_sig_tmp$Genes))) {
    depmap_sig_tmp_filt = depmap_sig_tmp[!duplicated(depmap_sig_tmp),]
  }
  
  write.csv(depmap_sig_tmp_filt , paste0(filepath, "04.Results/drug/depmap/",Cancername, "_depmap_sig_w_cf.csv"))
}
