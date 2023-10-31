library(tidyverse)
library(ggplot2)
library(ggridges)
library(reshape2)
library(readxl)
library(org.Hs.eg.db)

library(ggsignif)
library(ggrepel)
library(tidygraph)

library(ggraph)
library(showtext)

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

cancer_specific_critical_features = read.csv("~/nas/04.Results/critical_features/unique_short_long_critical_features_for_SW_ver2.csv")

# for fic
fig_path = "~/nas/04.Results/critical_features/cancer_specific_critical_features/"
if(!dir.exists(fig_path)){
  dir.create(fig_path)
  print(paste0("Created folder: ", fig_path))
} else {
  print(paste0("Folder already exists: ", fig_path))
}
setwd(fig_path)

## Read Dataset 

# num_CancerType = "10.TCGA-BLCA"
Cancerlist = dir(paste0(filepath, "/00.data/filtered_TCGA/"))
Cancerlist = Cancerlist[c(-7,-11,-12)]
colnames(cancer_specific_critical_features )  = c(colnames(cancer_specific_critical_features )[1:6] , gsub('[.]','',gsub('\\d','', Cancerlist)))
# num_CancerType = "04.TCGA-CESC"
for (num_CancerType in Cancerlist) {
  
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
 
  # call input
  # just extract only specific critical features 
  
  if (sum(cancer_specific_critical_features$which_cancer == CancerType) == 0) {
    next
  } else {
    
    cancer_critical_features = cancer_specific_critical_features[which(cancer_specific_critical_features$which_cancer == CancerType),]
    cancer_critical_features = cancer_critical_features[which(cancer_critical_features[,CancerType] != "common"),]
    total_link_genes = link_genes_filtered_df[which(link_genes_filtered_df$Pathway %in% cancer_critical_features$features),]$Genes
    total_single_genes = single_genes[which(single_genes$Pathway %in% cancer_critical_features$features),]$Genes
    interested_gene_set = c(total_link_genes,total_single_genes)
  }
  
  nodes = c()
  edges = data.frame(matrix(ncol = 2))
  colnames(edges) = c("from","to")

  for (critical_features in cancer_critical_features$variable) {
    nodes = unique(c(nodes , paste0("P",strsplit(critical_features, "P")[[1]][2]) , paste0("P",strsplit(critical_features, "P")[[1]][3])))
    
    if (is.na(strsplit(critical_features, "P")[[1]][3] )) {
      edges[critical_features,"from"] = paste0("P",strsplit(critical_features, "P")[[1]][2]) 
      edges[critical_features,"to"] = paste0("P",strsplit(critical_features, "P")[[1]][2])
    } else {
      edges[critical_features,"from"] = paste0("P",strsplit(critical_features, "P")[[1]][2]) 
      edges[critical_features,"to"] = paste0("P",strsplit(critical_features, "P")[[1]][3])
    }
    
  }
  edges = edges[-1,]
  rownames(edges) = NULL
  
  # 
  
  png(filename = paste0(CancerType,"_raw_network_critical_features.png"),
      width = 10, height = 10,  units = "cm" ,pointsize = 12,
      bg = "white", res = 1200, family = "")
  
  raw_network = ggraph(edges %>% as_tbl_graph(), layout = "fr") + # 레이아웃
    geom_edge_link(color = "gray50",          # 엣지 색깔
                   alpha = 0.5) +             # 엣지 명암
    
    geom_node_point(color = "lightcoral",     # 노드 색깔
                    size = 5) +               # 노드 크기
    
    geom_node_text(aes(label = name),         # 텍스트 표시
                   repel = T,                 # 노드밖 표시
                   size = 5,                  # 텍스트 크기
                   family = "nanumgothic") +  # 폰트
    
    theme_graph()                             # 배경 삭제
  
  print(raw_network)
  
  dev.off()
  
  png(filename = paste0(CancerType,"_centrality_network_critical_features.png"),
      width = 10, height = 10,  units = "cm" ,pointsize = 12,
      bg = "white", res = 1200, family = "")
  
  central_network = ggraph( edges %>%
                              as_tbl_graph(directed = F) %>%
                              mutate(centrality = centrality_degree(),        # 연결 중심성
                                     group = as.factor(group_louvain())), layout = "fr") + # 레이아웃
    geom_edge_link(color = "gray50",          # 엣지 색깔
                   alpha = 0.5) +             # 엣지 명암
    
    geom_node_point(aes(size = centrality,    # 노드 크기
                        color = group),       # 노드 색깔
                    show.legend = F) +        # 범례 삭제
    scale_size(range = c(3, 10)) +      
    
    geom_node_text(aes(label = name),         # 텍스트 표시
                   repel = T,                 # 노드밖 표시
                   size = 5,                  # 텍스트 크기
                   family = "nanumgothic") +  # 폰트
    
    theme_graph()                             # 배경 삭제
  
  print(central_network)
  
  dev.off()
  
}
