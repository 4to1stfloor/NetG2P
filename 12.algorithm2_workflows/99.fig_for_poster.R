library(readxl)
library(tidyverse)
filepath = "/home/seokwon/nas/"
ref_path = paste0(filepath, "99.reference/")

Cancerlist = dir(paste0(filepath, "/00.data/filtered_TCGA/"))
surv_total_results = read_xlsx("~/nas/04.Results/Total_results_survpval2.xlsx")
# setwd("~/nas/04.Results/short_long/ttest")
Cancerlist = Cancerlist[c(-7,-11,-12)]
total_features = list()

# for all
for (num_CancerType in Cancerlist) {
  
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  Cancernum = gsub('[.-]','',gsub('[a-zA-Z]','', num_CancerType))
  # call input
  
  cancer_bf = read.csv(paste0(filepath,"04.Results/bestfeatures/",CancerType, "_critical_features.csv"))
  
  # cancer_bf_cut = cancer_bf[1:surv_total_results[which(surv_total_results$CancerType == CancerType),]$num_of_features,]
  total_features[[CancerType]] = cancer_bf$variable
  
}

total_spe = data.frame()
num_CancerType = "19.TCGA-LIHC"
for (num_CancerType in Cancerlist) {
  
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  
  cancer_cf = total_features[[CancerType]]
 
  cancer_slc = read.csv(paste0("~/nas/04.Results/short_long/ttest_common/",CancerType,"_critical_features_short_long_common.csv"))
  tmp_spe = cancer_slc[which(cancer_slc$variable %in% cancer_cf),]
  
  tmp_sep_filt = tmp_spe %>% dplyr::select(variable, classification) %>% 
    mutate(long_rat = case_when( classification == "long" ~ 1,.default = 0)) %>% 
    mutate(short_rat = case_when( classification == "short" ~ 1,.default = 0)) %>%
    mutate(common_rat = case_when( classification == "common" ~ 1,.default = 0)) %>% 
    dplyr::select(variable, long_rat, short_rat, common_rat) %>%
    mutate(which_cancer = CancerType) 
  
  total_spe = rbind(total_spe , tmp_sep_filt)
  
  remove(tmp_sep_filt)
  
}



cancer_critical_features = total_spe %>% filter(common_rat != 1)

# total_link_genes = link_genes_filtered_df[which(link_genes_filtered_df$Pathway %in% cancer_critical_features$features),]$Genes
# total_single_genes = single_genes[which(single_genes$Pathway %in% cancer_critical_features$features),]$Genes
# interested_gene_set = c(total_link_genes,total_single_genes)


long_cf = cancer_critical_features %>% filter(long_rat == 1)
short_cf = cancer_critical_features %>% filter(short_rat == 1)
common_cf = total_spe %>% filter(common_rat == 1)

nodes = c()
edges = data.frame(matrix(ncol = 2))
colnames(edges) = c("from","to")

for (critical_features in short_cf$variable) {
  
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

nodes = c()
long_edges = data.frame(matrix(ncol = 2))
colnames(long_edges) = c("from","to")

for (critical_features in long_cf$variable) {
  
  nodes = unique(c(nodes , paste0("P",strsplit(critical_features, "P")[[1]][2]) , paste0("P",strsplit(critical_features, "P")[[1]][3])))
  
  if (is.na(strsplit(critical_features, "P")[[1]][3] )) {
    long_edges[critical_features,"from"] = paste0("P",strsplit(critical_features, "P")[[1]][2]) 
    long_edges[critical_features,"to"] = paste0("P",strsplit(critical_features, "P")[[1]][2])
  } else {
    long_edges[critical_features,"from"] = paste0("P",strsplit(critical_features, "P")[[1]][2]) 
    long_edges[critical_features,"to"] = paste0("P",strsplit(critical_features, "P")[[1]][3])
  }
  
}
long_edges = long_edges[-1,]
rownames(long_edges) = NULL

nodes = c()
common_edges = data.frame(matrix(ncol = 2))
colnames(common_edges) = c("from","to")

for (critical_features in common_cf$variable) {

  nodes = unique(c(nodes , paste0("P",strsplit(critical_features, "P")[[1]][2]) , paste0("P",strsplit(critical_features, "P")[[1]][3])))
  
  if (is.na(strsplit(critical_features, "P")[[1]][3] )) {
    common_edges[critical_features,"from"] = paste0("P",strsplit(critical_features, "P")[[1]][2]) 
    common_edges[critical_features,"to"] = paste0("P",strsplit(critical_features, "P")[[1]][2])
  } else {
    common_edges[critical_features,"from"] = paste0("P",strsplit(critical_features, "P")[[1]][2]) 
    common_edges[critical_features,"to"] = paste0("P",strsplit(critical_features, "P")[[1]][3])
  }
  
}
common_edges = common_edges[-1,]
rownames(common_edges) = NULL

edges$group = "short"
long_edges$group = "long"
common_edges$group = "common"

# total_edges = rbind(edges, long_edges)
total_edges = rbind(edges, long_edges,common_edges)

# library(GGally)
# library(network)
# library(sna)
library(ggplot2)
library(igraph)

count_node = data.frame()

# spe_nodes = "P18"
# spe_nodes = "P24"
for (spe_nodes in unique(c(total_edges$from, total_edges$to))) {
  
  tmp_edges = total_edges %>% 
    filter(from == spe_nodes | to == spe_nodes) %>% 
    filter(from != to)
  if (nrow(tmp_edges) == 0 ) {
    next
  }
  if (length(names(which(table(tmp_edges$group) == max(table(tmp_edges$group))))) != 1) {
    node_name = tmp_edges$group[1]
  } else {
    node_name = names(which(table(tmp_edges$group) == max(table(tmp_edges$group))))
  }
 
  count_node = rbind(count_node, data.frame(spe_nodes = spe_nodes, node_name = node_name))
}


count_node$node_name <- ifelse(is.na(count_node$node_name), "self", count_node$node_name)

i_total_g = graph.data.frame(total_edges, directed = FALSE)

for (a in names(V(i_total_g))) {

  V(i_total_g)[which(names(V(i_total_g)) == a)]$grp = count_node[which(count_node$spe_nodes == a ),]$node_name

}

library(ggforce)
library(igraph)
library(graphlayouts)
sort(betweenness(i_total_g),decreasing = F)
names(which(degree(i_total_g)== max(betweenness(i_total_g))))

# V(i_total_g)$grp
ggraph(i_total_g, layout = "nicely") + # 레이아웃
  geom_edge_link(aes(color = group),          # 엣지 색깔
                  alpha = 0.5) +             # 엣지 명암
  
  geom_node_point(aes(color = grp),     # 노드 색깔
                  size = 5) +               # 노드 크기
  
  geom_node_text(aes(label = name),         # 텍스트 표시
                 repel = T,                 # 노드밖 표시
                 size = 5) +  
  # scale_color_manual(values = c("short" = "blue", "common" = "red", "long" = "green")) +
  # scale_color_manual(values = c("short" = "#AD002AFF", "common" = "#925E9FFF", "long" = "#00468BFF")) +
  scale_color_manual(values = c("short" = "#AFCDFF", "common" = "#FBBEB9", "long" = "#93E1AB")) +
  theme_graph()                             # 배경 삭제
