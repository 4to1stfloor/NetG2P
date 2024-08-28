library(readxl)
library(tidyverse)
library(svglite)
library(ggplot2)
library(ggsignif)
library(openxlsx)
library(ggpubr)
library(ggsignif)
library(ggraph)
library(tidygraph)

filepath = "/home/seokwon/nas/"
ref_path = paste0(filepath, "99.reference/")

Cancerlist = dir(paste0(filepath, "/00.data/filtered_TCGA/"))

link_genes = readRDS(paste0(ref_path, "/KEGG_pathway_shared_genes.rds"))
single_genes = readRDS(paste0(ref_path, "/Kegg_pathway_genes.rds"))
link_genes$n_genes = NULL
link_genes_filtered = link_genes[which(link_genes$shared_genes != ""),]

link_genes_filtered <- link_genes_filtered %>%
  mutate(shared_genes = strsplit(shared_genes, ",")) %>%
  unnest(cols = shared_genes)
link_genes_filtered_df = as.data.frame(link_genes_filtered)
colnames(link_genes_filtered_df) = colnames(single_genes)

Cancerlist = Cancerlist[c(2,3,5,6,10)]
num_CancerType = "11.TCGA-STAD"

for (num_CancerType in Cancerlist) {
  
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  Cancername = gsub('TCGA-' , '', CancerType)
  
  # call input

  gc_TCGA = readRDS(paste0("~/nas/04.Results/short_long/", CancerType,"_critical_features_short_long.rds"))
 
  # network for 
  cf_sl = read.csv(paste0(filepath, "04.Results/short_long/ttest_common/",CancerType,"_critical_features_short_long_common.csv"))
  
  edges = data.frame(matrix(ncol = 4))
  colnames(edges) = c("from","to","group","variable")
  # critical_features = "P16"
  
  for (critical_features in cf_sl %>% pull(variable)) {
    # nodes = unique(c(nodes , paste0("P",strsplit(critical_features, "P")[[1]][2]) , paste0("P",strsplit(critical_features, "P")[[1]][3])))
    
    if (is.na(strsplit(critical_features, "P")[[1]][3] )) {
      edges[critical_features,"from"] = paste0("P",strsplit(critical_features, "P")[[1]][2])
      edges[critical_features,"to"] = paste0("P",strsplit(critical_features, "P")[[1]][2])
      edges[critical_features,"group"] = cf_sl %>% filter(variable == critical_features) %>% select(classification) %>% pull()
      edges[critical_features,"variable"] = critical_features
      
    } else {
      edges[critical_features,"from"] = paste0("P",strsplit(critical_features, "P")[[1]][2]) 
      edges[critical_features,"to"] = paste0("P",strsplit(critical_features, "P")[[1]][3])
      edges[critical_features,"group"] = cf_sl %>% filter(variable == critical_features) %>% select(classification) %>% pull()
      edges[critical_features,"variable"] = critical_features
      
    }
    
  }
  
  edges = edges[-1,]
  # spe_features_drug
  rownames(edges) = NULL
  
  edges = left_join(edges, cf_sl , by = "variable") %>% select(from, to, group, variable , classification)
  
  set.seed(1234)
  
  critical_network = edges %>%
    as_tbl_graph() %N>% 
    left_join(., cf_sl %>% select(variable, classification, relative_importance), by = c("name" = "variable")) %E>%
    left_join(., cf_sl %>% select(variable, relative_importance), by = c("variable" = "variable")) %N>%
    mutate(classification = ifelse(is.na(classification), "", classification),
           color = case_when(classification == "short"  ~ "exist",
                             classification == "long"  ~ "exist",
                             classification == "common"  ~ "exist",
                             .default = "non"),
           node_weight = relative_importance,
           node_weight = ifelse(is.na(node_weight), 0.1, node_weight)) %E>%
    mutate(color = case_when(classification == "short"  ~ "exist",
                             classification == "long"  ~ "exist",
                             classification == "common"  ~ "exist",
                             .default = "non"),
           edge_weight = relative_importance) %>% 
    ggraph(layout = "nicely") +
    geom_edge_link(aes(color = color,
                       edge_width = edge_weight)) +
    scale_edge_width(range = c(0.5,2)) +
    # scale_edge_alpha(range = c(0.3,1)) +
    geom_node_point(aes(color = color,
                        size = node_weight)) +               # 노드 크기
    geom_node_text(aes(label = name),         # 텍스트 표시
                   repel = T,                 # 노드밖 표시
                   size = 5) +  
    scale_size(range = c(3, 8)) +
    # scale_color_manual(values = c("short" = "blue", "common" = "red", "long" = "green")) +
    scale_color_manual(values = c("exist" = "black",
                                  "non"= "grey")) +
    # scale_color_manual(values = c("short" = "#E41A1C", "long" = "#4DAF4A")) +
    scale_edge_colour_manual(values =c("exist" = "black",
                                       "non"= "grey")) +
    # geom_segment(aes(x = -3, y = -1.5, xend = 0.1, yend = 0.1),
    #                           arrow = arrow(length = unit(0.1, "cm"))) +
    theme_graph() 

  setwd("~/nas/04.Results/critical_features/network/")
  if (Cancername %in% c("BLCA", "OV","UCEC")) {
    ggsave(critical_network , filename = paste0(Cancername,"_critical_network.svg"), width = 40, height = 40)
  } else {
    ggsave(critical_network , filename = paste0(Cancername,"_critical_network.svg"))
  }


  # lihc = readRDS(paste0("~/nas/04.Results/short_long/", CancerType,"_critical_features_short_long_with_drug.rds"))
}

# critical_network = edges %>%
#   as_tbl_graph() %N>% 
#   left_join(., cf_sl %>% select(variable, classification), by = c("name" = "variable")) %N>%
#   mutate(classification = ifelse(is.na(classification), "", classification)) %N>%
#   mutate(color = case_when(classification == "short"  ~ "exist",
#                            classification == "long"  ~ "exist",
#                            classification == "common"  ~ "exist",
#                            .default = "")) %E>%
#   mutate(color = case_when(classification == "short"  ~ "exist",
#                            classification == "long"  ~ "exist",
#                            classification == "common"  ~ "exist",
#                            .default = "")) %E>%
#   mutate(edge_weight = case_when(color == "short" ~ 1,
#                                  color == "long" ~ 1,
#                                  .default = 0.3)) %N>% 
#   mutate(node_alpha = case_when(classification %in% c("short", "long","common") ~ 1,
#                                 .default = 0.5)) %E>%
#   mutate(edge_alpha = case_when(classification %in% c("short", "long","common") ~ 1,
#                                 .default = 0.5)) %>% 
#   ggraph(layout = "nicely") +
#   geom_edge_link(aes(
#     edge_alpha = edge_alpha,
#     edge_width = edge_weight)) +
#   scale_edge_width(range = c(0.5,1.5)) +   
#   scale_edge_alpha(range = c(0.3,1)) +
#   geom_node_point(aes(
#     alpha = node_alpha),
#     # color = "red",            # 노드 색깔
#     size = 5) +               # 노드 크기
#   geom_node_text(aes(label = name),         # 텍스트 표시
#                  repel = T,                 # 노드밖 표시
#                  size = 5) +  
#   # scale_color_manual(values = c("short" = "blue", "common" = "red", "long" = "green")) +
#   # scale_color_manual(values = c("short" = "#E41A1C", 
#   #                               "long" = "#4DAF4A",
#   #                               "common"= "grey")) +
#   # # scale_color_manual(values = c("short" = "#E41A1C", "long" = "#4DAF4A")) +
#   # scale_edge_colour_manual(values =c("short" = "#E41A1C", 
#   #                                    "long" = "#4DAF4A",
#   #                                    "common"= "grey")) +
#   # geom_segment(aes(x = -3, y = -1.5, xend = 0.1, yend = 0.1),
#   #                           arrow = arrow(length = unit(0.1, "cm"))) +
#   theme_graph() 
