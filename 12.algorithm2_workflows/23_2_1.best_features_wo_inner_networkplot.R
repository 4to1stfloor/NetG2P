library(dplyr)
library(RColorBrewer)
library(readxl)
library(openxlsx)
library(igraph)
library(ggforce)
library(ggplot2)
library(graphlayouts)
library(ggraph)

filepath = "/home/seokwon/nas/"
ref_path = paste0(filepath, "99.reference/")

Cancerlist = dir(paste0(filepath, "/00.data/filtered_TCGA/"))
surv_total_results = read_xlsx("~/nas/04.Results/Total_results_survpval2.xlsx")
# setwd("~/nas/04.Results/short_long/ttest")

total_features = list()
# Cancerlist = Cancerlist[!Cancerlist %in% c("04.TCGA-CESC","18.TCGA-LUAD","26.TCGA-LUSC","29.TCGA-LGG" ,"30.TCGA-BRCA","34.TCGA-COADREAD")]
# Cancerlist = Cancerlist[!Cancerlist %in% c("10.TCGA-BLCA","11.TCGA-STAD", "19.TCGA-LIHC","24.TCGA-OV" ,"32.TCGA-UCEC" ,"35.TCGA-KIDNEY" )]

# for all
for (num_CancerType in Cancerlist) {
  
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  Cancernum = gsub('[.-]','',gsub('[a-zA-Z]','', num_CancerType))
  # call input
  
  cancer_bf = read.csv(paste0(filepath,"04.Results/bestfeatures/",CancerType, "_best_features.csv"))
  
  cancer_bf_cut = cancer_bf[1:surv_total_results[which(surv_total_results$CancerType == CancerType),]$num_of_features,]
  total_features[[CancerType]] = cancer_bf_cut$variable
  
}

total_common = list()

second_Cancerlist = Cancerlist

for (num_CancerType in Cancerlist) {
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  second_Cancerlist = second_Cancerlist[-1]
  
  for (second_cancertype in second_Cancerlist) {
    second_CancerType = gsub('[.]','',gsub('\\d','', second_cancertype))
    
    total_common[[paste0(CancerType,'-',second_CancerType)]] = intersect(total_features[[CancerType]] , total_features[[second_CancerType]])
  }
  
}

for (num_CancerType in Cancerlist) {
  
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  Cancernum = gsub('[.-]','',gsub('[a-zA-Z]','', num_CancerType))
  # call input
  
  cancer_bf = read.csv(paste0(filepath,"04.Results/bestfeatures/",CancerType, "_best_features.csv"))
  
  cancer_bf_cut = cancer_bf[1:surv_total_results[which(surv_total_results$CancerType == CancerType),]$num_of_features,]
  cancer_bf_cut$variable
  
  nodes = paste0(CancerType,"_",cancer_bf_cut$variable )
  graph <- make_full_graph(length(nodes), directed = FALSE)
  # graph <- make_ring(length(nodes), directed = FALSE)
  # graph <- make_empty_graph(length(nodes), directed = FALSE) 
  
  
  V(graph)$name = nodes
  
  if (CancerType == "TCGA-CESC") {
    total_network = graph
  }else {
    total_network = graph.union( total_network , graph )
  }
  
  # assign(paste0(CancerType,"_network"), graph )
  
}

second_Cancerlist = Cancerlist
for (num_CancerType in Cancerlist) {
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  second_Cancerlist = second_Cancerlist[-1]
  
  for (second_cancertype in second_Cancerlist) {
    
    second_CancerType = gsub('[.]','',gsub('\\d','', second_cancertype))
    
    for (commonfeature in total_common[[paste0(CancerType, "-",second_CancerType)]]) {
      # tmp_ind = grepl(commonfeature, V(total_network)[1:surv_total_results[which(surv_total_results$CancerType == "TCGA-CESC"),]$num_of_features]$name)
      tmp_edge1 = V(total_network)$name[which(grepl(CancerType , V(total_network)$name))]
      tmp_edge2 = V(total_network)$name[which(grepl(second_CancerType , V(total_network)$name))]
      
      tmp_ind1 = which(commonfeature == gsub(".*_", "",tmp_edge1))
      tmp_ind2 = which(commonfeature == gsub(".*_", "",tmp_edge2))
      
      # tmp_network = add_edges(tmp_network , c(V(total_network)[1:surv_total_results[which(surv_total_results$CancerType == "TCGA-CESC"),]$num_of_features]$name[tmp_ind1] , 
      #                                         V(total_network)[1:surv_total_results[which(surv_total_results$CancerType == "TCGA-BLCA"),]$num_of_features]$name[tmp_ind2] ))
      total_network = add_edges(total_network, c(tmp_edge1[tmp_ind1] , 
                                                 tmp_edge2[tmp_ind2] ))
      # print(second_CancerType)
    }
    
  }
  
}

for (num_CancerType in Cancerlist) {
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  second_Cancerlist = second_Cancerlist[-1]
  
  for (second_cancertype in second_Cancerlist) {
    
    second_CancerType = gsub('[.]','',gsub('\\d','', second_cancertype))
    
    for (commonfeature in total_common[[paste0(CancerType, "-",second_CancerType)]]) {
      # tmp_ind = grepl(commonfeature, V(total_network)[1:surv_total_results[which(surv_total_results$CancerType == "TCGA-CESC"),]$num_of_features]$name)
      tmp_edge1 = V(total_network)$name[which(grepl(CancerType , V(total_network)$name))]
      tmp_edge2 = V(total_network)$name[which(grepl(second_CancerType , V(total_network)$name))]
      
      tmp_ind1 = which(commonfeature == gsub(".*_", "",tmp_edge1))
      tmp_ind2 = which(commonfeature == gsub(".*_", "",tmp_edge2))
      
      # tmp_network = add_edges(tmp_network , c(V(total_network)[1:surv_total_results[which(surv_total_results$CancerType == "TCGA-CESC"),]$num_of_features]$name[tmp_ind1] , 
      #                                         V(total_network)[1:surv_total_results[which(surv_total_results$CancerType == "TCGA-BLCA"),]$num_of_features]$name[tmp_ind2] ))
      total_network = add_edges(total_network, c(tmp_edge1[tmp_ind1] , 
                                                 tmp_edge2[tmp_ind2] ))
      # print(second_CancerType)
    }
    
  }
  
}



V(total_network)$grp =  as.character(c(
  rep("TCGA-CESC", surv_total_results$num_of_features[1]),
  rep("TCGA-BLCA", surv_total_results$num_of_features[2]),
  rep("TCGA-STAD", surv_total_results$num_of_features[3]),
  rep("TCGA-LUAD", surv_total_results$num_of_features[4]),
  rep("TCGA-LIHC", surv_total_results$num_of_features[5]),
  rep("TCGA-OV", surv_total_results$num_of_features[6]),
  rep("TCGA-LUSC", surv_total_results$num_of_features[7]),
  rep("TCGA-LGG", surv_total_results$num_of_features[8]),
  rep("TCGA-BRCA", surv_total_results$num_of_features[9]),
  rep("TCGA-UCEC", surv_total_results$num_of_features[10]),
  rep("TCGA-COADREAD", surv_total_results$num_of_features[11]),
  rep("TCGA-KIDNEY", surv_total_results$num_of_features[12])
))
total_common
# Simplify the graph
total_network <- simplify(total_network)

# Convert edge data to a data frame
edge_data <- get.data.frame(total_network, what = "edges")

marked_edges_data = data.frame()
for (num in 1:length(names(total_common))) {

  
  tmp_split_common = unlist(strsplit(names(total_common)[num], "-", fixed = TRUE))
  first_cancer = paste(tmp_split_common[1], tmp_split_common[2], sep = "-")
  second_cancer <- paste(tmp_split_common[3], tmp_split_common[4], sep = "-")

  for (dup_feature in total_common[[num]]) {
    first_feature = paste0(first_cancer,"_",dup_feature)
    second_feature = paste0(second_cancer,"_",dup_feature)
    
    tmp_common_edge = edge_data[which(edge_data$from == first_feature),][which(edge_data[which(edge_data$from == first_feature),]$to == second_feature),]
    marked_edges_data = rbind(marked_edges_data, tmp_common_edge)
  } 
  
}

marked_edges_data <- data.frame(marked_edges_data, row.names = NULL)

hidden_edges_index = as.numeric(rownames(edge_data)[!rownames(edge_data) %in% rownames(marked_edges_data)])
total_modify_network = delete_edges(total_network, hidden_edges_index)

bb <- layout_as_backbone(total_network, keep = 0.1)
bb2 <- layout_as_backbone(total_network, keep = 0.8)

## need to test 
# E(total_modify_network)$col = F
# E(total_modify_network)$col[bb$backbone] = T

# ggraph(total_modify_network,
#        layout = "manual",
#        x = bb$xy[, 1],
#        y = bb$xy[, 2]) +
#   geom_edge_link0( aes(alpha = col), width = 0.2) +
#   geom_node_point(aes(fill = grp), shape = 21, size = 3) +
#   scale_color_brewer(palette = "Set3") +
#   scale_fill_brewer(palette = "Set3") +
#   scale_edge_color_manual(values = c(rgb(0, 0, 0, 0.3), rgb(0, 0, 0, 1))) +
#   theme_graph()

# fin
ggraph(total_modify_network,
       layout = "manual",
       x = bb$xy[, 1],
       y = bb$xy[, 2]) +
  geom_edge_link0( aes(alpha = 0.5), width = 0.2) +
  geom_node_point(aes(fill = grp), shape = 21, size = 5) +
  scale_color_brewer(palette = "Set3") +
  scale_fill_brewer(palette = "Set3") +
  scale_edge_color_manual(values = c(rgb(0, 0, 0, 0.3), rgb(0, 0, 0, 1))) +
  theme_graph()

ggraph(total_modify_network,
       layout = "manual",
       x = bb2$xy[, 1],
       y = bb2$xy[, 2]) +
  geom_edge_link0( aes(alpha = 0.5), width = 0.2) +
  geom_node_point(aes(fill = grp), shape = 21, size = 5) +
  scale_color_brewer(palette = "Set3") +
  scale_fill_brewer(palette = "Set3") +
  scale_edge_color_manual(values = c(rgb(0, 0, 0, 0.3), rgb(0, 0, 0, 1))) +
  theme_graph()

# test
degree_score <- degree(total_modify_network, mode = 'total')
betw_score <- betweenness(total_modify_network)


# edge color
ggraph(total_network,
       layout = "manual",
       x = bb$xy[, 1],
       y = bb$xy[, 2]) +
  
  geom_edge_link0(
    aes(filter = node1.grp == "TCGA-CESC" & node2.grp == "TCGA-CESC"),
    alpha = 0.3,
    edge_linewidth = 0.1,
    edge_colour = "red"
  ) +
  
  geom_edge_link0(
    aes(filter = node1.grp == "TCGA-BLCA" & node2.grp == "TCGA-BLCA"),
    alpha = 0.3,
    edge_linewidth = 0.1,
    edge_colour = "blue"
  ) +
  
  geom_edge_link0(
    aes(filter = node1.grp == "TCGA-STAD" & node2.grp == "TCGA-STAD"),
    alpha = 0.3,
    edge_linewidth = 0.1,
    edge_colour = "yellow"
  ) +
  geom_edge_link0(
    aes(filter = node1.grp == "TCGA-LUAD" & node2.grp == "TCGA-LUAD"),
    alpha = 0.3,
    edge_linewidth = 0.1,
    edge_colour = "gold"
  ) +
  geom_edge_link0(
    aes(filter = node1.grp == "TCGA-LIHC" & node2.grp == "TCGA-LIHC"),
    alpha = 0.3,
    edge_linewidth = 0.1,
    edge_colour = "green"
  ) +
  geom_edge_link0(
    aes(filter = node1.grp == "TCGA-OV" & node2.grp == "TCGA-OV"),
    alpha = 0.3,
    edge_linewidth = 0.1,
    edge_colour = "purple"
  ) +
  geom_edge_link0(
    aes(filter = node1.grp == "TCGA-LUSC" & node2.grp == "TCGA-LUSC"),
    alpha = 0.3,
    edge_linewidth = 0.1,
    edge_colour = "orange"
  ) +
  geom_edge_link0(
    aes(filter = node1.grp == "TCGA-LGG" & node2.grp == "TCGA-LGG"),
    alpha = 0.3,
    edge_linewidth = 0.1,
    edge_colour = "aliceblue"
  ) +
  geom_edge_link0(
    aes(filter = node1.grp == "TCGA-BRCA" & node2.grp == "TCGA-BRCA"),
    alpha = 0.3,
    edge_linewidth = 0.1,
    edge_colour = "darkkhaki"
  ) +
  geom_edge_link0(
    aes(filter = node1.grp == "TCGA-UCEC" & node2.grp == "TCGA-UCEC"),
    alpha = 0.3,
    edge_linewidth = 0.1,
    edge_colour = "darkslategray2"
  ) +
  geom_edge_link0(
    aes(filter = node1.grp == "TCGA-COADREAD" & node2.grp == "TCGA-COADREAD"),
    alpha = 0.3,
    edge_linewidth = 0.1,
    edge_colour = "burlywood1"
  ) +
  geom_edge_link0(
    aes(filter = node1.grp == "TCGA-KIDNEY" & node2.grp == "TCGA-KIDNEY"),
    alpha = 0.3,
    edge_linewidth = 0.1,
    edge_colour = "brown3"
  ) +
  
  geom_edge_link0(
    aes(filter = (node1.grp != node2.grp)),
    alpha = 0.3,
    edge_linewidth = 0.1,
    edge_colour = "black"
  ) +
  
  geom_node_point(aes(fill = grp), shape = 21, size = 5) +
  scale_color_brewer(palette = "Set3") +
  scale_fill_brewer(palette = "Set3") +
  scale_edge_color_manual(values = c(rgb(0, 0, 0, 0.3), rgb(0, 0, 0, 1))) +
  theme_graph()

ggraph(total_network,
       layout = "manual",
       x = bb2$xy[, 1],
       y = bb2$xy[, 2]) +
  
  geom_edge_link0(
    aes(filter = node1.grp == "TCGA-CESC" & node2.grp == "TCGA-CESC"),
    alpha = 0.3,
    edge_linewidth = 0.1,
    edge_colour = "red"
  ) +

  geom_edge_link0(
    aes(filter = node1.grp == "TCGA-BLCA" & node2.grp == "TCGA-BLCA"),
    alpha = 0.3,
    edge_linewidth = 0.1,
    edge_colour = "blue"
  ) +
  
  geom_edge_link0(
    aes(filter = node1.grp == "TCGA-STAD" & node2.grp == "TCGA-STAD"),
    alpha = 0.3,
    edge_linewidth = 0.1,
    edge_colour = "yellow"
  ) +
  geom_edge_link0(
    aes(filter = node1.grp == "TCGA-LUAD" & node2.grp == "TCGA-LUAD"),
    alpha = 0.3,
    edge_linewidth = 0.1,
    edge_colour = "gold"
  ) +
  geom_edge_link0(
    aes(filter = node1.grp == "TCGA-LIHC" & node2.grp == "TCGA-LIHC"),
    alpha = 0.3,
    edge_linewidth = 0.1,
    edge_colour = "green"
  ) +
  geom_edge_link0(
    aes(filter = node1.grp == "TCGA-OV" & node2.grp == "TCGA-OV"),
    alpha = 0.3,
    edge_linewidth = 0.1,
    edge_colour = "purple"
  ) +
  geom_edge_link0(
    aes(filter = node1.grp == "TCGA-LUSC" & node2.grp == "TCGA-LUSC"),
    alpha = 0.3,
    edge_linewidth = 0.1,
    edge_colour = "orange"
  ) +
  geom_edge_link0(
    aes(filter = node1.grp == "TCGA-LGG" & node2.grp == "TCGA-LGG"),
    alpha = 0.3,
    edge_linewidth = 0.1,
    edge_colour = "aliceblue"
  ) +
  geom_edge_link0(
    aes(filter = node1.grp == "TCGA-BRCA" & node2.grp == "TCGA-BRCA"),
    alpha = 0.3,
    edge_linewidth = 0.1,
    edge_colour = "darkkhaki"
  ) +
  geom_edge_link0(
    aes(filter = node1.grp == "TCGA-UCEC" & node2.grp == "TCGA-UCEC"),
    alpha = 0.3,
    edge_linewidth = 0.1,
    edge_colour = "darkslategray2"
  ) +
  geom_edge_link0(
    aes(filter = node1.grp == "TCGA-COADREAD" & node2.grp == "TCGA-COADREAD"),
    alpha = 0.3,
    edge_linewidth = 0.1,
    edge_colour = "burlywood1"
  ) +
  geom_edge_link0(
    aes(filter = node1.grp == "TCGA-KIDNEY" & node2.grp == "TCGA-KIDNEY"),
    alpha = 0.3,
    edge_linewidth = 0.1,
    edge_colour = "brown3"
  ) +
  
  geom_edge_link0(
    aes(filter = (node1.grp != node2.grp)),
    alpha = 0.6,
    edge_linewidth = 0.1,
    edge_colour = "black"
  ) +
  
  geom_node_point(aes(fill = grp), shape = 21, size = 5) +
  scale_color_brewer(palette = "Set3") +
  scale_fill_brewer(palette = "Set3") +
  scale_edge_color_manual(values = c(rgb(0, 0, 0, 0.3), rgb(0, 0, 0, 1))) +
  theme_graph()
