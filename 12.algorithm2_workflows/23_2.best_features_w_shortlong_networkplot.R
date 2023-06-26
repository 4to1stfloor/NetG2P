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


V(total_network)$grp =  as.character(c(rep("TCGA-CESC", surv_total_results$num_of_features[1]),
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

# Simplify the graph
total_network <- simplify(total_network)

bb <- layout_as_backbone(total_network, keep = 0.1)
# bb <- layout_as_backbone(total_network)

E(total_network)$col <- F
E(total_network)$col[bb$backbone] <- T

# Plot the graph with ggraph using the backbone layout
ggraph(total_network,
       layout = "manual",
       x = bb$xy[, 1],
       y = bb$xy[, 2]) +
  geom_edge_link0(aes(col = col), width = 0.2) +
  geom_node_point(aes(fill = grp), shape = 21, size = 3) +
  geom_mark_hull(
    aes(x, y, group = grp, fill = grp),
    concavity = 4,
    expand = unit(2, "mm"),
    alpha = 0.25
  ) +
  scale_color_brewer(palette = "Set3") +
  scale_fill_brewer(palette = "Set3") +
  scale_edge_color_manual(values = c(rgb(0, 0, 0, 0.3), rgb(0, 0, 0, 1))) +
  theme_graph()
  # + theme(legend.position = "none")


