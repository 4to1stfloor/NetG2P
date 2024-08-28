h2o.shutdown(prompt = F)

library(ggraph)
library(igraph)
library(graphlayouts)
library(wesanderson)
library(tidyverse)
library(readxl)
library(graphlayouts)
library(colorspace)
library(h2o)
library(caret)
library(tidygraph)

filepath = "/home/seokwon/"
ref_path = paste0(filepath, "99.reference/")
Cancerlist = dir(paste0(filepath, "/00.data/"))
# Cancerlist = Cancerlist[7:9]
Cancerlist = Cancerlist[c(-2,-5,-15)]
folder_name = "h2o_bias_pval_dual_cut_50"

localH2O = h2o.init(ip = "localhost", port = sample(x=00000:65536,size=1,replace=F), 
                    startH2O = TRUE,min_mem_size = "400G",
                    nthreads = 96,enable_assertions = FALSE)
pathway = as.data.frame(readxl::read_xlsx(paste0(ref_path ,"kegg_gene_set_w_cancer_hallmarks.xlsx")))

find_max_model <- function(vec) {
  vec <- vec[!grepl(".csv", vec)]  # remove elements with ".csv" extension
  max_model <- vec[which.max(as.numeric(str_sub(vec, -5,-1) ))]
  return(max_model)
}


for (num_CancerType in Cancerlist) {
  main.path_tc = paste0(filepath, "00.data/", num_CancerType, "/", folder_name)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  dual_cut_50 = list.files(paste0(main.path_tc), pattern = "model")
  best_model_folder = find_max_model(dual_cut_50)
  best_model = h2o.loadModel(paste0(main.path_tc,"/",best_model_folder,"/", gsub('.{6}$', '', best_model_folder)))
  best_model_vari_meta = as.data.frame(h2o.varimp(best_model))
  
  # links
  features_links = data.frame(matrix(ncol = 4))
  colnames(features_links) = c("from", "to","importance","rank")
  
  links_tmp = data.frame(matrix(ncol = 4))
  colnames(links_tmp) = c("from", "to","importance","rank")
  
  for (node_ids in best_model_vari_meta$variable) {
    if (length(unlist(str_split(node_ids, pattern = ""))) < 5) {
      from_pathway = paste0("P",str_split(node_ids,pattern = "P", simplify = T)[1,2])
      to_pathway = paste0("P",str_split(node_ids,pattern = "P", simplify = T)[1,2])
      links_tmp$from = from_pathway
      links_tmp$to = to_pathway
      features_links = rbind(features_links,links_tmp)
      remove(from_pathway,to_pathway)
      
    } else {
      from_pathway = paste0("P",str_split(node_ids,pattern = "P", simplify = T)[1,2])
      to_pathway = paste0("P",str_split(node_ids,pattern = "P", simplify = T)[1,3])
      links_tmp$from = from_pathway
      links_tmp$to = to_pathway
      features_links = rbind(features_links,links_tmp)
      remove(from_pathway,to_pathway)
    }
    
  } 
  
  features_links = features_links[-1,]
  rownames(features_links) = c(1:nrow(features_links))
  
  if (nrow(features_links) != length(best_model_vari_meta$variable)) {
    print("Something Wrong!!")
    break
  }
  
  # don't need
  features_links$rank = c(nrow(features_links):1)
  features_links$importance = best_model_vari_meta$relative_importance
  
  tmp_importance= preProcess(as.data.frame(features_links$importance), method=c("range"))
  
  features_links = cbind(features_links , as.data.frame(predict(tmp_importance, as.data.frame(features_links$importance))))
  features_links$min_max = features_links[,5]
  features_links$`features_links$importance` =NULL
  
  # top_features_each$`top_features_each$importance` =NULL
  # top_features_links$min_max = round(top_features_links$min_max,2)
  
  # nodes
  features_nodes = data.frame(matrix(ncol = 2, nrow = length(unique(c(features_links$from,features_links$to)))))
  colnames(features_nodes) = c("name", "origin")
  
  features_nodes$name = unique(c(features_links$from,features_links$to))
  
  for (node_num in 1:length(features_nodes$name)) {
    features_nodes[node_num,]$origin= pathway[which(pathway$No. == str_split(features_nodes$name,pattern = "P")[[node_num]][2]),]$`Kegg gene set` 
  }
  
  features_network <- graph_from_data_frame(d=features_links,
                                            vertices=features_nodes,
                                            directed=FALSE)
  # figure
  
  dir.create(paste0(main.path_tc,"/network"))
  fig_path = paste0(main.path_tc,"/network")
  setwd(fig_path)
  png(filename = paste0(CancerType, "_circle_network.png"),
      width = 1500, height = 1500, units = "px", pointsize = 12,
      bg = "white", res = NA, family = "")
  
  plot_circle = ggraph(features_network, layout = 'circle')+
    geom_edge_link(aes(width=importance,alpha=importance))+
    geom_node_point(aes(color=origin),size=7)+
    scale_size_continuous(range=c(2,7))+
    theme_void()+
    theme(legend.position='right')+
    geom_node_text(aes(label=name),
                   vjust=2.7)+
    scale_fill_discrete_sequential(palette='Red-Blue')+
    guides(col=guide_legend(ncol=1), fill=guide_legend(ncol=1))
  
  print(plot_circle)
  dev.off()
  
  png(filename = paste0(CancerType, "_circle_wo_network.png"),
      width = 1500, height = 1500, units = "px", pointsize = 12,
      bg = "white", res = NA, family = "")
  
  plot_circle = ggraph(features_network, layout = 'circle')+
    geom_edge_link(aes(color='gray50', alpha=.2))+
    geom_node_point(aes(color=origin),size=7)+
    scale_size_continuous(range=c(2,7))+
    theme_void()+
    theme(legend.position='right')+
    geom_node_text(aes(label=name),
                   vjust=2.7)+
    scale_fill_discrete_sequential(palette='Red-Blue')+
    guides(col=guide_legend(ncol=1), fill=guide_legend(ncol=1))
  
  print(plot_circle)
  dev.off()
  
  png(filename = paste0(CancerType, "_nicely_network.png"),
      width = 1500, height = 1500, units = "px", pointsize = 12,
      bg = "white", res = NA, family = "")
  
  plot_nicely = ggraph(features_network, layout = 'nicely')+
    geom_edge_link(aes(width=importance,alpha=importance))+
    geom_node_point(aes(color=origin),size=7)+
    scale_size_continuous(range=c(2,7))+
    theme_void()+
    theme(legend.position='right')+
    geom_node_text(aes(label=name),
                   vjust=2.7)+
    scale_fill_discrete_sequential(palette='Red-Blue')+
    guides(col=guide_legend(ncol=1), fill=guide_legend(ncol=1))
  
  print(plot_nicely)
  dev.off()
  
  png(filename =  paste0(CancerType, "_GC_network.png"),
      width = 1500, height = 1500, units = "px", pointsize = 12,
      bg = "white", res = NA, family = "")
  features_network_sim = igraph::simplify(features_network, remove.loops = TRUE, remove.multiple = FALSE)
  wc = cluster_edge_betweenness(features_network_sim)
  plot_wc = plot(wc, features_network_sim)
  print(plot_wc)
  
  dev.off()
  
  remove(plot_circle,plot_nicely, plot_wc, features_network_sim, features_network)
}
