library(ggraph)
library(igraph)
library(graphlayouts)
library(wesanderson)
library(tidyverse)
library(readxl)
library(graphlayouts)
library(colorspace)
library(caret)
library(tidygraph)
library(correlation)

filepath = "/home/seokwon/nas/"
ref_path = paste0(filepath, "99.reference/")
Cancerlist = dir(paste0(filepath, "/00.data/filtered_TCGA"))

surv_total_results = read_xlsx("~/nas/04.Results/Total_results_survpval2.xlsx")

# for fic
fig_path = "~/nas/04.Results/short_long/network/"
if(!dir.exists(fig_path)){
  dir.create(fig_path)
  print(paste0("Created folder: ", fig_path))
} else {
  print(paste0("Folder already exists: ", fig_path))
}
setwd(fig_path)
# pathway = as.data.frame(readxl::read_xlsx(paste0(ref_path ,"kegg_gene_set_w_cancer_hallmarks.xlsx")))
# num_CancerType = "04.TCGA-CESC"

for (num_CancerType in Cancerlist) {
  
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  Cancernum = gsub('[.-]','',gsub('[a-zA-Z]','', num_CancerType))
  # call input
  
  cancer_bf = read.csv(paste0(filepath,"04.Results/bestfeatures/",CancerType, "_best_features.csv"))
  short_long = readRDS(paste0(filepath,"04.Results/short_long/",CancerType,"_best_features_short_long.rds"))
  
  cancer_bf_cut = cancer_bf[1:surv_total_results[which(surv_total_results$CancerType == CancerType),]$num_of_features,]
  
  tmp_cancer_nodes = cancer_bf_cut$variable 
  
  short_group <- short_long %>%
    filter(cluster == "short") %>%
    select(-vitalstatus, -duration, -status, -cluster)
  
  long_group <- short_long %>%
    filter(cluster == "long") %>%
    select(-vitalstatus, -duration, -status, -cluster)
  
  merge_group <- short_long %>%
    select(-vitalstatus, -duration, -status, -cluster)
  
  short_correleation = as.data.frame(correlation::correlation(short_group,
                                                              include_factors = TRUE,
                                                              method = "auto")) %>% filter(p < 0.05)
  short_correleation$cluster = "short"
  long_correleation = as.data.frame(correlation::correlation(long_group,
                                                             include_factors = TRUE,
                                                             method = "auto")) %>% filter(p < 0.05)
  long_correleation$cluster = "long"
  

  # Step 1: Classify short pairs (present in short_correlation but not in long_correlation)
  short_pairs <- anti_join(short_correleation , long_correleation, by = c("Parameter1", "Parameter2"))
  
  # Step 2: Classify long pairs (present in long_correlation but not in short_correlation)
  long_pairs <- anti_join(long_correleation , short_correleation, by = c("Parameter1", "Parameter2"))
  
  # Step 3: Classify common pairs (present in both data frames)
  common_pairs <- inner_join(short_correleation , long_correleation, by = c("Parameter1", "Parameter2"))
  
  common_pairs_short <- common_pairs %>%
    group_by(Parameter1, Parameter2) %>%
    filter(r.x > r.y) %>%
    ungroup() %>%
    select(Parameter1, Parameter2, r = r.x, CI.x,CI_low.x,CI_high.x ,t.x, df_error.x, p.x, Method.x , n_Obs.x ,cluster.x)
  
  common_pairs_long <- common_pairs %>%
    group_by(Parameter1, Parameter2) %>%
    filter(r.y > r.x) %>%
    ungroup() %>%
    select(Parameter1, Parameter2, r = r.y, CI.y,CI_low.y,CI_high.y ,t.y, df_error.y, p.y, Method.y , n_Obs.y ,cluster.y)
  
  colnames(common_pairs_short) = colnames(short_pairs)
  colnames(common_pairs_long) = colnames(long_pairs)
  
  filtered_correlation = rbind(common_pairs_short,short_pairs,common_pairs_long,long_pairs)
  
  filtered_correlation <- arrange(filtered_correlation, desc(r))
  
  if (nrow(filtered_correlation) > 50) {
    filtered_correlation = filtered_correlation[1:50,]
  }

  #
  tmp_links = filtered_correlation %>% select(Parameter1, Parameter2 , r , cluster)
  colnames(tmp_links) = c("from" , "to" , "value", "cluster")
  
  tmp_short_links = short_correleation %>% select(Parameter1, Parameter2, r)
  colnames(tmp_short_links) = c("from" , "to" , "value")
  
  tmp_long_links = long_correleation %>% select(Parameter1, Parameter2, r)
  colnames(tmp_long_links) = c("from" , "to" , "value")
  
  cancer_sl_network <- graph_from_data_frame(d=tmp_links,
                                             vertices=tmp_cancer_nodes,
                                             directed=FALSE)
  
  cancer_short_network <- graph_from_data_frame(d=tmp_short_links,
                                             vertices=tmp_cancer_nodes,
                                             directed=FALSE)
  
  cancer_long_network <- graph_from_data_frame(d=tmp_long_links,
                                             vertices=tmp_cancer_nodes,
                                             directed=FALSE)
  
  # merge
  png(filename = paste0(CancerType, "_merge_network.png"),
      width = 35, height = 35,  units = "cm" ,pointsize = 12,
      bg = "white", res = 1200, family = "")
  
  cluster_colors <- c("short" = "red", "long" = "darkgreen")
  
  plort_merge = ggraph(cancer_sl_network, layout = 'circle')+
    geom_edge_link(aes(width=value,alpha=value,color = cluster))+
    geom_node_point(size=7)+
    scale_size_continuous(range=c(2,7))+
    theme_void()+
    theme(legend.position='right')+
    geom_node_text(aes(label=name),
                   vjust=2.7)+
    scale_edge_color_manual(values = cluster_colors) +
    guides(col=guide_legend(ncol=1), fill=guide_legend(ncol=1))
  
  print(plort_merge)
  
  dev.off()
  # short
  
  png(filename = paste0(CancerType, "_short_network.png"),
      width = 35, height = 35,  units = "cm" ,pointsize = 12,
      bg = "white", res = 1200, family = "")
  
  plot_circle_short = ggraph(cancer_short_network, layout = 'circle')+
    geom_edge_link(aes(width=value,alpha=value),color = "red")+
    geom_node_point(size=7)+
    scale_size_continuous(range=c(2,7))+
    theme_void()+
    theme(legend.position='right')+
    geom_node_text(aes(label=name),
                   vjust=2.7)+
    # scale_fill_discrete_sequential(palette='Red-Blue')+
    guides(col=guide_legend(ncol=1), fill=guide_legend(ncol=1))
  
  print(plot_circle_short)
  
  dev.off()
  
  # long
  
  png(filename = paste0(CancerType, "_long_network.png"),
      width = 35, height = 35,  units = "cm" ,pointsize = 12,
      bg = "white", res = 1200, family = "")

  plot_circle_long = ggraph(cancer_long_network, layout = 'circle')+
    geom_edge_link(aes(width=value,alpha=value ),color = "darkgreen")+
    geom_node_point(size=7)+
    scale_size_continuous(range=c(2,7))+
    theme_void()+
    theme(legend.position='right')+
    geom_node_text(aes(label=name),
                   vjust=2.7)+
    # scale_fill_discrete_sequential(palette='Red-Blue')+
    guides(col=guide_legend(ncol=1), fill=guide_legend(ncol=1))
  
  print(plot_circle_long)
  
  dev.off()
  
  remove(plot_circle_long,plot_circle_short, plort_merge, cancer_long_network, cancer_short_network, cancer_sl_network, filtered_correlation)
}
