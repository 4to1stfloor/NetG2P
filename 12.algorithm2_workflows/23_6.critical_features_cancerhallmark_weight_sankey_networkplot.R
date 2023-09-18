## need critical features and cancer hallmark
## sankeyNetwork
library(networkD3)
library(magrittr)
library(stringr)
library(data.table)
library(pheatmap)
library(readxl)
library(dplyr)
library(readxl)
library(tidyverse)

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

surv_total_results = read_xlsx("~/nas/04.Results/Total_results_survpval2.xlsx")
cancerhallmark = read_xlsx(paste0(ref_path, "kegg_gene_set_w_cancer_hallmarks_edit.xlsx"))  

# for fic
fig_path = paste0(filepath,"04.Results/sankey_cancerhallmark/short_long")
if(!dir.exists(fig_path)){
  dir.create(fig_path)
  print(paste0("Created folder: ", fig_path))
} else {
  print(paste0("Folder already exists: ", fig_path))
}
setwd(fig_path)

# for all
for (num_CancerType in Cancerlist) {
  
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  
  # call input
  
  # t.test 한거
  # short_long_features = read_xlsx(paste0(filepath, "04.Results/short_long/",CancerType, "_best_features_short_long.xlsx"))
  bf_short_long = readRDS(paste0(filepath, "04.Results/short_long/",CancerType,"_best_features_short_long.rds"))
  best_features_importance = readRDS(paste0(main.path_tc, "/h2o_bias_pval_dual_cut_50/network/",CancerType,"_best_features_links.rds"))
  cancer_bf = read.csv(paste0(filepath,"04.Results/bestfeatures/",CancerType, "_best_features.csv"))

  # 일단 short long 합쳐서 진행 -> short long ttest전 features
  # cut 100
  # if (surv_total_results[which(surv_total_results$CancerType == CancerType),]$num_of_features > 100) {
  #   cancer_bf_cut = cancer_bf[1:100,]
  #   cancer_bf_cut$minmax = best_features_importance[1:100,]$min_max
  # }else {
  #   cancer_bf_cut = cancer_bf[1:surv_total_results[which(surv_total_results$CancerType == CancerType),]$num_of_features,]
  #   cancer_bf_cut$minmax = best_features_importance[1:surv_total_results[which(surv_total_results$CancerType == CancerType),]$num_of_features,]$min_max
  # }
  
  cancer_bf_cut = cancer_bf[1:surv_total_results[which(surv_total_results$CancerType == CancerType),]$num_of_features,]
  cancer_bf_cut$minmax = best_features_importance[1:surv_total_results[which(surv_total_results$CancerType == CancerType),]$num_of_features,]$min_max
  
  cancer_logp = bf_short_long[,cancer_bf_cut$variable]
  
  # # 일단 unique gene으로 해봄 
  # total_link_genes = link_genes_filtered_df[which(link_genes_filtered_df$Pathway %in% cancer_bf_cut$variable),]$Genes
  # total_single_genes = single_genes[which(single_genes$Pathway %in% cancer_bf_cut$variable),]$Genes
  # total_bf_genes = c(total_link_genes,total_single_genes)
  # total_uni_bf_genes = unique(total_bf_genes)
  
  split_elements <- strsplit(cancer_bf_cut$variable, "P")
  divide_features = c()
  # Extract the first and second parts if the element has two "P" 
  for (i in 1:length(split_elements)) {
    if (length(split_elements[[i]]) > 2) {
      divide_features = c(divide_features, paste0("P", split_elements[[i]][2]))
      divide_features = c(divide_features, paste0("P", split_elements[[i]][3]))
    } else {
      divide_features = c(divide_features, paste0("P", split_elements[[i]][2]))
    }
  }
  
  divide_features = unique(divide_features)
  cancerhallmark_filt = as.data.frame(cancerhallmark[which(cancerhallmark$pathway %in% divide_features),])
  
  tmp_cancer_nodes = as.data.frame(c(cancer_bf_cut$variable , unique(cancerhallmark_filt$no_cancer_hallmark)))
  # tmp_cancer_nodes = as.data.frame(c(total_uni_bf_genes,cancer_bf_cut$variable , unique(cancerhallmark_filt$no_cancer_hallmark)))
  colnames(tmp_cancer_nodes) = "name"
  
  tmp_cancer_links = data.frame(matrix(ncol = 3))
  colnames(tmp_cancer_links) = c("source", "target", "value")
  
  for (tmp_cancer_feature in colnames(cancer_logp)) {
    
    tmp_features <- strsplit(tmp_cancer_feature, "P")
    if (length(tmp_features[[1]]) > 2) {
      
      first_pathwaylink = data.frame(source = tmp_cancer_feature , 
                                     target = cancerhallmark_filt[which(cancerhallmark_filt$pathway == paste0("P", tmp_features[[1]][2])),]$no_cancer_hallmark,
                                     value = as.numeric(colSums(cancer_logp)[tmp_cancer_feature]) * 
                                       cancer_bf_cut[which(cancer_bf_cut$variable == tmp_cancer_feature),]$minmax *
                                       cancerhallmark_filt[which(cancerhallmark_filt$pathway == paste0("P", tmp_features[[1]][2])),]$weight)
      
      second_pathwaylink = data.frame(source =  tmp_cancer_feature , 
                                      target = cancerhallmark_filt[which(cancerhallmark_filt$pathway == paste0("P", tmp_features[[1]][3])),]$no_cancer_hallmark,
                                      value = as.numeric(colSums(cancer_logp)[tmp_cancer_feature]) *
                                        cancer_bf_cut[which(cancer_bf_cut$variable == tmp_cancer_feature),]$minmax *
                                        cancerhallmark_filt[which(cancerhallmark_filt$pathway == paste0("P", tmp_features[[1]][3])),]$weight)
      
      tmp_edge_df = rbind(first_pathwaylink,second_pathwaylink)
      
      # duplication mulitiple method
      
      if (length(unique(tmp_edge_df$target)) != length(tmp_edge_df$target) ) {
        dup_target = tmp_edge_df[which(duplicated(tmp_edge_df$target)),]
        
        for (dup_tar in dup_target$target) {
          dup_filt = tmp_edge_df[which(tmp_edge_df$target == dup_tar),]
          # multiple_num = length(dup_filt$target)
          dup_filt = dup_filt[1,]
          # dup_filt$value = dup_filt$value * multiple_num
          tmp_edge_df = tmp_edge_df[which(!tmp_edge_df$target %in% dup_tar),]
          tmp_edge_df = rbind(tmp_edge_df, dup_filt)
        }
      }
      
    } else {
      
      tmp_edge_df = data.frame(source = tmp_cancer_feature , 
                               target = cancerhallmark_filt[which(cancerhallmark_filt$pathway == paste0("P", tmp_features[[1]][2])),]$no_cancer_hallmark,
                               value = as.numeric(colSums(cancer_logp)[tmp_cancer_feature]) * 
                                 cancer_bf_cut[which(cancer_bf_cut$variable == tmp_cancer_feature),]$minmax *
                                 cancerhallmark_filt[which(cancerhallmark_filt$pathway == paste0("P", tmp_features[[1]][2])),]$weight)
    }
    
    tmp_cancer_links = rbind(tmp_cancer_links, tmp_edge_df)
    
    remove(tmp_edge_df ,first_pathwaylink,second_pathwaylink)
    
  }
  
  tmp_cancer_links = tmp_cancer_links[-1,]
  rownames(tmp_cancer_links) = NULL
  
  total_features_short_long = data.frame()
  tmp_judge = data.frame(matrix(ncol = 2))
  colnames(tmp_judge) =  c("source","cluster")
  for (critical_feature in cancer_bf_cut$variable) {
    if (sum(bf_short_long[which(bf_short_long$cluster == "long"),][,critical_feature]) >
        sum(bf_short_long[which(bf_short_long$cluster == "short"),][,critical_feature])) {
      tmp_judge$source = critical_feature
      tmp_judge$cluster = "long"
      
    }else {
      tmp_judge$source = critical_feature
      tmp_judge$cluster = "short"
    }
    
    total_features_short_long = rbind(total_features_short_long, tmp_judge)
  }
 
  tmp_cancer_links_w_sl = merge(tmp_cancer_links,total_features_short_long)

  # transfer num of node
  for (tmp_source in unique(tmp_cancer_links_w_sl$source )) {
    source_num = as.numeric(rownames(tmp_cancer_nodes)[which(tmp_cancer_nodes$name == tmp_source)])
    tmp_cancer_links_w_sl[which(tmp_cancer_links_w_sl$source == tmp_source),]$source = source_num
    
  }
  
  for (tmp_target in unique(tmp_cancer_links_w_sl$target )) {
    target_num = as.numeric(rownames(tmp_cancer_nodes)[which(tmp_cancer_nodes$name == tmp_target)])
    tmp_cancer_links_w_sl[which(tmp_cancer_links_w_sl$target == tmp_target),]$target = target_num
    
  }
  
  tmp_cancer_links_w_sl$source = as.numeric(tmp_cancer_links_w_sl$source) -1
  tmp_cancer_links_w_sl$target = as.numeric(tmp_cancer_links_w_sl$target) -1
  
 
  
  sankeyNetwork_li = list(nodes = tmp_cancer_nodes , links = tmp_cancer_links_w_sl)
  
  
  
  sankeyNetwork(Links = sankeyNetwork_li$links, Nodes = sankeyNetwork_li$nodes, Source = "source",
                Target = "target", Value = "value", NodeID = "name",LinkGroup = "cluster",
                fontSize = 12) %>%
    saveNetwork(file = paste0(CancerType,"_sanketnet_ch_bf_short_long.html"), FALSE)
  
}  




# 
# 
# # Read data
# URL <- paste0("https://cdn.rawgit.com/christophergandrud/networkD3/","master/JSONdata/energy.json")
# 
# # Convert to list format
# Energy <- jsonlite::fromJSON("https://raw.githubusercontent.com/apache/incubator-echarts/master/test/data/energy.json")
# Energy <- jsonlite::fromJSON(URL)
# 
# class(Energy)
# # Create Sankey Chart and save to file          
# sankeyNetwork(Links = Energy$links, Nodes = Energy$nodes, Source = "source", Target = "target", Value = "value", 
#               NodeID = "name", units = "TWh", fontSize = 12, nodeWidth = 30) %>%
#   saveNetwork(file = "Net2.html", FALSE)
# 
# 
# 
# sankeyNetwork(Links = Energy$links, Nodes = Energy$nodes, Source = "source",
#               Target = "target", Value = "value", NodeID = "name",
#               units = "TWh", fontSize = 12, nodeWidth = 30)
