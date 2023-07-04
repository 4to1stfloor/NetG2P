library(circlize)
library(tidyverse)
library(RColorBrewer)

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
fig_path = paste0(filepath,"04.Results/circos_cancerhallmark/")
if(!dir.exists(fig_path)){
  dir.create(fig_path)
  print(paste0("Created folder: ", fig_path))
} else {
  print(paste0("Folder already exists: ", fig_path))
}
setwd(fig_path)
Cancerlist

num_CancerType = "35.TCGA-KIDNEY"  
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
  
  circlize.tbl <- tmp_cancer_links_w_sl
  circlize.tbl <- circlize.tbl %>% mutate(sender = target, receiver = source) %>%
    dplyr::select(sender, receiver, everything()) %>% arrange(sender, desc(value)) %>%
    mutate(source_color = 'black')
  
  sender_order = circlize.tbl$sender %>% unique() # %>% sort()
  names(sender_order) = rep("sender",length(sender_order))
  receiver_order = circlize.tbl$receiver %>% unique() %>% sort() %>%
    lapply(function(receiver_oi) {
      circos_links_n = circlize.tbl %>% dplyr::filter(receiver == receiver_oi) %>% 
        group_by(receiver) %>% summarize(n = n()) %>% 
        ungroup()
      receivers = circlize.tbl %>% dplyr::filter(receiver == receiver_oi) %>% 
        inner_join(circos_links_n) %>% 
        dplyr::arrange(-n, sender) %>% 
        dplyr::distinct(receiver)
    }) %>% unlist()
  
  order = c(sender_order, receiver_order)
  
  sender_color <- set_names(brewer.pal(10, name='Paired'), sender_order %>% unique() %>% sort())
  receiver_color <- set_names(circlize.tbl$source_color, circlize.tbl$receiver)
  
  grid_col = c(sender_color, receiver_color)
  
  # circlize.tbl$ori_value = circlize.tbl$value
  
  circlize.tbl$value[circlize.tbl$value > quantile(circlize.tbl$value, probs = 0.75) + 1.5*IQR(circlize.tbl$value)] <- 
    round(quantile(circlize.tbl$value, probs = 0.75) + 1.5*IQR(circlize.tbl$value),1)
  
  transparency = circlize.tbl %>% mutate(value = (value - min(value))/(max(value)-min(value))) %>%
    mutate(transparency = 1-value) %>% .$transparency
  
  
  circos.clear()
  # circos.par(start.degree = 45,
  #            gap.degree=0.5,
  #            gap.after = c("CH9" = 1.5, "P9P54" = 1.5))
  
  circos.par(start.degree = -30,
             gap.degree= 0.5,
             gap.after = c(rep(1, length(sender_order)-1), 10, rep(1, length(receiver_order)-1), 10))
  
  # circos.par(start.degree = -55,
  #            gap.degree = 0.5)
  
  chordDiagram(circlize.tbl, directional = -1, order=order,
               grid.col = grid_col,
               transparency = transparency, 
               diffHeight = 0.05, 
               target.prop.height = mm_h(4),
               direction.type = c("diffHeight", "arrows"),
               link.arr.type = "big.arrow",
               # big.gap = 3,
               # link.arr.type = "curved",
               annotationTrack = "grid", 
               preAllocateTracks = list(track.height = 0.15))
  
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
                facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 1)
  }, bg.border = NA) #
  title(paste0(CancerType,"_critical_features_to_hallmark"))
  
  
}  
