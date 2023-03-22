# NETWORK upload and carculate centrality
library(igraph)
library(pheatmap)
library(survival) 
library(survminer) 

filepath = "/home/seokwon/nas/"
ref_path = paste0(filepath, "99.reference/")
Cancerlist = dir(paste0(filepath, "/00.data/filtered_TCGA/"))

folder_name = "h2o_bias_pval_dual_cut_50/network/"

swap_func <- function(row) {
  if (as.numeric(gsub("P", "", row[1])) > as.numeric(gsub("P", "", row[2]))) {
    return(c(row[2], row[1]))
  } else {
    return(row)
  }
}

for (num_CancerType in Cancerlist) {
  
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  network_path = paste0(main.path_tc, "/", folder_name)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
 
  cancer_network = readRDS(paste0(network_path,CancerType,"_best_features_network.rds"))

  features_network_sim = igraph::simplify(cancer_network, remove.loops = TRUE, remove.multiple = FALSE)
  wc = cluster_edge_betweenness(features_network_sim)
  
  # graph largest group
  sizes <- wc$membership %>% table()
  largest_cluster <- names(sizes)[which.max(sizes)]
  
  largest_cluster_vertices <- V(features_network_sim)[wc$membership == largest_cluster]
  
  # fig_folder_create
  if (!file.exists(paste0(main.path_tc,"/cluster_fig2"))) {
    dir.create(paste0(main.path_tc,"/cluster_fig2"))  
  } else {
    message("Folder already exists. Continuing...")
  }

  fig_path = paste0(main.path_tc,"/cluster_fig2")
  setwd(fig_path)
  
  png(filename = paste0(CancerType,"_largest_group.png"),
      width = 1000, height = 1000, units = "px", pointsize = 12,
      bg = "white", res = NA, family = "")
  
  large_group = plot(subgraph(features_network_sim, largest_cluster_vertices))
  
  print(large_group)
  dev.off()
  
  # get features from networks 
  group_count <- table(wc$membership)
  largest_group <- which.max(group_count)
  subgraph = induced_subgraph(features_network_sim, which(wc$membership == largest_group))
  
  edge_list <- get.edgelist(subgraph)
  edge_list = as.data.frame(edge_list)
  edge_list_swapped <- as.data.frame(t(apply(edge_list, 1, swap_func)))
  
  # convert the result back to a matrix
  edge_list_swapped <- matrix(edge_list_swapped, ncol = 2, byrow = TRUE)
  reverse_edges <- data.frame( from=edge_list_swapped[,2], to=edge_list_swapped[,1])
  combined_vertices <- cbind(reverse_edges[,2], reverse_edges[,1])
  combined_vertices <- data.frame(combined_vertices, combined = paste0(combined_vertices[,1], combined_vertices[,2], sep = ""))

  # call input
  annotate_best_features =  read.csv(paste0(filepath,"13.analysis/",CancerType,"_best_features.csv"))
  duration_log_df = readRDS(paste0(main.path_tc, "/", CancerType,"_dual_add_duration_log.rds"))
  
  # results_pval
  result_surv_pval = data.frame(matrix(ncol = 1))
  colnames(result_surv_pval) = "pval"
  
  best_features_df = duration_log_df[,combined_vertices$combined]
  
  best_features_df$vitalstatus = duration_log_df$vitalstatus
  best_features_df$duration = duration_log_df$duration
  
  best_features_df$status = NA
  best_features_df$status[which(best_features_df$vitalstatus == "Dead")] = 1
  best_features_df$status[which(best_features_df$vitalstatus == "Alive")] = 0
 
  for (last_num in 2:length(combined_vertices$combined)) {
    
    out = pheatmap((best_features_df[,1:last_num] > -log(0.05))*1 , cluster_cols = T,
                   cluster_rows = T, labels_cols = "",
                   show_rownames = T)
    print(paste0(CancerType , "_start : ",last_num ))
    print("----------------------------------------")
    # print(table(cutree(out$tree_row, 2)))
    
    tmp_pheat_cut = as.data.frame (cutree(out$tree_row, 2) , out[["tree_row"]][["labels"]])
    colnames(tmp_pheat_cut) = "cluster"
    best_features_df$cluster = tmp_pheat_cut$cluster
    fit = survfit(Surv(duration, status) ~ cluster, data = best_features_df)
    
    png(filename = paste0("cut",last_num ,"_cluster_test.png"),
        width = 1000, height = 1000, units = "px", pointsize = 12,
        bg = "white", res = NA, family = "")
    
    plot_roc = ggsurvplot(fit, data = best_features_df, risk.table = TRUE,
                          palette = "jco", pval = TRUE, surv.median.line = "hv", xlab = "days")
    print(plot_roc)
    dev.off()
    
    
    result_surv_pval[last_num,"pval"] = surv_pvalue(fit)$pval
    remove(tmp_pheat_cut,fit,plot_roc,combined_vertices)
  }
  remove(best_features_df,annotate_best_features,duration_log_df )
  write.csv(result_surv_pval, paste0(CancerType, "_result_survpval.csv"))
  
}
