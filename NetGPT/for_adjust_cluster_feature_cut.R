phyper = readRDS(paste0(main.path_tc, "/phyper_enrichment_54_KEGG_net_GC_20220117.rds"))
filt_phyper = readRDS(file = paste0(main.path_tc,"/",CancerType,"_pathwaylink_all_0or1.rds"))

top_5_with_adjust_cluster = readRDS(file = paste0(main.path_tc ,CancerType , "_",top_num , "_matrix_origin_with_cluster.rds")) 
top_5_with_adjust_cluster
origin_for_matrix = as.data.frame(filt_phyper)
origin_for_matrix_01 = origin_for_matrix[,-length(origin_for_matrix)]

origin_for_matrix_01_filt = origin_for_matrix_01
origin_for_matrix_01_filt$cluster = NA
for (top_patient in rownames(top_5_with_adjust_cluster)) {
  origin_for_matrix_01_filt[which(rownames(origin_for_matrix_01_filt) == top_patient),]$cluster = top_5_with_adjust_cluster[top_patient,]$cluster
}

origin_matrix_01_filt_wo_cluster = origin_for_matrix_01_filt[,-ncol(origin_for_matrix_01_filt)]

features = rep(colnames(origin_matrix_01_filt_wo_cluster),each = 4)
clusters =  rep(c("cluster1","cluster1","cluster2","cluster2"), time = 644)
exist = rep(c("yes","no"), time = 644 *2)
feature_1_count_adjust_cluster = data.frame(features,clusters,exist)

feature_1_count_adjust_cluster$counts = NA

cluster_tmp1 = as.data.frame(colSums(origin_for_matrix_01_filt[which(origin_for_matrix_01_filt$cluster == 1),])[-645])
colnames(cluster_tmp1) = "counts"
cluster_tmp2 = as.data.frame(colSums(origin_for_matrix_01_filt[which(origin_for_matrix_01_filt$cluster == 2),])[-645])
colnames(cluster_tmp2) = "counts"

feature_1_count_adjust_cluster[which(feature_1_count_adjust_cluster$clusters == "cluster1" & feature_1_count_adjust_cluster$exist == "yes"),]$counts = cluster_tmp1$counts
feature_1_count_adjust_cluster[which(feature_1_count_adjust_cluster$clusters == "cluster2" & feature_1_count_adjust_cluster$exist == "yes"),]$counts = cluster_tmp2$counts

feature_1_count_adjust_cluster[which(feature_1_count_adjust_cluster$clusters == "cluster1" & feature_1_count_adjust_cluster$exist == "no"),]$counts = 
  nrow(origin_for_matrix_01_filt[which(origin_for_matrix_01_filt$cluster == 1),]) - cluster_tmp1$counts

feature_1_count_adjust_cluster[which(feature_1_count_adjust_cluster$clusters == "cluster2" & feature_1_count_adjust_cluster$exist == "no"),]$counts = 
  nrow(origin_for_matrix_01_filt[which(origin_for_matrix_01_filt$cluster == 2),]) - cluster_tmp2$counts

feature_1_count_adjust_cluster

# fisher test

feature_pval = data.frame(matrix(nrow = length(unique(feature_1_count_adjust_cluster$features))))
rownames(feature_pval) = unique(feature_1_count_adjust_cluster$features)
colnames(feature_pval) = "pval"
if (sum(is.na(feature_1_count_adjust_cluster$count)) == 0) {
  for (features_tmp in unique(feature_1_count_adjust_cluster$features)) {
    tmp_tab = feature_1_count_adjust_cluster[which(feature_1_count_adjust_cluster$features == features_tmp),][,2:4]
    
    result_fisher = fisher.test(xtabs(counts~clusters+exist,data = tmp_tab), workspace = 2e8)
    feature_pval[features_tmp,] =result_fisher$p.value
    remove(tmp_tab)
  }
  
} else {
  print("It has NA")
  }

feature_pval_specific = subset(feature_pval, pval < 0.05)
dl_vari_meta = h2o.varimp(best_dl)
dl_vari_feature = dl_vari_meta[which(dl_vari_meta$variable %in% rownames(feature_pval_specific)),][,1:2]
feature_pval_specific$variable = rownames(feature_pval_specific)
feature_pval_specific_vari = merge(feature_pval_specific,dl_vari_feature, by = "variable")

feature_pval_specific_vari_order = data.frame()
for (order_id in rownames(feature_pval_specific)) {
  feature_pval_specific_vari_order = rbind(feature_pval_specific_vari_order,feature_pval_specific_vari[which(feature_pval_specific_vari$variable == order_id),])
}

rownames(feature_pval_specific_vari_order) = feature_pval_specific_vari_order$variable
feature_pval_specific_vari_order = feature_pval_specific_vari_order[,2:3]
feature_pval_specific_vari_order$negative_log10 = -log10(feature_pval_specific_vari_order$pval) 
feature_pval_specific_vari_order$negative_loge = -log(feature_pval_specific_vari_order$pval) 

ggplot(data=feature_pval_specific_vari_order, aes(x=relative_importance, y=negative_log10, label = rownames(feature_pval_specific_vari_order))) + 
     geom_point(shape=19, size=3, colour="red") + # shape 19: solid circle
     geom_text_repel()+ theme_classic() 

ggplot(data=feature_pval_specific_vari_order, aes(x=relative_importance, y=negative_loge, label = rownames(feature_pval_specific_vari_order))) + 
  geom_point(shape=19, size=3, colour="red") + # shape 19: solid circle
  geom_text_repel()+ theme_classic() 
