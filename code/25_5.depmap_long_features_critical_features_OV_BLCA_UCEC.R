#making heatmaps
#Short/long analysis

library(openxlsx)
crit.sl.features = read.xlsx("~/nas/04.Results/critical_features/most_common_with_short_long_critical_features.xlsx")

remove.cancer = c("TCGA-COADREAD", "TCGA-LUSC")
crit.sl.features = crit.sl.features[,!colnames(crit.sl.features) %in% remove.cancer]


#making short long common balance
crit.sl.features$short.rat = NA
crit.sl.features$long.rat = NA
crit.sl.features$common.rat = NA

for (j in 1:nrow(crit.sl.features)) {
  query.cancer.type = crit.sl.features$which_cancer[j]
  query.cancer = sort(unlist(strsplit(query.cancer.type, ", ")))
  query.cancer.clean= query.cancer[!query.cancer %in% remove.cancer]
  
  crit.sl.features$count[j] = length(query.cancer.clean)
  crit.sl.features$which_cancer[j] = paste(query.cancer.clean, collapse = ", ")
  
  crit.sl.features$short.rat[j] = length(grep("short", crit.sl.features[j,])) / crit.sl.features$count[j]
  crit.sl.features$long.rat[j] = length(grep("long", crit.sl.features[j,])) / crit.sl.features$count[j]
  crit.sl.features$common.rat[j] = length(grep("common", crit.sl.features[j,])) / crit.sl.features$count[j]
}

crit.sl.features = crit.sl.features[order(crit.sl.features$count, decreasing = T),]
common.features = crit.sl.features[crit.sl.features$count > 1,]




#make a heatmap
################################################################################
common.feat.mat = common.features[,grep("TCGA", colnames(common.features))]
rownames(common.feat.mat) = common.features$features

common.feat.mat[common.feat.mat == "long"] = -1
common.feat.mat[common.feat.mat == "common"] = 0
common.feat.mat[common.feat.mat == "short"] = 1

common.feat.mat = as.matrix(common.feat.mat)
common.n.mat = apply(common.feat.mat, 2, as.numeric)
rownames(common.n.mat) = rownames(common.feat.mat)

common.n.raw = common.n.mat
#masking NAs
#dont do it too often
common.n.mat[is.na(common.n.mat)] = 0
#removing with more than 8 zeros
zero.test = rowSums(common.n.mat != 0)

common.nn.mat = common.n.mat[zero.test > 2,]

library(pheatmap)
pheatmap(common.n.mat, cluster_rows = T, cluster_cols = T)
pheatmap(common.nn.mat, cluster_rows = T, cluster_cols = T)


#
loss.count = rowSums(common.nn.mat[,colnames(common.nn.mat) %in% c("TCGA-LIHC", "TCGA-KIDNEY")])
loss.feat = names(loss.count)[abs(loss.count) > 1]

gain.count = rowSums(common.nn.mat[,colnames(common.nn.mat) %in% c("TCGA-STAD", "TCGA-OV", "TCGA-BLCA", "TCGA-UCEC")])
gain.feat = names(gain.count)[abs(gain.count) > 1]

common.filt = common.nn.mat[rownames(common.nn.mat) %in% unique(c(loss.feat, gain.feat)),]
test = pheatmap(common.filt, cluster_rows = T, cluster_cols = T)


cancer.order = c("TCGA-LIHC", "TCGA-KIDNEY", "TCGA-BRCA", "TCGA-LUAD", "TCGA-LGG", "TCGA-CESC", "TCGA-STAD", "TCGA-OV", "TCGA-BLCA", "TCGA-UCEC")
#extracting gene order
gene.order = test$tree_row$labels[test$tree_row$order]


common.n.mat = apply(common.feat.mat, 2, as.numeric)
rownames(common.n.mat) = rownames(common.feat.mat)

common.custom.mat = common.n.mat[match(gene.order, rownames(common.n.mat)), match(cancer.order, colnames(common.n.mat))]

common.cus.pheat = pheatmap(common.custom.mat, cluster_rows = F, cluster_cols = F)
################################################################################




#making 3long graph
feat.sel = rowSums(common.n.mat[, colnames(common.n.mat) %in% c("TCGA-OV", "TCGA-BLCA", "TCGA-UCEC")])
feat.sel = rowSums(common.n.mat[, colnames(common.n.mat) %in% c("TCGA-OV", "TCGA-BLCA", "TCGA-UCEC")])

feat.3l = names(feat.sel)[which(feat.sel == -3)]

#making edgelist
feat.edge = matrix(data = 0, nrow = length(feat.3l), ncol = 2)
for (i in 1:length(feat.3l)) {
  query = strsplit(feat.3l[i], "P")
  query.list = query[[1]]
  query.list = query.list[query.list != ""]
  query.list = paste("P", query.list, sep = "")
  
  feat.edge[i,1] = query.list[1]
  if (length(query.list) == 2) {
    feat.edge[i,2] = query.list[2]
  } else {
    feat.edge[i,2] = query.list[1]
  }
}


library(igraph)
graph.3l = graph_from_edgelist(feat.edge, directed = F)
graph.3ls = simplify(graph.3l)

library(ggraph)
edge.param = geom_edge_link(colour='#cfcbc4', edge_width=0.5)
label.gene.type = NULL
gene.type = NULL
vertex.degree = degree(graph.3ls)
deg.type = NULL


print.feed.cnv = ggraph(graph = graph.3ls)   +
  edge.param +
  label.gene.type +
  geom_node_point(aes(shape = gene.type, size = vertex.degree, color = deg.type)) +
  theme_void()

print(print.feed.cnv)


####
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

## Read Dataset 
meta.tbl <- read.csv(file = paste0(ref_path,'DepMap/Model.csv'), 
                     sep=',', header=T, fill=T)

ge.tbl <- read.csv(file = paste0(ref_path, 'DepMap/CRISPRGeneEffect.csv'), 
                   sep=',', header=T, fill=T, row.names=1)

# for fic
fig_path = "~/nas/04.Results/drug/depmap/"
if(!dir.exists(fig_path)){
  dir.create(fig_path)
  print(paste0("Created folder: ", fig_path))
} else {
  print(paste0("Folder already exists: ", fig_path))
}
setwd(fig_path)

transform_string <- function(s) {
  return(paste("GOBP_", gsub(" ", "_", toupper(s)), sep = ""))
}



has_self_loops <- any(is.loop(graph.3ls))

if (has_self_loops) {
  # Get the indices of the self-loops
  self_loop_indices <- which(is.loop(graph.3ls))
  
  # Remove the self-loops
  graph.3ls_re <- delete.edges(graph.3ls, self_loop_indices)
}

plot(graph.3ls_re)
str_split(E(graph.3ls_re)[1],"--")

incident(graph.3ls_re, "P18", mode = "all")

total_edges = data.frame()
for (max_degree in names(degree(graph.3ls_re)[degree(graph.3ls_re) == max(sort(degree(graph.3ls_re)))])) {
  tmp_edges <- as_edgelist(graph.3ls_re)[incident(graph.3ls_re, max_degree, mode = "all"),]
  total_edges = rbind(total_edges,tmp_edges)
}


# Extract numeric part from V1 and V2
numeric_part <- gsub("^P", "", total_edges$V1)
numeric_part2 <- gsub("^P", "", total_edges$V2)

total_edges$Combined <- ifelse(as.numeric(numeric_part) <= as.numeric(numeric_part2), 
                               paste(total_edges$V1, total_edges$V2, sep = ""),
                               paste(total_edges$V2, total_edges$V1, sep = ""))


depmap_common_genes <- unique(link_genes_filtered_df[which(link_genes_filtered_df$Pathway %in% total_edges$Combined),]$Genes)

# num_CancerType = "04.TCGA-CESC"
Cancerlist = c("10.TCGA-BLCA" , "24.TCGA-OV" , "32.TCGA-UCEC")
for (num_CancerType in Cancerlist) {
  
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  ## Select cells you want
  cell_oi <- meta.tbl %>% 
    #subset(GrowthPattern != 'Organoid' & DepmapModelType == 'COAD') %>%
    #subset(GrowthPattern != 'Organoid' & DepmapModelType == 'COAD' & PrimaryOrMetastasis == 'Primary') %>%
    # subset(GrowthPattern != 'Organoid' & DepmapModelType == 'BLCA' & PrimaryOrMetastasis == 'Primary' & SampleCollectionSite == 'Colon') %>%
    subset(GrowthPattern != 'Organoid' & 
             PrimaryOrMetastasis == 'Primary' &
             OncotreeLineage %in% unique(meta_for_TCGA[which(meta_for_TCGA$TCGA == CancerType),]$oncotree) &
             OncotreeCode %in% unique(meta_for_TCGA[which(meta_for_TCGA$TCGA %in% CancerType),]$oncotreecode)) %>%
    dplyr::select(ModelID) %>% unlist(use.names=F)
    
    
    ## Subset dataset
    ge.sub <- ge.tbl[rownames(ge.tbl) %in% cell_oi, c('cell_id', depmap_common_genes)] 
    ge.interest <- melt(ge.sub, id.vars = 'cell_id')
    
    gene.filter.plot <- as.data.frame(matrix(nrow = length(depmap_common_genes), ncol = 2))
    colnames(gene.filter.plot) <- c("Genes", "mean_depmap")
    
    gene.filter.plot$Genes <- depmap_common_genes
    # Mean(depmap score)
    gene.filter.plot$mean_depmap <- sapply(1:length(depmap_common_genes), 
                                           function(x) -mean(dplyr::filter(ge.interest, variable == depmap_common_genes[x])$value) ) 
    
    gene_filtere_edit = gene.filter.plot %>% filter(mean_depmap > 0.5)
    
    ge_edit_specific = ge.interest[which(ge.interest$variable %in% gene_filtere_edit$Genes),]
    
    png(filename = paste0(CancerType,"_depmap_specific_diff_long_cri_fe.png"),
        width = 35, height = 35,  units = "cm" ,pointsize = 12,
        bg = "white", res = 1200, family = "")
    
    ## Kenel Density Curve
    total_out = ggplot(ge_edit_specific, aes(x = value, y = variable, fill = after_stat(x))) +
      geom_density_ridges_gradient(scale = 2, gradient_lwd = .5, color = "black",
                                   jittered_points = TRUE,
                                   position = position_points_jitter(width = 0.05, height = 0),
                                   point_shape = '|', point_size = 3, point_alpha = 1, alpha = 0.7) +
      scale_fill_viridis_c(option = "plasma", name = "Gene Effect") +
      geom_vline(xintercept = -1, color='darkred')+
      geom_vline(xintercept = c(1,0,-2), color = 'darkgrey') +
      geom_vline(xintercept = c(-2.5, -1.5, -0.5, 0.5), color = 'darkgrey', linetype = 'dashed') +
      
      ylab('Gene') +
      xlab('Gene Effect \n(Chronos score, 23Q2)') +
      theme_bw()
    
    print(total_out)
    
    dev.off()
}
  

cell_oi <- meta.tbl %>% 
  #subset(GrowthPattern != 'Organoid' & DepmapModelType == 'COAD') %>%
  #subset(GrowthPattern != 'Organoid' & DepmapModelType == 'COAD' & PrimaryOrMetastasis == 'Primary') %>%
  # subset(GrowthPattern != 'Organoid' & DepmapModelType == 'BLCA' & PrimaryOrMetastasis == 'Primary' & SampleCollectionSite == 'Colon') %>%
  subset(GrowthPattern != 'Organoid' & 
           PrimaryOrMetastasis == 'Primary' &
           OncotreeLineage %in% unique(meta_for_TCGA[which(meta_for_TCGA$TCGA %in% c("TCGA-OV", "TCGA-BLCA", "TCGA-UCEC")),]$oncotree) &
           OncotreeCode %in% unique(meta_for_TCGA[which(meta_for_TCGA$TCGA %in% c("TCGA-OV", "TCGA-BLCA", "TCGA-UCEC")),]$oncotreecode)) %>%
  dplyr::select(ModelID) %>% unlist(use.names=F)

## Subset dataset
ge.sub <- ge.tbl[rownames(ge.tbl) %in% cell_oi, c('cell_id', depmap_common_genes)] 
ge.interest <- melt(ge.sub, id.vars = 'cell_id')

gene.filter.plot <- as.data.frame(matrix(nrow = length(depmap_common_genes), ncol = 2))
colnames(gene.filter.plot) <- c("Genes", "mean_depmap")

gene.filter.plot$Genes <- depmap_common_genes
# Mean(depmap score)
gene.filter.plot$mean_depmap <- sapply(1:length(depmap_common_genes), 
                                       function(x) -mean(dplyr::filter(ge.interest, variable == depmap_common_genes[x])$value) ) 

gene_filtere_edit = gene.filter.plot %>% filter(mean_depmap > 0.5)

ge_edit_specific = ge.interest[which(ge.interest$variable %in% gene_filtere_edit$Genes),]

png(filename = paste0("OV_UCEC_BLCA_depmap_specific_diff_long_cri_fe.png"),
    width = 35, height = 35,  units = "cm" ,pointsize = 12,
    bg = "white", res = 1200, family = "")

## Kenel Density Curve
total_out = ggplot(ge_edit_specific, aes(x = value, y = variable, fill = after_stat(x))) +
  geom_density_ridges_gradient(scale = 2, gradient_lwd = .5, color = "black",
                               jittered_points = TRUE,
                               position = position_points_jitter(width = 0.05, height = 0),
                               point_shape = '|', point_size = 3, point_alpha = 1, alpha = 0.7) +
  scale_fill_viridis_c(option = "plasma", name = "Gene Effect") +
  geom_vline(xintercept = -1, color='darkred')+
  geom_vline(xintercept = c(1,0,-2), color = 'darkgrey') +
  geom_vline(xintercept = c(-2.5, -1.5, -0.5, 0.5), color = 'darkgrey', linetype = 'dashed') +
  
  ylab('Gene') +
  xlab('Gene Effect \n(Chronos score, 23Q2)') +
  theme_bw()

print(total_out)

dev.off()
