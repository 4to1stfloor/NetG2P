library(stringr)
library(data.table)
library(pheatmap)
library(readxl)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(dplyr)
library(readxl)
library(tidyverse)
library(rrvgo)
library(treemapify)

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
sce_path = "/mnt/gluster_server/data/raw/TCGA_data/00.data/"

surv_total_results = read_xlsx("~/nas/04.Results/Total_results_survpval2.xlsx")

# for fic
fig_path = paste0(filepath,"04.Results/GOenrichment_test/treemap_edit/")
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
  
  short_long_features = read_xlsx(paste0(filepath, "04.Results/short_long/",CancerType, "_best_features_short_long.xlsx"))
  # bf_short_long = readRDS(paste0(filepath, "04.Results/short_long/",CancerType,"_best_features_short_long.rds"))
  # duration_log_df = readRDS(paste0(main.path_tc, "/", CancerType,"_dual_add_duration_log.rds"))
  
  # 일단 unique gene으로 해봄 
  # total_link_genes = link_genes_filtered_df[which(link_genes_filtered_df$Pathway %in% short_long_features$variable),]$Genes
  # total_single_genes = single_genes[which(single_genes$Pathway %in% short_long_features$variable),]$Genes
  # total_bf_genes = c(total_link_genes,total_single_genes)
  
  short_features = short_long_features[which(short_long_features$classification == "short"),]$variable 
  long_features = short_long_features[which(short_long_features$classification == "long"),]$variable 
  
  short_gene = c()
  long_gene = c()
  
  # short
  for (sf in short_features) {
    count <- str_count(sf, "P")
    if (count == 1) {
      tmp_short_gene = single_genes[which(single_genes$Pathway == sf),]$Genes
    } else {
      tmp_short_gene = link_genes_filtered_df[which(link_genes_filtered_df$Pathway == sf),]$Genes
    }
    short_gene = c(short_gene, tmp_short_gene)
    
  }
  
  # long
  for (lf in long_features) {
    count <- str_count(lf, "P")
    if (count == 1) {
      tmp_long_gene = single_genes[which(single_genes$Pathway == lf),]$Genes
    } else {
      tmp_long_gene = link_genes_filtered_df[which(link_genes_filtered_df$Pathway == lf),]$Genes
    }
    long_gene = c(long_gene, tmp_long_gene)
    
  }
  
  if (length(long_gene) != 0 && length(short_gene) != 0 ) {
    short_gene_en = AnnotationDbi::select(org.Hs.eg.db, short_gene, 'ENTREZID', 'SYMBOL')[
      which(!is.na(AnnotationDbi::select(org.Hs.eg.db, short_gene, 'ENTREZID', 'SYMBOL')$ENTREZID)),]
    short_gene_en <- data.frame(short_gene_en, row.names = NULL)
    
    long_gene_en = AnnotationDbi::select(org.Hs.eg.db, long_gene, 'ENTREZID', 'SYMBOL')[
      which(!is.na(AnnotationDbi::select(org.Hs.eg.db, long_gene, 'ENTREZID', 'SYMBOL')$ENTREZID)),]
    long_gene_en <- data.frame(long_gene_en, row.names = NULL)
    
    short_enrichGO = enrichGO(short_gene_en$ENTREZID ,
                              OrgDb = org.Hs.eg.db ,
                              keyType = "ENTREZID" , 
                              ont = "BP" ,
                              pvalueCutoff  = 0.05,
                              qvalueCutoff = 0.01)
    
    long_enrichGO = enrichGO(long_gene_en$ENTREZID ,
                             OrgDb = org.Hs.eg.db ,
                             keyType = "ENTREZID" , 
                             ont = "BP" ,
                             pvalueCutoff  = 0.05,
                             qvalueCutoff = 0.01 )
    
    # qval cut -> 
    short_enrichGO_df = as.data.frame(short_enrichGO)
    long_enrichGO_df = as.data.frame(long_enrichGO)
    
    simMatrix_short <- calculateSimMatrix(short_enrichGO_df$ID,
                                          orgdb=org.Hs.eg.db,
                                          ont="BP",
                                          method="Rel")
    
    simMatrix_long <- calculateSimMatrix(long_enrichGO_df$ID,
                                         orgdb=org.Hs.eg.db,
                                         ont="BP",
                                         method="Rel")
    
    scores_short <- setNames(-log10(short_enrichGO_df$qvalue), short_enrichGO_df$ID)
    scores_long <- setNames(-log10(long_enrichGO_df$qvalue), long_enrichGO_df$ID)
    
    
    
    if (sum(is.na(scores_short)) != 0) {
      scores_short[which(is.na(scores_short))] = 0
    }
    
    if (sum(is.na(scores_long)) != 0) {
      scores_long[which(is.na(scores_long))] = 0
    }
    
    if (sum(scores_short) == 0) {
      reducedTerms_short <- reduceSimMatrix(simMatrix_short,
                                            threshold=0.7,
                                            orgdb="org.Hs.eg.db")
    } else {
      reducedTerms_short <- reduceSimMatrix(simMatrix_short,
                                            threshold=0.7,
                                            scores_short,
                                            orgdb="org.Hs.eg.db")
    }
    
    if (sum(scores_long) == 0) {
      reducedTerms_long <- reduceSimMatrix(simMatrix_long,
                                           threshold=0.7,
                                           orgdb="org.Hs.eg.db")
    } else {
      reducedTerms_long <- reduceSimMatrix(simMatrix_long,
                                           scores_long,
                                           threshold=0.7,
                                           orgdb="org.Hs.eg.db")
    }
    
    score_short_head = head(unique(reducedTerms_short[order(reducedTerms_short$score,decreasing = T),]$parentTerm), n= 3)
    score_long_head = head(unique(reducedTerms_long[order(reducedTerms_long$score,decreasing = T),]$parentTerm), n= 3)
    
    count_short_head = names(tail(sort(table(reducedTerms_short$parentTerm)), n= 3))
    count_long_head = names(tail(sort(table(reducedTerms_long$parentTerm)), n= 3))
    
    total_short_imp = unique(c(total_short_imp, score_short_head, count_short_head))
    total_long_imp = unique(c(total_long_imp, score_long_head, count_long_head))
    
    reducedTerms_short_edit = reducedTerms_short %>% 
      mutate(color = case_when( parentTerm %in% total_short_imp ~ "#CC0000",
                                .default = "#999999")) %>% 
      mutate(subgroup_font_size = case_when( parentTerm %in% total_short_imp ~ 30, .default = 10))%>% 
      mutate(subgroup_font_fontface = case_when( parentTerm %in% total_short_imp ~ "bold",.default = "plain"))
    
    reducedTerms_long_edit = reducedTerms_long %>% 
      mutate(color = case_when( parentTerm %in% total_long_imp ~ "#339900",
                                .default = "#999999")) %>% 
      mutate(subgroup_font_size = case_when( parentTerm %in% total_long_imp ~ 30, .default = 10))%>% 
      mutate(subgroup_font_fontface = case_when( parentTerm %in% total_long_imp ~ "bold",.default = "plain"))
    
    # fiq

    png(filename = paste0(CancerType,"_tree_short_0.7_wo_ttest_not_exp_BP.png"),
        width = 35, height = 35,  units = "cm" ,pointsize = 12,
        bg = "white", res = 1200, family = "")
    
    tree_mat_short =  ggplot(reducedTerms_short_edit, aes(area = score, fill = color, 
                                                          label = term, # country
                                                          subgroup = parentTerm)) + # continent
      scale_fill_identity() +
      # 1. Draw country borders and fill colors
      geom_treemap(colour = "#FFFFFFDD", start = "topleft") +
      # 2. Draw continent borders
      geom_treemap_subgroup_border(
        # colour = "#00000080", 
        start = "topleft", 
        colour = "gray100",
        size = 5) +
      # 3. Print continent text
      
      geom_treemap_subgroup_text(place = "centre", colour = "black",
                                 reflow = T,
                                 fontface = reducedTerms_short_edit$subgroup_font_fontface,
                                 # min.size = 0
                                 start = "topleft",
                                 grow = F,
                                 size= reducedTerms_short_edit$subgroup_font_size) +
      # 4. Print country text
      geom_treemap_text(colour = "#FFFFFFDD", 
                        place = "centre",
                        alpha = 0.5,
                        grow = F,
                        reflow = T,
                        start = "topleft") +
      # option
      theme(legend.position = 0)
    
    print(tree_mat_short)
    dev.off()
    
    png(filename = paste0(CancerType,"_tree_long_0.7_wo_ttest_not_exp_BP.png"),
        width = 35, height = 35,  units = "cm" ,pointsize = 12,
        bg = "white", res = 1200, family = "")
    
    tree_mat_long =  ggplot(reducedTerms_long_edit, aes(area = score, fill = color, 
                                                        label = term, # country
                                                        subgroup = parentTerm)) + # continent
      scale_fill_identity() +
      # 1. Draw country borders and fill colors
      geom_treemap(colour = "#FFFFFFDD", start = "topleft") +
      # 2. Draw continent borders
      geom_treemap_subgroup_border(
        # colour = "#00000080", 
        start = "topleft", 
        colour = "gray100",
        size = 5) +
      # 3. Print continent text
      
      geom_treemap_subgroup_text(place = "centre", colour = "black",
                                 reflow = T,
                                 fontface = reducedTerms_long_edit$subgroup_font_fontface,
                                 # min.size = 0
                                 start = "topleft",
                                 grow = F,
                                 size= reducedTerms_long_edit$subgroup_font_size) +
      # 4. Print country text
      geom_treemap_text(colour = "#FFFFFFDD", 
                        place = "centre",
                        alpha = 0.5,
                        grow = F,
                        reflow = T,
                        start = "topleft") +
      # option
      theme(legend.position = 0)
    
    print(tree_mat_long)
    dev.off()
    
    total_gene_li = list(short =short_gene_en$ENTREZID , long = long_gene_en$ENTREZID)
    
    total_enrich <- compareCluster(geneCluster   = total_gene_li, 
                                   ont           = "BP",
                                   OrgDb         = org.Hs.eg.db,
                                   readable      = TRUE,
                                   pvalueCutoff  = 0.05,
                                   fun = "enrichGO")
    
    png(filename = paste0(CancerType,"_dotplot_wo_ttest_not_exp_BP.png"),
        width = 35, height = 35,  units = "cm" ,pointsize = 12,
        bg = "white", res = 1200, family = "")
    
    total_dotplot = dotplot(total_enrich)
    print(total_dotplot)
    dev.off()
    
  } else if (length(long_gene) == 0 ) {
    short_gene_en = AnnotationDbi::select(org.Hs.eg.db, short_gene, 'ENTREZID', 'SYMBOL')[
      which(!is.na(AnnotationDbi::select(org.Hs.eg.db, short_gene, 'ENTREZID', 'SYMBOL')$ENTREZID)),]
    short_gene_en <- data.frame(short_gene_en, row.names = NULL)
    
    short_enrichGO = enrichGO(short_gene_en$ENTREZID ,
                              OrgDb = org.Hs.eg.db ,
                              keyType = "ENTREZID" , 
                              ont = "BP" ,
                              pvalueCutoff  = 0.05,
                              qvalueCutoff = 0.01)
    
    short_enrichGO_df = as.data.frame(short_enrichGO)
    
    simMatrix_short <- calculateSimMatrix(short_enrichGO_df$ID,
                                          orgdb=org.Hs.eg.db,
                                          ont="BP",
                                          method="Rel")
    
    scores_short <- setNames(-log10(short_enrichGO_df$qvalue), short_enrichGO_df$ID)
    
    if (sum(is.na(scores_short)) != 0) {
      scores_short[which(is.na(scores_short))] = 0
    }
    
    if (sum(scores_short) == 0) {
      reducedTerms_short <- reduceSimMatrix(simMatrix_short,
                                            threshold=0.7,
                                            orgdb="org.Hs.eg.db")
    } else {
      reducedTerms_short <- reduceSimMatrix(simMatrix_short,
                                            scores_short,
                                            threshold=0.7,
                                            orgdb="org.Hs.eg.db")
    }
    
    score_short_head = head(unique(reducedTerms_short[order(reducedTerms_short$score,decreasing = T),]$parentTerm), n= 3)
    
    count_short_head = names(tail(sort(table(reducedTerms_short$parentTerm)), n= 3))
    
    total_short_imp = unique(c(total_short_imp, score_short_head, count_short_head))
    
    reducedTerms_short_edit = reducedTerms_short %>% 
      mutate(color = case_when( parentTerm %in% total_short_imp ~ "#CC0000",
                                .default = "#999999")) %>% 
      mutate(subgroup_font_size = case_when( parentTerm %in% total_short_imp ~ 30, .default = 10))%>% 
      mutate(subgroup_font_fontface = case_when( parentTerm %in% total_short_imp ~ "bold",.default = "plain"))
    
    # fiq
    
    png(filename = paste0(CancerType,"_tree_short_0.7_wo_ttest_not_exp_BP.png"),
        width = 35, height = 35,  units = "cm" ,pointsize = 12,
        bg = "white", res = 1200, family = "")
    
    tree_mat_short =  ggplot(reducedTerms_short_edit, aes(area = score, fill = color, 
                                                          label = term, # country
                                                          subgroup = parentTerm)) + # continent
      scale_fill_identity() +
      # 1. Draw country borders and fill colors
      geom_treemap(colour = "#FFFFFFDD", start = "topleft") +
      # 2. Draw continent borders
      geom_treemap_subgroup_border(
        # colour = "#00000080", 
        start = "topleft", 
        colour = "gray100",
        size = 5) +
      # 3. Print continent text
      
      geom_treemap_subgroup_text(place = "centre", colour = "black",
                                 reflow = T,
                                 fontface = reducedTerms_short_edit$subgroup_font_fontface,
                                 # min.size = 0
                                 start = "topleft",
                                 grow = F,
                                 size= reducedTerms_short_edit$subgroup_font_size) +
      # 4. Print country text
      geom_treemap_text(colour = "#FFFFFFDD", 
                        place = "centre",
                        alpha = 0.5,
                        grow = F,
                        reflow = T,
                        start = "topleft") +
      # option
      theme(legend.position = 0)
    
    print(tree_mat_short)
    dev.off()
    
    
  } else if (length(short_gene) == 0 ) {
    long_gene_en = AnnotationDbi::select(org.Hs.eg.db, long_gene, 'ENTREZID', 'SYMBOL')[
      which(!is.na(AnnotationDbi::select(org.Hs.eg.db, long_gene, 'ENTREZID', 'SYMBOL')$ENTREZID)),]
     
    long_enrichGO = enrichGO(long_gene_en$ENTREZID ,
                             OrgDb = org.Hs.eg.db ,
                             keyType = "ENTREZID" , 
                             ont = "BP" ,
                             pvalueCutoff  = 0.05,
                             qvalueCutoff = 0.01 )
    
    long_enrichGO_df = as.data.frame(long_enrichGO)
    
    simMatrix_long <- calculateSimMatrix(long_enrichGO_df$ID,
                                         orgdb=org.Hs.eg.db,
                                         ont="BP",
                                         method="Rel")
    
    scores_long <- setNames(-log10(long_enrichGO_df$qvalue), long_enrichGO_df$ID)
    
    if (sum(is.na(scores_long)) != 0) {
      scores_long[which(is.na(scores_long))] = 0
    }
    
    if (sum(scores_long) == 0) {
      reducedTerms_long <- reduceSimMatrix(simMatrix_long,
                                           threshold=0.7,
                                           orgdb="org.Hs.eg.db")
    } else {
      reducedTerms_long <- reduceSimMatrix(simMatrix_long,
                                           scores_long,
                                           threshold=0.7,
                                           orgdb="org.Hs.eg.db")
    }
    
    score_long_head = head(unique(reducedTerms_long[order(reducedTerms_long$score,decreasing = T),]$parentTerm), n= 3)
    count_long_head = names(tail(sort(table(reducedTerms_long$parentTerm)), n= 3))
    total_long_imp = unique(c(total_long_imp, score_long_head, count_long_head))
    
    reducedTerms_long_edit = reducedTerms_long %>% 
      mutate(color = case_when( parentTerm %in% total_long_imp ~ "#339900",
                                .default = "#999999")) %>% 
      mutate(subgroup_font_size = case_when( parentTerm %in% total_long_imp ~ 30, .default = 10))%>% 
      mutate(subgroup_font_fontface = case_when( parentTerm %in% total_long_imp ~ "bold",.default = "plain"))
    
    # fiq
    
    png(filename = paste0(CancerType,"_tree_long_0.7_wo_ttest_not_exp_BP.png"),
        width = 35, height = 35,  units = "cm" ,pointsize = 12,
        bg = "white", res = 1200, family = "")
    
    tree_mat_long =  ggplot(reducedTerms_long_edit, aes(area = score, fill = color, 
                                                        label = term, # country
                                                        subgroup = parentTerm)) + # continent
      scale_fill_identity() +
      # 1. Draw country borders and fill colors
      geom_treemap(colour = "#FFFFFFDD", start = "topleft") +
      # 2. Draw continent borders
      geom_treemap_subgroup_border(
        # colour = "#00000080", 
        start = "topleft", 
        colour = "gray100",
        size = 5) +
      # 3. Print continent text
      
      geom_treemap_subgroup_text(place = "centre", colour = "black",
                                 reflow = T,
                                 fontface = reducedTerms_long_edit$subgroup_font_fontface,
                                 # min.size = 0
                                 start = "topleft",
                                 grow = F,
                                 size= reducedTerms_long_edit$subgroup_font_size) +
      # 4. Print country text
      geom_treemap_text(colour = "#FFFFFFDD", 
                        place = "centre",
                        alpha = 0.5,
                        grow = F,
                        reflow = T,
                        start = "topleft") +
      # option
      theme(legend.position = 0)
    
    print(tree_mat_long)
    dev.off()
    
  } else {
    print("I don't know")
  }
  

}  


