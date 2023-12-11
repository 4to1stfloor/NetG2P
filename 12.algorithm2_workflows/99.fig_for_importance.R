library(tidyverse)
library(ggplot2)
library(ggridges)
library(reshape2)
library(readxl)

filepath = "/home/seokwon/nas/"
ref_path = paste0(filepath, "99.reference/")

name_pathway = read_xlsx(paste0(ref_path,"pathway_name.xlsx"))

Cancerlist = dir(paste0(filepath, "/00.data/filtered_TCGA/"))

# 
# # for fic
# fig_path = "~/nas/04.Results/drug/depmap/"
# if(!dir.exists(fig_path)){
#   dir.create(fig_path)
#   print(paste0("Created folder: ", fig_path))
# } else {
#   print(paste0("Folder already exists: ", fig_path))
# }
# setwd(fig_path)

# top5 importance 및 이름 heatmap 
# 
# num_CancerType = "04.TCGA-CESC"

total_top5_df = data.frame()
for (num_CancerType in Cancerlist) {
  
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  ModelType = gsub('TCGA-','', CancerType)
  cancer_bf = read.csv(paste0("~/nas/04.Results/short_long/ttest_common/",CancerType,"_critical_features_short_long_common.csv"))
  
  cancer_bf_top5 = cancer_bf %>% head(, n= 5) %>% select(variable, minmax)

  
  # short
  each_top5_df = data.frame()
  # top5 = "P36P51"
  for (top5 in cancer_bf_top5$variable) {
    count <- str_count(top5, "P")

    if (count == 1) {
      tmp_pathway = data.frame( features = top5,
                                pathway = top5,
                                pathway_name = name_pathway[which(name_pathway$num_pathway == top5),]$Kegg_pathway_name,
                                importance = cancer_bf_top5[which(cancer_bf_top5$variable == top5),]$minmax)
      
    } else {
      first_path = paste0("P",str_split(top5, "P")[[1]][2])
      second_path = paste0("P",str_split(top5, "P")[[1]][3])
      
      tmp_pathway_first = data.frame( features = top5, 
                                      pathway = first_path,
                                      pathway_name = name_pathway[which(name_pathway$num_pathway == first_path),]$Kegg_pathway_name,
                                      importance = cancer_bf_top5[which(cancer_bf_top5$variable == top5),]$minmax)
      
      tmp_pathway_second = data.frame( features = top5, 
                                       pathway = second_path,
                                       pathway_name = name_pathway[which(name_pathway$num_pathway == second_path),]$Kegg_pathway_name,
                                       importance = cancer_bf_top5[which(cancer_bf_top5$variable == top5),]$minmax)
      
      tmp_pathway = rbind(tmp_pathway_first,tmp_pathway_second)
    }
    
    each_top5_df = rbind(each_top5_df, tmp_pathway)
    
  }
  
  each_top5_df$cancer_type = CancerType
  
  total_top5_df = rbind(total_top5_df, each_top5_df)
 
  
}

table(total_top5_df$pathway)[table(total_top5_df$pathway) == max(table(total_top5_df$pathway))]

pheatmap::pheatmap(total_top5_df %>% select(features,importance))
fig_total_top5_df = total_top5_df

tmp = as.data.frame(table(total_top5_df$pathway))
fig_total_top5_df$count <- tmp$Freq[match(fig_total_top5_df$pathway, tmp$Var1)]
fig_total_top5_df$cancer_type = gsub("TCGA-", "",fig_total_top5_df$cancer_type)
library(pheatmap)

fig_total_top5_df$rownames = paste0(fig_total_top5_df$features, "_", fig_total_top5_df$pathway_name)
rownames(fig_total_top5_df) = fig_total_top5_df$rownames

annotation_df <- data.frame(cancer_type = fig_total_top5_df$cancer_type,
                            count = fig_total_top5_df$count)

rownames(annotation_df) <- rownames(fig_total_top5_df)
library(RColorBrewer)
# Create a named color vector for the unique values of vital_status

num_types = brewer.pal(length(unique(fig_total_top5_df$cancer_type)), "Paired")
col_types = setNames(num_types, unique(fig_total_top5_df$cancer_type))

num_count = brewer.pal(length(unique(fig_total_top5_df$count)), "Set1")
col_count = setNames(num_count, unique(fig_total_top5_df$count))

cancertypes_colors = list(cancer_type = col_types , count = col_count)

pheatmap::pheatmap(fig_total_top5_df %>% select(importance),
                   annotation_row = annotation_df,
                   annotation_colors = cancertypes_colors,
                   # block.size = c(2, 2),
                   cluster_cols = F,
                   cluster_rows = F)
