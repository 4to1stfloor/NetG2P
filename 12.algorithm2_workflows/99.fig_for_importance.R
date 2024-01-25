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
  
  cancer_bf_top5 = cancer_bf %>% head(n= 5) %>% select(variable, minmax)

  
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

fig_total_top5_df = total_top5_df

tmp = as.data.frame(table(total_top5_df$pathway))
fig_total_top5_df$count <- tmp$Freq[match(fig_total_top5_df$pathway, tmp$Var1)]
fig_total_top5_df$cancer_type = gsub("TCGA-", "",fig_total_top5_df$cancer_type)
library(pheatmap)

fig_total_top5_df$rownames = paste0(fig_total_top5_df$pathway_name, " ; ", fig_total_top5_df$features)
rownames(fig_total_top5_df) = fig_total_top5_df$rownames

fig_total_top5_df

annotation_df <- data.frame(cancer_type = fig_total_top5_df$cancer_type,
                            count = fig_total_top5_df$count)

rownames(annotation_df) <- rownames(fig_total_top5_df)
library(RColorBrewer)
# Create a named color vector for the unique values of vital_status

color_map <- c(
  "UCEC" =  "#E64B35FF",
  "BRCA" = "#4DBBD5FF",
  "LGG" = "#00A087FF",
  "LUSC" = "#B09C85FF",
  "OV" = "#3C5488FF",
  "LUAD" = "#F39B7FFF",
  "LIHC" = "#8491B4FF",
  "STAD" = "#91D1C2FF",
  "BLCA" = "#7E6148FF",
  "CESC" = "#DC0000FF",
  "COADREAD" = "#FFFF99",
  "KIDNEY" = "#B15928"
)


# num_types = brewer.pal(length(unique(fig_total_top5_df$cancer_type)), "Paired")
# col_types = setNames(num_types, unique(fig_total_top5_df$cancer_type))
col_types = color_map[unique(fig_total_top5_df$cancer_type)]

num_count = brewer.pal(length(unique(fig_total_top5_df$count)), "Greys")
col_count = setNames(num_count, unique(fig_total_top5_df$count))

cancertypes_colors = list(cancer_type = col_types , count = col_count)
Colors = brewer.pal(9, "Blues")
tmp_heatmap = pheatmap::pheatmap(fig_total_top5_df %>% select(importance),
                   annotation_row = annotation_df,
                   annotation_colors = cancertypes_colors,
                   cellwidth = 10,
                   # block.size = c(2, 2),
                   color = Colors,
                   cluster_cols = F,
                   cluster_rows = F)



library(gridExtra)

t = pheatmap::pheatmap(fig_total_top5_df %>% 
                         filter(cancer_type %in% c("CESC","BLCA", "STAD", "LUAD", "LIHC", "OV"))%>% 
                         select(importance),
                       annotation_row = annotation_df,
                       annotation_colors = cancertypes_colors,
                       cellwidth = 10,
                       # block.size = c(2, 2),
                       color = Colors,
                       cluster_cols = F,
                       cluster_rows = F)

p = pheatmap::pheatmap(fig_total_top5_df %>% 
                         filter(cancer_type %in% c("LUSC","LGG", "BRCA", "UCEC", "COADREAD", "KIDNEY"))%>% 
                         select(importance),
                       annotation_row = annotation_df,
                       annotation_colors = cancertypes_colors,
                       cellwidth = 10,
                       # block.size = c(2, 2),
                       color = Colors,
                       cluster_cols = F,
                       cluster_rows = F,)
library(grid)

print(p)
p$gtable$grobs[[5]]$gp = gpar(fontface = "bold")
grid = grid.arrange(arrangeGrob(grobs= list(t[[4]],p[[4]]) , ncol = 2))

ggsave(file = "importance_edit.svg", grid, width=13, height=10, device = svg)


tmp_p = fig_total_top5_df %>% 
  filter(cancer_type %in% c("LUSC","LGG", "BRCA", "UCEC", "COADREAD", "KIDNEY"))%>% 
  select(importance)

tmp_t = fig_total_top5_df %>% 
  filter(cancer_type %in% c("CESC","BLCA", "STAD", "LUAD", "LIHC", "OV"))%>% 
  select(importance)

write.csv(as.data.frame(rownames(tmp_p)), "right_name.csv")
write.csv(as.data.frame(rownames(tmp_t)), "left_name.csv")
write.csv(as.data.frame(unique(fig_total_top5_df$cancer_type)), "cancer_name.csv")

