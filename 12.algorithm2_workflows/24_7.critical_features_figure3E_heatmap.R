#####

library(dplyr)
library(RColorBrewer)
library(readxl)
library(openxlsx)
library(igraph)
library(ggforce)
library(ggplot2)
library(graphlayouts)
library(ggraph)
library(igraph)
library(tidyr)

filepath = "/home/seokwon/nas/"
ref_path = paste0(filepath, "99.reference/")

Cancerlist = dir(paste0(filepath, "/00.data/filtered_TCGA/"))
Cancerlist = Cancerlist[c(-11,-12)]
# setwd("~/nas/04.Results/short_long/ttest")

setwd("~/nas/04.Results/short_long/quantile")
total_features = data.frame()

# for all
for (num_CancerType in Cancerlist) {
  
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  Cancernum = gsub('[.-]','',gsub('[a-zA-Z]','', num_CancerType))
  # call input
  
  cancer_bf_sl = read.csv(paste0(filepath,"04.Results/short_long/ttest_common/",CancerType, "_critical_features_short_long_common.csv"))
  
  tmp_features = data.frame(features = cancer_bf_sl$variable, sl = cancer_bf_sl$classification)
  tmp_features = tmp_features %>% mutate(value = case_when(sl == "long" ~ 5,
                                                 sl == "short" ~ -5,
                                                 sl == "common" ~ 1,
                                                 .default = NA))
  
  if (CancerType == "TCGA-CESC") {
    total_features = tmp_features %>% select(-sl)
  } else {
    total_features = merge(total_features , tmp_features %>% select(features, value) , by = "features", all = TRUE)
  }
  
}
library(RColorBrewer)
library(viridis)
colnames(total_features) = c("features" ,gsub('[.]','',gsub('\\d','', Cancerlist)) )
total_features[is.na(total_features)] = 0
rownames(total_features) = total_features$features
total_features$features = NULL

tmp = pheatmap(total_features, 
         # color = colorRampPalette(c('#C0392B','whitesmoke','#2471A3'))(11),
         color = c("#C0392B",
                   "#CA5E53",
                   "#D5847B",
                   "#DFA9A4",
                   "#EACFCC",
                   "#F5F5F5",
                   "#CBDAE4",
                   "#C0DEBE",
                   "#A6D2A2",
                   "#8CC787",
                   "#72BC6C"),
         show_rownames = F,
         # border_color = "black",
         # cutree_cols = 2,
         clustering_method = "ward.D2")
ggsave(file="supple_figure3E.svg", plot=tmp, width=10, height=10)

# ####
total_features_wo_common = total_features
total_features_wo_common[total_features_wo_common == 1] = 0
total_features_wo_common = total_features_wo_common[which(rowSums(total_features_wo_common) != 0),]

# pheatmap(total_features_wo_common, 
#          # color = colorRampPalette(c("red", "white", "navy"))(50),
#          color = colorRampPalette(c('#C0392B','grey','#2471A3'))(50),
#          cutree_cols = 2,
#          clustering_method = "ward.D2")
# 
# pheatmap(cor(total_features), 
#          # color = colorRampPalette(c("red", "white", "navy"))(50),
#          color = colorRampPalette(c('#C0392B','whitesmoke','#2471A3'))(40),
#          # cutree_cols = 2,
#          border_color = "white",
#          clustering_method = "average")

## only three long cancer 

total_features_3l = total_features %>% 
  select(`TCGA-UCEC`,`TCGA-OV`,`TCGA-BLCA`)
total_features_3l

values <- c(5, 1, 0, -5)
combinations <- expand.grid(replicate(3, values, simplify = FALSE))

sorted_combinations <- combinations %>%
  mutate(Sum = rowSums(.)) %>%
  arrange(desc(Sum))

total_features_3l_reorder = data.frame()

for (rownum in 1:nrow(sorted_combinations)) {
  tmp_condition = total_features_3l %>% filter(`TCGA-UCEC` == sorted_combinations[rownum,"Var1"],
                                               `TCGA-OV` == sorted_combinations[rownum,"Var2"], 
                                               `TCGA-BLCA` == sorted_combinations[rownum,"Var3"])
  
  total_features_3l_reorder = rbind(total_features_3l_reorder,tmp_condition)
}

long = pheatmap(total_features_3l_reorder %>% select(`TCGA-BLCA`,`TCGA-OV`,`TCGA-UCEC`), 
         # color = colorRampPalette(c("red", "white", "navy"))(50),
         # color = colorRampPalette(c('#C0392B', "#F5F5F5","#4DAF4A"))(11),
         color = c("#C0392B",
                   "#CA5E53",
                   "#D5847B",
                   "#DFA9A4",
                   "#EACFCC",
                   "#F5F5F5",
                   "#CBDAE4",
                   "#C0DEBE",
                   "#A6D2A2",
                   "#8CC787",
                   "#72BC6C"),
         # cutree_cols = 2,
         cluster_rows = F,
         cluster_cols = F,
         border_color = "white",
         clustering_method = "ward.D")

ggsave(file="figure3E.svg", plot=long, width=10, height=10)

# ###
# wo_3L
total_features_wo_common_wo_3l = total_features_wo_common %>%
  select(-`TCGA-UCEC`,-`TCGA-OV`,-`TCGA-BLCA`)
total_features_wo_common_wo_3l = total_features_wo_common_wo_3l[which(rowSums(total_features_wo_common_wo_3l) != 0),]

pheatmap(total_features_wo_common_wo_3l,
         # color = colorRampPalette(c("red", "white", "navy"))(50),
         color = c("#C0392B",
                   "#CA5E53",
                   "#D5847B",
                   "#DFA9A4",
                   "#EACFCC",
                   "#F5F5F5",
                   "#CBDAE4",
                   "#C0DEBE",
                   "#A6D2A2",
                   "#8CC787",
                   "#72BC6C"),
         # cutree_cols = 2,
         cluster_rows = T,
         border_color = "white",
         clustering_method = "complete")

ggsave(file="figure3E.svg", plot=long, width=10, height=10)
# 
# ### 
# total_features_wo_3l = total_features %>% select(-`TCGA-UCEC`,-`TCGA-OV`,-`TCGA-BLCA`)
# total_features_wo_3l = total_features_wo_3l[which(rowSums(total_features_wo_3l) != 0),]
# pheatmap(total_features_wo_3l, 
#          # color = colorRampPalette(c("red", "white", "navy"))(50),
#          color = colorRampPalette(c('#C0392B','grey','#2471A3'))(50),
#          cutree_cols = 2,
#          clustering_method = "ward.D2")
# ?pheatmap
# "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC)