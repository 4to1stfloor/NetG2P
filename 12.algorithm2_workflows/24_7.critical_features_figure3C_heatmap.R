#####
library(pheatmap)
library(dplyr)
library(RColorBrewer)
library(readxl)
library(openxlsx)
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
    total_features = tmp_features %>% dplyr::select(-sl)
  } else {
    total_features = merge(total_features , tmp_features %>% dplyr::select(features, value) , by = "features", all = TRUE)
  }
  
}

cal_total_features = total_features
library(RColorBrewer)
library(viridis)
colnames(total_features) = c("features" ,gsub('TCGA-','',gsub('[.]','',gsub('\\d','', Cancerlist))) )
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
                   "#D3DFE5",
                   "#C0DEBE",
                   "#A6D2A2",
                   "#8CC787",
                   "#72BC6C"),
         show_rownames = F,
         # border_color = "black",
         # cutree_cols = 2,
         clustering_method = "ward.D2")
ggsave(file="figure3E_total.svg", plot=tmp, width=10, height=10)

#####
tmp_total = total_features
values <- c(5, 1, 0, -5)
colnames(total_features)
combinations <- expand.grid(replicate(ncol(tmp_total), values, simplify = FALSE))
colnames(combinations) = colnames(total_features)
sorted_combinations <- combinations %>%
  mutate(Sum = rowSums(.)) %>%
  arrange(desc(Sum))
sorted_combinations

###
tmp_long = tmp_total %>%
  filter_all(all_vars(. %in% c(0, 5)))

values <- c(5, 0)
combinations <- expand.grid(replicate(ncol(tmp_long), values, simplify = FALSE))
colnames(combinations) = colnames(tmp_long)
sorted_combinations <- combinations %>%
  mutate(Sum = rowSums(.)) %>%
  arrange(desc(Sum))

tmp_long_reorder = data.frame()

for (rownum in 1:nrow(sorted_combinations)) {
  tmp_condition = tmp_long %>% filter(CESC == sorted_combinations[rownum,"CESC"],
                                       BLCA == sorted_combinations[rownum,"BLCA"],
                                       STAD == sorted_combinations[rownum,"STAD"],
                                       LUAD == sorted_combinations[rownum,"LUAD"],
                                       LIHC == sorted_combinations[rownum,"LIHC"],
                                       OV == sorted_combinations[rownum,"OV"],
                                       LUSC == sorted_combinations[rownum,"LUSC"],
                                       LGG == sorted_combinations[rownum,"LGG"],
                                       BRCA == sorted_combinations[rownum,"BRCA"],
                                       UCEC == sorted_combinations[rownum,"UCEC"]
                                       )
  
  tmp_long_reorder = rbind(tmp_long_reorder,tmp_condition)
}

##
tmp_long_common <- tmp_total %>%
  filter_all(all_vars(. %in% c(0, 5, 1))) %>%
  filter((rowSums(. == 5) > 0) & (rowSums(. == 1) > 0) & (rowSums(. == 0) > 0))

values <- c(5, 1, 0)
combinations <- expand.grid(replicate(ncol(tmp_long_common), values, simplify = FALSE))
colnames(combinations) = colnames(tmp_long_common)
sorted_combinations <- combinations %>%
  mutate(Sum = rowSums(.)) %>%
  arrange(desc(Sum))

tmp_long_common_reorder = data.frame()

for (rownum in 1:nrow(sorted_combinations)) {
  tmp_condition = tmp_long_common %>% filter(CESC == sorted_combinations[rownum,"CESC"],
                                       BLCA == sorted_combinations[rownum,"BLCA"],
                                       STAD == sorted_combinations[rownum,"STAD"],
                                       LUAD == sorted_combinations[rownum,"LUAD"],
                                       LIHC == sorted_combinations[rownum,"LIHC"],
                                       OV == sorted_combinations[rownum,"OV"],
                                       LUSC == sorted_combinations[rownum,"LUSC"],
                                       LGG == sorted_combinations[rownum,"LGG"],
                                       BRCA == sorted_combinations[rownum,"BRCA"],
                                       UCEC == sorted_combinations[rownum,"UCEC"]
  )
  
  tmp_long_common_reorder = rbind(tmp_long_common_reorder,tmp_condition)
}

##
tmp_common = tmp_total %>%
  filter_all(all_vars(. %in% c(0, 1)))

values <- c( 1, 0)
combinations <- expand.grid(replicate(ncol(tmp_common), values, simplify = FALSE))
colnames(combinations) = colnames(tmp_common)
sorted_combinations <- combinations %>%
  mutate(Sum = rowSums(.)) %>%
  arrange(desc(Sum))

tmp_common_reorder = data.frame()

for (rownum in 1:nrow(sorted_combinations)) {
  tmp_condition = tmp_common %>% filter(CESC == sorted_combinations[rownum,"CESC"],
                                       BLCA == sorted_combinations[rownum,"BLCA"],
                                       STAD == sorted_combinations[rownum,"STAD"],
                                       LUAD == sorted_combinations[rownum,"LUAD"],
                                       LIHC == sorted_combinations[rownum,"LIHC"],
                                       OV == sorted_combinations[rownum,"OV"],
                                       LUSC == sorted_combinations[rownum,"LUSC"],
                                       LGG == sorted_combinations[rownum,"LGG"],
                                       BRCA == sorted_combinations[rownum,"BRCA"],
                                       UCEC == sorted_combinations[rownum,"UCEC"]
  )
  
  tmp_common_reorder = rbind(tmp_common_reorder,tmp_condition)
}

##

tmp_common_short = tmp_total %>%
  filter_all(all_vars(. %in% c(0, -5, 1))) %>%
  filter((rowSums(. == -5) > 0) & (rowSums(. == 1) > 0) & (rowSums(. == 0) > 0))

values <- c(-5, 1, 0)
combinations <- expand.grid(replicate(ncol(tmp_common_short), values, simplify = FALSE))
colnames(combinations) = colnames(tmp_common_short)
sorted_combinations <- combinations %>%
  mutate(Sum = rowSums(.)) %>%
  arrange(desc(Sum))

tmp_common_short_reorder = data.frame()

for (rownum in 1:nrow(sorted_combinations)) {
  tmp_condition = tmp_common_short %>% filter(CESC == sorted_combinations[rownum,"CESC"],
                                       BLCA == sorted_combinations[rownum,"BLCA"],
                                       STAD == sorted_combinations[rownum,"STAD"],
                                       LUAD == sorted_combinations[rownum,"LUAD"],
                                       LIHC == sorted_combinations[rownum,"LIHC"],
                                       OV == sorted_combinations[rownum,"OV"],
                                       LUSC == sorted_combinations[rownum,"LUSC"],
                                       LGG == sorted_combinations[rownum,"LGG"],
                                       BRCA == sorted_combinations[rownum,"BRCA"],
                                       UCEC == sorted_combinations[rownum,"UCEC"]
  )
  
  tmp_common_short_reorder = rbind(tmp_common_short_reorder,tmp_condition)
}

##
tmp_short = tmp_total %>%
  filter_all(all_vars(. %in% c(0, -5)))

values <- c(-5, 0)
combinations <- expand.grid(replicate(ncol(tmp_short), values, simplify = FALSE))
colnames(combinations) = colnames(tmp_short)
sorted_combinations <- combinations %>%
  mutate(Sum = rowSums(.)) %>%
  arrange(desc(Sum))

tmp_short_reorder = data.frame()

for (rownum in 1:nrow(sorted_combinations)) {
  tmp_condition = tmp_short %>% filter(CESC == sorted_combinations[rownum,"CESC"],
                                              BLCA == sorted_combinations[rownum,"BLCA"],
                                              STAD == sorted_combinations[rownum,"STAD"],
                                              LUAD == sorted_combinations[rownum,"LUAD"],
                                              LIHC == sorted_combinations[rownum,"LIHC"],
                                              OV == sorted_combinations[rownum,"OV"],
                                              LUSC == sorted_combinations[rownum,"LUSC"],
                                              LGG == sorted_combinations[rownum,"LGG"],
                                              BRCA == sorted_combinations[rownum,"BRCA"],
                                              UCEC == sorted_combinations[rownum,"UCEC"]
  )
  
  tmp_short_reorder = rbind(tmp_short_reorder,tmp_condition)
}

###
tmp_last = tmp_total %>%
  filter(!rownames(.) %in% c(rownames(tmp_long),
                             rownames(tmp_long_common),
                             rownames(tmp_common),
                             rownames(tmp_common_short),
                             rownames(tmp_short)))

values <- c(-5, 0 ,1 , 5)
combinations <- expand.grid(replicate(ncol(tmp_last), values, simplify = FALSE))
colnames(combinations) = colnames(tmp_last)
sorted_combinations <- combinations %>%
  mutate(Sum = rowSums(.)) %>%
  arrange(desc(Sum))

tmp_last_reorder = data.frame()

for (rownum in 1:nrow(sorted_combinations)) {
  tmp_condition = tmp_last %>% filter(CESC == sorted_combinations[rownum,"CESC"],
                                       BLCA == sorted_combinations[rownum,"BLCA"],
                                       STAD == sorted_combinations[rownum,"STAD"],
                                       LUAD == sorted_combinations[rownum,"LUAD"],
                                       LIHC == sorted_combinations[rownum,"LIHC"],
                                       OV == sorted_combinations[rownum,"OV"],
                                       LUSC == sorted_combinations[rownum,"LUSC"],
                                       LGG == sorted_combinations[rownum,"LGG"],
                                       BRCA == sorted_combinations[rownum,"BRCA"],
                                       UCEC == sorted_combinations[rownum,"UCEC"]
                                      )
  
  tmp_last_reorder = rbind(tmp_last_reorder,tmp_condition)
}

tmp_total_reorder_filt = rbind(tmp_long_reorder,tmp_long_common_reorder, tmp_common_reorder, tmp_common_short_reorder, tmp_short_reorder, tmp_last_reorder)

tmp = pheatmap(tmp_total_reorder_filt, 
               # color = colorRampPalette(c('#C0392B','whitesmoke','#2471A3'))(11),
               color = c("#C0392B",
                         "#CA5E53",
                         "#D5847B",
                         "#DFA9A4",
                         "#EACFCC",
                         "#F5F5F5",
                         "#D3DFE5",
                         "#C0DEBE",
                         "#A6D2A2",
                         "#8CC787",
                         "#72BC6C"),
               show_rownames = F,
               cluster_rows = F,
               # border_color = "black",
               # cutree_cols = 2,
               clustering_method = "ward.D2")
ggsave(file="figure3E_total.svg", plot=tmp, width=10, height=10)

####

values <- c( 5,1, 0 ,-5)
combinations <- expand.grid(replicate(ncol(tmp_total), values, simplify = FALSE))
colnames(combinations) = colnames(tmp_total)
sorted_combinations <- combinations %>%
  mutate(Sum = rowSums(.)) %>%
  arrange(desc(Sum))

tmp_total_reorder = data.frame()

for (rownum in 1:nrow(sorted_combinations)) {
  tmp_condition = tmp_total %>% filter(CESC == sorted_combinations[rownum,"CESC"],
                                       BLCA == sorted_combinations[rownum,"BLCA"],
                                       STAD == sorted_combinations[rownum,"STAD"],
                                       LUAD == sorted_combinations[rownum,"LUAD"],
                                       LIHC == sorted_combinations[rownum,"LIHC"],
                                       OV == sorted_combinations[rownum,"OV"],
                                       LUSC == sorted_combinations[rownum,"LUSC"],
                                       LGG == sorted_combinations[rownum,"LGG"],
                                       BRCA == sorted_combinations[rownum,"BRCA"],
                                       UCEC == sorted_combinations[rownum,"UCEC"]
                                       )

  tmp_total_reorder = rbind(tmp_total_reorder,tmp_condition)
}

tmp = pheatmap(tmp_total_reorder, 
               # color = colorRampPalette(c('#C0392B','whitesmoke','#2471A3'))(11),
               color = c("#C0392B",
                         "#CA5E53",
                         "#D5847B",
                         "#DFA9A4",
                         "#EACFCC",
                         "#F5F5F5",
                         "#D3DFE5",
                         "#C0DEBE",
                         "#A6D2A2",
                         "#8CC787",
                         "#72BC6C"),
               show_rownames = F,
               cluster_rows = F,
               # border_color = "black",
               # cutree_cols = 2,
               clustering_method = "ward.D2")
ggsave(file="figure3E_total.svg", plot=tmp, width=10, height=10)

## 값조정해서 다시해야함
tmp_total[tmp_total == 5] = 3
tmp_total[tmp_total == -5] = -20

tmp_total_reorder = tmp_total %>%
  mutate(Sum = rowSums(.)) %>%
  arrange(desc(Sum))

tmp_total_reorder[tmp_total_reorder == 3] = 5
tmp_total_reorder[tmp_total_reorder == -20] = -5

tmp = pheatmap(tmp_total_reorder %>% select(-Sum), 
               # color = colorRampPalette(c('#C0392B','whitesmoke','#2471A3'))(11),
               color = c("#C0392B",
                         "#CA5E53",
                         "#D5847B",
                         "#DFA9A4",
                         "#EACFCC",
                         "#F5F5F5",
                         "#D3DFE5",
                         "#C0DEBE",
                         "#A6D2A2",
                         "#8CC787",
                         "#72BC6C"),
               show_rownames = F,
               cluster_rows = F,
               # border_color = "black",
               # cutree_cols = 2,
               clustering_method = "ward.D2")
ggsave(file="figure3C_total.svg", plot=tmp, width=10, height=10)

# ####
# total_features_wo_common = total_features
# total_features_wo_common[total_features_wo_common == 1] = 0
# total_features_wo_common = total_features_wo_common[which(rowSums(total_features_wo_common) != 0),]

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
                   "#D3DFE5",
                   "#C0DEBE",
                   "#A6D2A2",
                   "#8CC787",
                   "#72BC6C"),
         # cutree_cols = 2,
         cluster_rows = F,
         cluster_cols = F,
         border_color = "white",
         clustering_method = "ward.D")

ggsave(file="figure3E_long.svg", plot=long, width=10, height=10)

# ###
# wo_3L
total_features_wo_3l = total_features %>%
  select(-`TCGA-UCEC`,-`TCGA-OV`,-`TCGA-BLCA`)
total_features_wo_3l = total_features_wo_3l[which(rowSums(total_features_wo_3l) != 0),]

short = pheatmap(total_features_wo_3l,
         # color = colorRampPalette(c("red", "white", "navy"))(50),
         color = c("#C0392B",
                   "#CA5E53",
                   "#D5847B",
                   "#DFA9A4",
                   "#EACFCC",
                   "#F5F5F5",
                   "#D3DFE5",
                   "#C0DEBE",
                   "#A6D2A2",
                   "#8CC787",
                   "#72BC6C"),
         # cutree_cols = 2,
         cluster_rows = T,
         border_color = "white",
         clustering_method = "ward.D")

ggsave(file="figure3E_short.svg", plot=short, width=10, height=10)
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


# calculate ratio of features
colnames(cal_total_features) = c("features" ,gsub('[.]','',gsub('\\d','', Cancerlist)) )
rownames(cal_total_features) = cal_total_features$features
cal_total_features$features = NULL

node_classification = data.frame()
for (tmp_name in colnames(cal_total_features)) {
  tmp_df = cal_total_features %>% 
    dplyr::select(tmp_name) %>% dplyr::filter(!is.na(.))

  long_ratio = sum(tmp_df == 5) / nrow(tmp_df)
  common_ratio = sum(tmp_df == 1) / nrow(tmp_df)
  short_ratio = sum(tmp_df == -5) / nrow(tmp_df)
  
  tmp_classi = data.frame(long = long_ratio, short = short_ratio, common = common_ratio)
  rownames(tmp_classi) = tmp_name
  node_classification = rbind(node_classification,tmp_classi)
  
}


#### sort by long 
####


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
