library(ggplot2)
library(dplyr)
library(ggpubr)
library(RColorBrewer)
filepath = "/home/seokwon/nas/"
Cancerlist = dir(paste0(filepath, "/00.data/filtered_TCGA/"))
Cancername = c()
for (num_CancerType in Cancerlist) {
  
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  Cancername = c(Cancername,strsplit(num_CancerType, "-")[[1]][2])

}

hsize = 4
Cancername = Cancername[c(-7,-8,-9,-11,-12)]
df_cancername = as.data.frame(Cancername)

df_cancername$value = 1

# Compute percentages
df_cancername$fraction <- df_cancername$value / sum(df_cancername$value)

# Compute the cumulative percentages (top of each rectangle)
df_cancername$ymax <- cumsum(df_cancername$fraction)

# Compute the bottom of each rectangle
df_cancername$ymin <- c(0, head(df_cancername$ymax, n=-1))

# Compute label position
df_cancername$labelPosition <- (df_cancername$ymax + df_cancername$ymin) / 2

color_map <- c(
  "UCEC" =  "#E64B35FF",
  "BRCA" = "#4DBBD5FF",
  "LGG" = "#00A087FF",
  "OV" = "#3C5488FF",
  "LUAD" = "#F39B7FFF",
  "LIHC" = "#8491B4FF",
  "STAD" = "#91D1C2FF",
  "BLCA" = "#DC0000FF",
  "CESC" = "#7E6148FF"
)

ggplot(df_cancername, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill= Cancername, color = Cancername)) +
  geom_rect() +
  geom_text( x=4.8, aes(y=labelPosition, label=Cancername), size=5) + # x here controls label position (inner / outer)
  # scale_fill_brewer(palette=3) +
  # scale_color_brewer(palette=3) +
  scale_color_manual(values = color_map, aesthetics = "fill")+
  scale_color_manual(values = color_map, aesthetics = "color")+
  coord_polar(theta="y") +
  xlim(c(-1, 4)) +
  theme_void() +
  # scale_fill_manual(
  #   labels = df_cancername$Cancername,
  #   values = getPalette(colorcount)
  # ) +
  theme(legend.position = "none")


#################################

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
library(tidygraph)
library(tidyverse)

filepath = "/home/seokwon/nas/"
ref_path = paste0(filepath, "99.reference/")

Cancerlist = dir(paste0(filepath, "/00.data/filtered_TCGA/"))
surv_total_results = read_xlsx("~/nas/04.Results/Total_results_survpval2.xlsx")
# setwd("~/nas/04.Results/short_long/ttest")
Cancerlist = Cancerlist[c(-7,-11,-12)]
total_features = list()

# for all
for (num_CancerType in Cancerlist) {
  
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  Cancernum = gsub('[.-]','',gsub('[a-zA-Z]','', num_CancerType))
  # call input
  
  cancer_bf = read.csv(paste0(filepath,"04.Results/bestfeatures/",CancerType, "_critical_features.csv"))
  
  # cancer_bf_cut = cancer_bf[1:surv_total_results[which(surv_total_results$CancerType == CancerType),]$num_of_features,]
  total_features[[CancerType]] = cancer_bf$variable
  
}

most_common = read_xlsx("~/nas/04.Results/critical_features/most_common_with_short_long_critical_features.xlsx")
most_common_max = most_common[which(most_common$count == max(most_common$count)),]

# most_common_max_filt = most_common_max %>%
#   filter(nchar(features) > 3)

fit_features = surv_total_results %>% 
  filter(pval < 0.05 ) %>% 
  filter(bias == "balanced") %>%
  arrange(desc(num_of_features))

# c(fit_features$CancerType[4:9] , fit_features$CancerType[c(2,3)])
max_comb = list()
combn(c(fit_features$CancerType), 3)
for (n in 3:length(c(fit_features$CancerType))) {
  print(n)
  tmp_comb <- combn(c(fit_features$CancerType), n)
  
  for (nc in 1:ncol(tmp_comb)) {
    # print(nc)
    for (shared_cancer in tmp_comb[,nc]) {
      # print(shared_cancer)
      if (shared_cancer == tmp_comb[1,nc]) {
        tmp_common = total_features[[shared_cancer]]
        # print(tmp_common)
      } else {
        tmp_common = intersect(tmp_common , total_features[[shared_cancer]])
        # print(tmp_common)
      }
      
    }
    
    # print(tmp_common)
    
    if (is.null(tmp_common)) {
      next
    }
    
    if (length(max_comb) == 0 ) {
      max_comb[[paste(tmp_comb[,nc], collapse = ",")]] = tmp_common
      # print(max_comb)
      
    } else if (length(max_comb[[1]]) < length(tmp_common)) {
      # print(max_comb)
      # print(tmp_common)
      max_comb[[paste(tmp_comb[,nc], collapse = ",")]] = tmp_common
      # print(max_comb)
      max_comb[[1]] = NULL
      # print(max_comb)
      # pddddd
    }
    
  }
}

tmp_common = max_comb[[1]]

## if add TCGA-LIHC
# cancer_fit = fit_features$CancerType[1:4]
# for (shared_cancer in cancer_fit) {
# 
#   if (shared_cancer == cancer_fit[1]) {
#     tmp_common = total_features[[shared_cancer]]
#   } else {
#     tmp_common = intersect(tmp_common , total_features[[shared_cancer]])
#   }
# 
# }


# max_features_shared = table(unlist(str_split(most_common_max_filt$which_cancer, ", ")))
# for (shared_cancer in names(max_features_shared[max_features_shared == max(max_features_shared)])) {
# 
#   if (shared_cancer == names(max_features_shared[max_features_shared == max(max_features_shared)])[1]) {
#     tmp_common = total_features[[shared_cancer]]
#   } else {
#     tmp_common = intersect(tmp_common , total_features[[shared_cancer]])
#   }
# 
# }

nodes = c()
edges = data.frame(matrix(ncol = 2))
colnames(edges) = c("from","to")

for (critical_features in tmp_common) {
  nodes = unique(c(nodes , paste0("P",strsplit(critical_features, "P")[[1]][2]) , paste0("P",strsplit(critical_features, "P")[[1]][3])))
  
  if (is.na(strsplit(critical_features, "P")[[1]][3] )) {
    edges[critical_features,"from"] = paste0("P",strsplit(critical_features, "P")[[1]][2]) 
    edges[critical_features,"to"] = paste0("P",strsplit(critical_features, "P")[[1]][2])
  } else {
    edges[critical_features,"from"] = paste0("P",strsplit(critical_features, "P")[[1]][2]) 
    edges[critical_features,"to"] = paste0("P",strsplit(critical_features, "P")[[1]][3])
  }
  
}
edges = edges[-1,]
rownames(edges) = NULL

# 
setwd("~/nas/04.Results/critical_features/cancer_specific_critical_features/")

png(filename = paste0("raw_3large_features_cancer_common_critical_features.png"),
    width = 10, height = 10,  units = "cm" ,pointsize = 12,
    bg = "white", res = 1200, family = "")

raw_network = ggraph(edges %>% as_tbl_graph(), layout = "fr") + # 레이아웃
  geom_edge_link(color = "gray50",          # 엣지 색깔
                 alpha = 0.5) +             # 엣지 명암

  geom_node_point(color = "lightcoral",     # 노드 색깔
                  size = 5) +               # 노드 크기

  geom_node_text(aes(label = name),         # 텍스트 표시
                 repel = T,                 # 노드밖 표시
                 size = 5,                  # 텍스트 크기
                 family = "nanumgothic") +  # 폰트

  theme_graph()                             # 배경 삭제

print(raw_network)

dev.off()
# 
# png(filename = paste0(CancerType,"_centrality_6large_features_cancer_common_critical_features.png"),
#     width = 10, height = 10,  units = "cm" ,pointsize = 12,
#     bg = "white", res = 1200, family = "")
# 
# central_network = ggraph( edges %>%
#                             as_tbl_graph(directed = F) %>%
#                             mutate(centrality = centrality_degree(),        # 연결 중심성
#                                    group = as.factor(group_louvain())), layout = "fr") + # 레이아웃
#   geom_edge_link(color = "gray50",          # 엣지 색깔
#                  alpha = 0.5) +             # 엣지 명암
# 
#   geom_node_point(aes(size = centrality,    # 노드 크기
#                       color = group),       # 노드 색깔
#                   show.legend = F) +        # 범례 삭제
#   scale_size(range = c(3, 10)) +
# 
#   geom_node_text(aes(label = name),         # 텍스트 표시
#                  repel = T,                 # 노드밖 표시
#                  size = 5,                  # 텍스트 크기
#                  family = "nanumgothic") +  # 폰트
# 
#   theme_graph()                             # 배경 삭제
# 
# print(central_network)
# 
# dev.off()
# 
