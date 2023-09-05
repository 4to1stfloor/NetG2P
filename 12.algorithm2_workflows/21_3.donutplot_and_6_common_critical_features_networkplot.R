library(ggplot2)
library(dplyr)
library(ggpubr)
library(RColorBrewer)

Cancerlist = dir(paste0(filepath, "/00.data/filtered_TCGA/"))
Cancername = c()
for (num_CancerType in Cancerlist) {
  
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  Cancername = c(Cancername,strsplit(num_CancerType, "-")[[1]][2])

}

hsize = 4
Cancername = Cancername[c(-7,-8,-9,-11)]
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

colorcount <- nrow(df_cancername)
getPalette <- colorRampPalette(brewer.pal(nrow(df_cancername) , "Paired"))

p = ggplot(df_cancername, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Cancername)) +
  geom_rect() +
  geom_text( x=4.8, aes(y=labelPosition, label=Cancername), size=5) + # x here controls label position (inner / outer)
  # scale_fill_brewer(palette=3) +
  # scale_color_brewer(palette=3) +
  coord_polar(theta="y") +
  xlim(c(-1, 4)) +
  theme_void() +
  # scale_fill_manual(
  #   labels = df_cancername$Cancername,
  #   values = getPalette(colorcount)
  # ) +
  theme(legend.position = "none")

set_palette(p , "npg")

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

filepath = "/home/seokwon/nas/"
ref_path = paste0(filepath, "99.reference/")

Cancerlist = dir(paste0(filepath, "/00.data/filtered_TCGA/"))
surv_total_results = read_xlsx("~/nas/04.Results/Total_results_survpval2.xlsx")
# setwd("~/nas/04.Results/short_long/ttest")
Cancerlist = Cancerlist[c(-7,-11)]
total_features = list()

# for all
for (num_CancerType in Cancerlist) {
  
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  Cancernum = gsub('[.-]','',gsub('[a-zA-Z]','', num_CancerType))
  # call input
  
  cancer_bf = read.csv(paste0(filepath,"04.Results/bestfeatures/",CancerType, "_best_features.csv"))
  
  cancer_bf_cut = cancer_bf[1:surv_total_results[which(surv_total_results$CancerType == CancerType),]$num_of_features,]
  total_features[[CancerType]] = cancer_bf_cut$variable
  
}


total_common = list()
Cancername

for (n in c(2,3,5,6,9)) {
  print(n)
  if (n == 2) {
    tmp_common = total_features[[n]]
  } else {
    tmp_common = intersect(tmp_common , total_features[[n]])
  }
  
  
}

tmp_common

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

png(filename = paste0(CancerType,"_raw_6large_features_cancer_common_critical_features.png"),
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

png(filename = paste0(CancerType,"_centrality_6large_features_cancer_common_critical_features.png.png"),
    width = 10, height = 10,  units = "cm" ,pointsize = 12,
    bg = "white", res = 1200, family = "")

central_network = ggraph( edges %>%
                            as_tbl_graph(directed = F) %>%
                            mutate(centrality = centrality_degree(),        # 연결 중심성
                                   group = as.factor(group_louvain())), layout = "fr") + # 레이아웃
  geom_edge_link(color = "gray50",          # 엣지 색깔
                 alpha = 0.5) +             # 엣지 명암
  
  geom_node_point(aes(size = centrality,    # 노드 크기
                      color = group),       # 노드 색깔
                  show.legend = F) +        # 범례 삭제
  scale_size(range = c(3, 10)) +      
  
  geom_node_text(aes(label = name),         # 텍스트 표시
                 repel = T,                 # 노드밖 표시
                 size = 5,                  # 텍스트 크기
                 family = "nanumgothic") +  # 폰트
  
  theme_graph()                             # 배경 삭제

print(central_network)

dev.off()

