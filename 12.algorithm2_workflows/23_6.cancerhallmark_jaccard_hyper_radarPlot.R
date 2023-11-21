## need critical features and cancer hallmark
## sankeyNetwork
# library(networkD3)
library(magrittr)
library(stringr)
library(data.table)
library(pheatmap)
library(readxl)
library(dplyr)
library(readxl)
library(tidyverse)

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

surv_total_results = read_xlsx("~/nas/04.Results/Total_results_survpval2.xlsx")
cancerhallmark = read_xlsx(paste0(ref_path, "kegg_gene_set_w_cancer_hallmarks_edit.xlsx"))  

# for fic
fig_path = paste0(filepath,"04.Results/sankey_cancerhallmark/")
if(!dir.exists(fig_path)){
  dir.create(fig_path)
  print(paste0("Created folder: ", fig_path))
} else {
  print(paste0("Folder already exists: ", fig_path))
}
setwd(fig_path)
tmp_hallmark_links = data.frame()

total_hallmark_genes = list()

for (hall_num in unique(cancerhallmark$Cancer_hallmark_no.)) {
  each_hallmark = cancerhallmark %>% filter(Cancer_hallmark_no. == hall_num)
  hallmark_pathway = paste0("P" , each_hallmark$No.)
  hallmark_genes = single_genes %>% filter(Pathway %in% hallmark_pathway)

  total_hallmark_genes[[gsub(" ", "_",  unique(cancerhallmark[which((cancerhallmark$Cancer_hallmark_no. == hall_num)),]$Hallmark_Name))]] =  unique(hallmark_genes$Genes)

} 

max_length <- max(sapply(total_hallmark_genes, length))

total_hallmark_genes = lapply(total_hallmark_genes, function(x) {
  if (length(x) < max_length) {
    c(x, rep(NA, max_length - length(x)))
  } else {
    x
  }
})


df_total_hallmark_genes = do.call(data.frame, total_hallmark_genes)
# setwd("~/nas/99.reference/")
# write.xlsx(df_total_hallmark_genes, "gene_set_w_cancer_hallmarks_edit.xlsx")
num_CancerType = "10.TCGA-BLCA"
# for all
for (num_CancerType in Cancerlist) {
  
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  
  # call input
  
  # t.test 한거
  # short_long_features = read_xlsx(paste0(filepath, "04.Results/short_long/",CancerType, "_best_features_short_long.xlsx"))
  cf_sl = read.csv(paste0(filepath, "04.Results/short_long/ttest_common/",CancerType,"_critical_features_short_long_common.csv"))
 
  # cf = "P19P31"
  total_long_cf_genes = c()
  total_short_cf_genes = c()
  total_common_cf_genes = c()
  
  for (slc in unique(cf_sl$classification)) {
    if (slc == "short") {
      
      for (cf in cf_sl[which(cf_sl$classification == slc),]$variable ) {
        
        num_P <- nchar(cf) - nchar(gsub("P", "", cf))
        
        if (num_P == 1) {
          tmp_single_genes = single_genes[which(single_genes$Pathway == cf),]$Genes 
          total_short_cf_genes = c(total_short_cf_genes,tmp_single_genes)
        } else if (num_P == 2) {
          
          tmp_link_genes = link_genes_filtered_df[which(link_genes_filtered_df$Pathway == cf),]$Genes 
          total_short_cf_genes = c(total_short_cf_genes,tmp_link_genes)
        }
      }
      
    } else if (slc == "long") {
     
      for (cf in cf_sl[which(cf_sl$classification == slc),]$variable) {
        
        num_P <- nchar(cf) - nchar(gsub("P", "", cf))
        
        if (num_P == 1) {
          tmp_single_genes = single_genes[which(single_genes$Pathway == cf),]$Genes 
          total_long_cf_genes = c(total_long_cf_genes,tmp_single_genes)
        } else if (num_P == 2) {
          
          tmp_link_genes = link_genes_filtered_df[which(link_genes_filtered_df$Pathway == cf),]$Genes 
          total_long_cf_genes = c(total_long_cf_genes,tmp_link_genes)
        }
      }
      
    } else if (slc == "common") {
      
     
      for (cf in cf_sl[which(cf_sl$classification == slc),]$variable) {
        
        num_P <- nchar(cf) - nchar(gsub("P", "", cf))
        
        if (num_P == 1) {
          tmp_single_genes = single_genes[which(single_genes$Pathway == cf),]$Genes 
          total_common_cf_genes = c(total_common_cf_genes,tmp_single_genes)
        } else if (num_P == 2) {
          
          tmp_link_genes = link_genes_filtered_df[which(link_genes_filtered_df$Pathway == cf),]$Genes 
          total_common_cf_genes = c(total_common_cf_genes,tmp_link_genes)
        }
      }
    }
  }
  
  total_long_cf_genes = unique(total_long_cf_genes)
  total_short_cf_genes = unique(total_short_cf_genes)
  total_common_cf_genes = unique(total_common_cf_genes)
  
  df_total_hallmark_genes
  
  # total_long_cf_genes와 gene_set의 길이를 구합니다.
  total_long_cf_genes_length <- length(total_long_cf_genes)
  total_short_cf_genes_length <- length(total_short_cf_genes)
  total_common_cf_genes_length <- length(total_common_cf_genes)
  
  # 결과를 저장할 빈 리스트 생성
  long_p_values <- list()
  short_p_values <- list()
  common_p_values <- list()
  
  jaccard_similarity <- function(set1, set2) {
    intersection <- length(intersect(set1, set2))
    union <- length(unique(c(set1, set2)))
    return(intersection / union)
  }
  
  jaccard_df = data.frame()
  
  # 각 gene_set 열에 대한 hypergeometric test를 수행
  for (col_name in colnames(df_total_hallmark_genes)) {
    
    tmp_hallmark_links = df_total_hallmark_genes[,col_name][!is.na(df_total_hallmark_genes[,col_name])]

    shared_long_genes_len <- length(intersect(total_long_cf_genes, tmp_hallmark_links))
    shared_short_genes_len <- length(intersect(total_short_cf_genes_length, tmp_hallmark_links))
    shared_common_genes_len <- length(intersect(total_common_cf_genes_length, tmp_hallmark_links))
    
    tmp_df = data.frame( long = jaccard_similarity( total_long_cf_genes, tmp_hallmark_links),
                         short = jaccard_similarity( total_short_cf_genes, tmp_hallmark_links),
                         common = jaccard_similarity( total_common_cf_genes, tmp_hallmark_links))
    rownames(tmp_df) = col_name
    jaccard_df = rbind(jaccard_df,tmp_df)
    
    ## 굳굳
    long_p_value <- phyper(shared_long_genes_len, 
                      length(tmp_hallmark_links) , 
                      sum(!is.na(unlist(df_total_hallmark_genes)))- length(tmp_hallmark_links), 
                      total_long_cf_genes_length, 
                      lower.tail = TRUE)
    
    short_p_value <- phyper(shared_short_genes_len, 
                           length(tmp_hallmark_links) , 
                           sum(!is.na(unlist(df_total_hallmark_genes)))- length(tmp_hallmark_links), 
                           total_short_cf_genes_length, 
                           lower.tail = TRUE)
    
    common_p_value <- phyper(shared_common_genes_len, 
                           length(tmp_hallmark_links) , 
                           sum(!is.na(unlist(df_total_hallmark_genes)))- length(tmp_hallmark_links), 
                           total_common_cf_genes_length, 
                           lower.tail = TRUE)
    
    
    long_p_values[[col_name]] <- long_p_value
    short_p_values[[col_name]] <- short_p_value
    common_p_values[[col_name]] <- common_p_value
  }
  
  # p-값을 출력합니다.
  # print(long_p_values)
  # print(short_p_values)
  # print(common_p_values)
  # 
  # fig for radar 
  df_long <- data.frame(
    long = unname(unlist(long_p_values)),
    row.names = names(long_p_values)
  )
  df_short <- data.frame(
    short = unname(unlist(short_p_values)),
    row.names = names(short_p_values)
  )
  df_common <- data.frame(
    common = unname(unlist(common_p_values)),
    row.names = names(common_p_values)
  )
  
  total_p_values = cbind(df_long, df_short, df_common)
  total_p_values_edit = -log2(total_p_values)
  
  total_score = total_p_values_edit * jaccard_df
  
  min_max_normalize <- function(x) {
    (x - min(x)) / (max(x) - min(x))
  }
  
  total_score = total_score %>%
    mutate(across(everything(), min_max_normalize))
  
  total_p_values_t = as.data.frame(t(total_score))
  total_p_values_t$group = rownames(total_p_values_t)
  rownames(total_p_values_t) = NULL
  total_p_values_t = total_p_values_t %>% 
    select(group , everything())
  
  total_p_values_t$group = factor(total_p_values_t$group, levels = c("long", "short", "common"))
  
  total_p_values_t[is.na(total_p_values_t)] <- 0
  
  library(ggradar)
  setwd("~/nas/04.Results/cancerhallmark/")
  
  png(filename = paste0(CancerType, "_critical_features_jaccard_multiple_hyper_radarplot.png"),
      width = 30, height = 30,  units = "cm" ,pointsize = 12,
      bg = "white", res = 1200, family = "")

  radar_plot = total_p_values_t %>% 
    ggradar(
      values.radar = c(0, 0.5, 1),
      
      background.circle.colour = "white",
      legend.position = "bottom",
      axis.label.offset = 1.05,
      # grid.label.size = 13,  # Affects the grid annotations (0%, 50%, etc.)
      axis.label.size = 2.9, # Afftects the names of the variables
      group.point.size = 3 ,  # Simply the size of the point 
      group.colours = c(short = "#FF5A5F", long = "#007A87", common = "#9CA299")
    )
  
  print(radar_plot)
  
  dev.off()
  
}
# "#FF5A5F", "#FFB400", "#007A87", "#8CE071", "#7B0051",
# "#00D1C1", "#FFAA91", "#B4A76C", "#9CA299", "#565A5C", "#00A04B", "#E54C20"


## 
length(intersect(total_long_cf_genes, tmp_hallmark_links)) / length(total_long_cf_genes)

total_genes_length = sum(!is.na(df_total_hallmark_genes$Genome_Instability_and_Mutation) , 
                         !is.na(df_total_hallmark_genes$Sustaining_Proliferative_Signaling), 
                         !is.na(df_total_hallmark_genes$Inducing_Angiogenesis), 
                         !is.na(df_total_hallmark_genes$Resisting_Cell_Death),
                         !is.na(df_total_hallmark_genes$Immune_modulation),
                         !is.na(df_total_hallmark_genes$Deregulating_Cellular_Energetics),
                         !is.na(df_total_hallmark_genes$Activating_Invasion_and_Metastasis)
)

ratio_cancerhallmark = data.frame()


for (cancer_hallmark in colnames(df_total_hallmark_genes)) {
  tmp_ratio = sum(!is.na(df_total_hallmark_genes[,cancer_hallmark])) / total_genes_length
  tmp_df = data.frame(ca  ㅌㅋ` ncer_hallmark = tmp_ratio)
  colnames(tmp_df) = cancer_hallmark
  
  if (ncol(ratio_cancerhallmark) == 0) {
    ratio_cancerhallmark = tmp_df
  } else {
    ratio_cancerhallmark = cbind(ratio_cancerhallmark, tmp_df)
  }
  
}

df_ratio_hallmark = data.frame(Hallmark = names(ratio_cancerhallmark), Ratio = t(ratio_cancerhallmark), row.names = NULL)
ggplot(data = df_ratio_hallmark, aes(x = Hallmark, y = Ratio, fill = Hallmark)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = round(Ratio, 3), vjust = -0.5), size = 3) +
  labs(x = "Cancer Hallmarks", y = "Ratio", title = "Ratio of Cancer Hallmarks") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = rainbow(nrow(df_ratio_hallmark))) +
  guides(fill = FALSE)

df_hallmark_ratio = data.frame()

for (col_name in colnames(df_total_hallmark_genes)) {
  
  tmp_hallmark_links = df_total_hallmark_genes[,col_name][!is.na(df_total_hallmark_genes[,col_name])]
  
  short_genes_ratio <- length(intersect(total_short_cf_genes, tmp_hallmark_links)) / length(total_short_cf_genes)
  long_genes_ratio <- length(intersect(total_long_cf_genes, tmp_hallmark_links)) / length(total_long_cf_genes)
  common_genes_ratio <- length(intersect(total_common_cf_genes, tmp_hallmark_links)) / length(total_common_cf_genes)
  
  tmp_dataframe = data.frame(short_genes_ratio,long_genes_ratio,common_genes_ratio , row.names = col_name)
  df_hallmark_ratio = rbind(df_hallmark_ratio, tmp_dataframe)
  
  
}

df_hallmark_ratio <- data.frame(
  Hallmark = rownames(df_hallmark_ratio),
  short_genes_ratio = df_hallmark_ratio$short_genes_ratio,
  long_genes_ratio = df_hallmark_ratio$long_genes_ratio,
  common_genes_ratio = df_hallmark_ratio$common_genes_ratio
)

# 데이터를 "long" 형식으로 변환
df_long <- reshape2::melt(df_hallmark_ratio, id.vars = "Hallmark")

# 그래프 그리기
ggplot(data = df_long, aes(x = Hallmark, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge" , width = 0.9) +
  labs(x = "Cancer Hallmarks", y = "Ratio", title = "Ratio of TCGA-LIHC") +
  geom_text(data = df_long, aes(label = round(value, 3), vjust = -0.5), position = position_dodge(width=0.9), size = 3) +
  # geom_text(data = df_long, aes(x = Hallmark, y = short_genes_ratio, label = round(short_genes_ratio, 3), vjust = -0.5), position = position_dodge(width=0.9),size = 3) +
  # geom_text(data = df_long, aes(x = Hallmark, y = long_genes_ratio, label = round(long_genes_ratio, 3), vjust = -0.5), position = position_dodge(width=0.9),size = 3) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
