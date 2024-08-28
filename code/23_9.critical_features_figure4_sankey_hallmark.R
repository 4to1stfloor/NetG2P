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
library(circlize)

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
# cancerhallmark2 = read_xlsx(paste0(ref_path, "kegg_gene_set_w_cancer_hallmarks_edit.xlsx"))
cancerhallmark = read.csv(paste0(ref_path, "Zhang2020_CHG.csv"), sep = "\t")  
colnames(cancerhallmark)= gsub("\\.","_",colnames(cancerhallmark))
# 
jaccard_similarity <- function(set1, set2) {
  intersection <- length(intersect(set1, set2))
  union <- length(unique(c(set1, set2)))
  return(intersection / union)
}
min_max_normalize <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}
# for fic
fig_path = paste0(filepath,"04.Results/sankey_cancerhallmark/")
if(!dir.exists(fig_path)){
  dir.create(fig_path)
  print(paste0("Created folder: ", fig_path))
} else {
  print(paste0("Folder already exists: ", fig_path))
}
setwd(fig_path)

total_hallmark_genes_df = reshape2::melt(cancerhallmark, id.vars = NULL)
total_hallmark_genes_df$value = toupper(total_hallmark_genes_df$value)
colnames(total_hallmark_genes_df) = c("cancerhallmark_name" , "genes")

library(clusterProfiler)
# setwd("~/nas/99.reference/")
# write.xlsx(df_total_hallmark_genes, "gene_set_w_cancer_hallmarks_edit.xlsx")
# num_CancerType = "32.TCGA-UCEC"
Cancerlist = c("10.TCGA-BLCA","24.TCGA-OV" , "32.TCGA-UCEC")
total_score_df = data.frame()
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
  # total_common_cf_genes = c()
 
  for (slc in unique(cf_sl$classification)) {
    if (slc == "short") {
      # cf = "P20P28"
      for (cf in cf_sl[which(cf_sl$classification == slc),]$variable ) {
        
        num_P <- nchar(cf) - nchar(gsub("P", "", cf))
        
        if (num_P == 1) {
          tmp_single_genes = single_genes[which(single_genes$Pathway == cf),]$Genes 
          # print(enricher(gene = tmp_link_genes, TERM2GENE = total_hallmark_genes_df, minGSSize=1)@result)
          total_short_cf_genes = c(total_short_cf_genes,tmp_single_genes)
        } else if (num_P == 2) {
          
          tmp_link_genes = link_genes_filtered_df[which(link_genes_filtered_df$Pathway == cf),]$Genes 
          # print(enricher(gene = tmp_link_genes, TERM2GENE = total_hallmark_genes_df, minGSSize=1)@result)
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
      next
    }
  }

  # tmp_enrichr_df = data.frame(matrix(nrow = length(unique(total_hallmark_genes_df$cancerhallmark_name)), ncol = 3))
  tmp_enrichr_df = data.frame(matrix(nrow = length(unique(total_hallmark_genes_df$cancerhallmark_name)), ncol = 2))
  rownames(tmp_enrichr_df) = unique(total_hallmark_genes_df$cancerhallmark_name)
  # colnames(tmp_enrichr_df) = c("long", "short" , "common")
  colnames(tmp_enrichr_df) = c("long", "short")
  ## 굳굳
  if (!is.null(enricher(gene = total_long_cf_genes, TERM2GENE = total_hallmark_genes_df, minGSSize=1))) {
    long_res <- enricher(gene = total_long_cf_genes, TERM2GENE = total_hallmark_genes_df, minGSSize=1)@result
    long_res_filt = long_res[long_res$qvalue < 0.05,]
    
    if (!has_element(rownames(long_res_filt),"NA")) {
      for (des in long_res_filt$Description) {
        tmp_enrichr_df[des,"long"] = long_res_filt[which(long_res_filt$Description == des),]$qvalue
      }
    }
    
  } 
  
  if (!is.null(enricher(gene = total_short_cf_genes, TERM2GENE = total_hallmark_genes_df, minGSSize=1))) {
    short_res <- enricher(gene = total_short_cf_genes, TERM2GENE = total_hallmark_genes_df, minGSSize=1)@result
    short_res_filt = short_res[short_res$qvalue < 0.05,]
    
    if (!has_element(rownames(short_res_filt),"NA")) {
      for (des in short_res_filt$Description) {
        tmp_enrichr_df[des,"short"] = short_res_filt[which(short_res_filt$Description == des),]$qvalue
      }
    }
   
  } 
  
  # if (!is.null(enricher(gene = total_common_cf_genes, TERM2GENE = total_hallmark_genes_df, minGSSize=1))) {
  #   common_res <- enricher(gene = total_common_cf_genes, TERM2GENE = total_hallmark_genes_df, minGSSize=1)@result
  #   common_res_filt = common_res[common_res$qvalue < 0.05,]
  # 
  #   if (!has_element(rownames(common_res_filt),"NA")) {
  #     for (des in common_res_filt$Description) {
  #       tmp_enrichr_df[des,"common"] = common_res_filt[which(common_res_filt$Description == des),]$qvalue
  #     }
  #   }
  #  
  #   
  # } 

  tmp_enrichr_df[is.na(tmp_enrichr_df)] = 1

  total_p_values_edit = -log2(tmp_enrichr_df)
  
  jaccard_df = data.frame()
  
  # 각 gene_set 열에 대한 hypergeometric test를 수행
  for (col_name in unique(total_hallmark_genes_df$cancerhallmark_name)) {
   
    tmp_hallmark_links = total_hallmark_genes_df[which(total_hallmark_genes_df$cancerhallmark_name == col_name),]$genes

    # tmp_df = data.frame( long = jaccard_similarity( total_long_cf_genes, tmp_hallmark_links),
    #                      short = jaccard_similarity( total_short_cf_genes, tmp_hallmark_links),
    #                      common = jaccard_similarity( total_common_cf_genes, tmp_hallmark_links))
    tmp_df = data.frame( long = jaccard_similarity( total_long_cf_genes, tmp_hallmark_links),
                         short = jaccard_similarity( total_short_cf_genes, tmp_hallmark_links))
    
    rownames(tmp_df) = col_name
    jaccard_df = rbind(jaccard_df,tmp_df)
  }
  
  total_score = total_p_values_edit + jaccard_df
  
  total_score = total_score %>%
    mutate(across(everything(), min_max_normalize))
  
  total_score$hallmark = rownames(total_score)
  total_score$cancertype = gsub('TCGA-','',CancerType)
  total_score = total_score %>% select(hallmark, cancertype , everything())
  rownames(total_score) = NULL
  
  tmp_score_gather <- tidyr::gather(total_score, key="features", value="value", -hallmark:-cancertype)
  total_score_df = rbind(total_score_df,tmp_score_gather)
  
}

total_score_df[is.na(total_score_df$value),]$value = 0
circlize_tbl <- total_score_df 
circlize_tbl <- circlize_tbl %>%
  filter(cancertype != "COADREAD" &
           cancertype != "KIDNEY")

### ordering
# sender = cancer type
# reciver = features
# sender_order = circlize_tbl$cancertype %>% unique() # %>% sort()

sender_order = circlize_tbl %>% group_by(cancertype) %>%
  summarise(total_value = sum(value)) %>%
  arrange(total_value)
sender_order = sender_order$cancertype

# receiver_order = unique(circlize_tbl$hallmark)
receiver_order = circlize_tbl %>% group_by(hallmark) %>%
  summarise(total_value = sum(value)) %>%
  arrange(desc(total_value))

receiver_order = receiver_order$hallmark

order = c(receiver_order , sender_order)

### color
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
  "CESC" = "#DC0000FF"
)

sender_color <- set_names(color_map, c("UCEC" , 
                                       "BRCA",    
                                       "LGG" ,
                                       "LUSC" ,
                                       "OV",
                                       "LUAD",
                                       "LIHC",
                                       "STAD",
                                       "BLCA",
                                       "CESC"))

sender_color = sender_color[sender_order]

receiver_color <- set_names(rep("black", nrow(circlize_tbl)), circlize_tbl$hallmark)
receiver_color = receiver_color[receiver_order]

grid_col = c( receiver_color,sender_color)

slc_map = c("long" = "#72BC6C",
            "short" = "#C0392B",
            "common" = "#D3DFE5")
circlize_tbl = circlize_tbl %>% mutate(edges_color = case_when( features == "long" ~ "#72BC6C",
                                                                features == "short" ~ "#C0392B",
                                                                features == "common" ~ "#D3DFE5",
                                                                .default = NA))

svg("figure3G.svg")

circos.clear()
# circos.par(start.degree = 100,
#            gap.degree=0.5,
#            gap.after = c("UCEC" = 10, "P15" = 10))

circos.par(start.degree = 85,
           gap.degree= 0.5,
           gap.after = c( rep(2, length(receiver_order)-1), 10,rep(2, length(sender_order)-1), 10))

chordDiagram(circlize_tbl %>% 
               select(cancertype  , hallmark  , value  ), directional = 1, 
             order=order,
             grid.col = grid_col,
             col = circlize_tbl$edges_color,
             transparency = 0.6, 
             diffHeight = 0, 
             target.prop.height = mm_h(4),
             direction.type = c("diffHeight", "arrows"),
             link.arr.type = "big.arrow",
             # big.gap = 3,
             # link.arr.type = "curved",
             annotationTrack = "grid", 
             preAllocateTracks = list(track.height = 0.15))

circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 1)
}, bg.border = NA) #

title("Shared_critical_features")

# Close the graphics device
dev.off() 

##### supple

circlize_tbl_total = circlize_tbl %>% mutate(edges_color = color_map[cancertype])

svg("supple_figure_circos.svg")

circos.clear()
# circos.par(start.degree = 100,
#            gap.degree=0.5,
#            gap.after = c("UCEC" = 10, "P15" = 10))

circos.par(start.degree = 85,
           gap.degree= 0.5,
           gap.after = c( rep(2, length(receiver_order)-1), 10,rep(2, length(sender_order)-1), 10))

chordDiagram(circlize_tbl_total %>% 
               select(cancertype  , hallmark  , value  ), directional = 1, 
             order=order,
             grid.col = grid_col,
             col = circlize_tbl_total$edges_color,
             transparency = 0.6, 
             diffHeight = 0, 
             target.prop.height = mm_h(4),
             direction.type = c("diffHeight", "arrows"),
             link.arr.type = "big.arrow",
             # big.gap = 3,
             # link.arr.type = "curved",
             annotationTrack = "grid", 
             preAllocateTracks = list(track.height = 0.15))

circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 1)
}, bg.border = NA) #

# title("Shared_critical_features")

# Close the graphics device
dev.off() 




