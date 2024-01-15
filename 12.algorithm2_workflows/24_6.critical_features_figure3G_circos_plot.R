library(circlize)
library(tidyverse)
library(RColorBrewer)
library(networkD3)
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

Cancerlist = dir(paste0(filepath, "/00.data/filtered_TCGA/"))

# num_CancerType = "19.TCGA-LIHC"

####
Cancerlist = Cancerlist[c(-11,-12)]
total_features = list()

# for all
for (num_CancerType in Cancerlist) {
  
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  Cancernum = gsub('[.-]','',gsub('[a-zA-Z]','', num_CancerType))
  # call input
  
  # cancer_bf = read.csv(paste0(filepath,"04.Results/bestfeatures/",CancerType, "_critical_features.csv"))
  cancer_bf_sl = read.csv(paste0(filepath,"04.Results/short_long/ttest_common/",CancerType, "_critical_features_short_long_common.csv"))
  
  # cancer_bf_cut = cancer_bf[1:surv_total_results[which(surv_total_results$CancerType == CancerType),]$num_of_features,]
  Cancername = gsub('TCGA-','', CancerType)
  total_features[[Cancername]] = cancer_bf_sl$variable
  
}

total_comb = gsub('.TCGA-','',gsub('\\d','', t(combn(Cancerlist, 2))))
total_network_sl_df = data.frame()

for (i in 1:nrow(total_comb)) {
  first_cancer = total_comb[i,1]
  second_cancer = total_comb[i,2]
  # first_cancer = "CESC"
  # second_cancer = "UCEC"
  if (length(intersect(total_features[[first_cancer]], total_features[[second_cancer]])) != 0 ) {
    
    shared_features = intersect(total_features[[first_cancer]], total_features[[second_cancer]])
    
    first_cancer_sl = read.csv(paste0(filepath,"04.Results/short_long/ttest_common/",gsub('[.]','',gsub('\\d','',Cancerlist[grep(first_cancer, Cancerlist)])), "_critical_features_short_long_common.csv"))
    second_cancer_sl = read.csv(paste0(filepath,"04.Results/short_long/ttest_common/",gsub('[.]','',gsub('\\d','',Cancerlist[grep(second_cancer, Cancerlist)])), "_critical_features_short_long_common.csv"))
    
    first_cl = first_cancer_sl %>% 
      subset(., variable %in% shared_features) 
    
    second_cl = second_cancer_sl %>% 
      subset(., variable %in% shared_features)
    
    second_cl = second_cl %>% arrange(factor(variable, levels = first_cl$variable))
    
    if (all.equal(first_cl$variable ,second_cl$variable)) {
      equal_enriched = first_cl[which(first_cl$classification == second_cl$classification),]
    }
    
    if (has_element(equal_enriched$classification , "common")) {
      equal_enriched = equal_enriched %>% filter( classification != 'common')
    }
    
    if (nrow(equal_enriched) > 5) {
      equal_enriched = equal_enriched %>% slice(1:5)
    }
    # print(unique(equal_enriched$classification))
    if (length(unique(equal_enriched$classification)) == 1 ) {
      tmp_prop = prop.table(table(equal_enriched$classification))
      # print(tmp_prop)
      if (sum(tmp_prop > 0.8) == 1) {
        
        equal_enriched = equal_enriched %>% filter( classification == names(tmp_prop[tmp_prop > 0.8]))
        
      } else {
        equal_enriched = data.frame()
      }
    } else if (length(unique(equal_enriched$classification)) == 2) {
      # print(first_cancer)
      # print(second_cancer)
      equal_enriched$ratio = tmp_prop
      equal_enriched$classification_backup = equal_enriched$classification
      equal_enriched$classification = "mixed"
    } else {
      equal_enriched = equal_enriched
    }
    
    # print(unique(equal_enriched$classification))
    
    if (nrow(equal_enriched) != 0) {
      tmp_first_cir_df = equal_enriched %>% mutate(cancertype = first_cancer) %>%
        select(cancertype, variable , relative_importance, minmax, pval, classification)
      tmp_second_cir_df = equal_enriched %>% mutate(cancertype = second_cancer) %>%
        select(cancertype, variable , relative_importance, minmax, pval, classification)
    } else {
      next
    }
    tmp_cir_df = rbind(tmp_first_cir_df, tmp_second_cir_df)
    

    tmp_cir_df[which(tmp_cir_df$cancertype == first_cancer),]$minmax = first_cl[which(first_cl$variable %in% tmp_cir_df$variable),]$minmax
    tmp_cir_df[which(tmp_cir_df$cancertype == second_cancer),]$minmax = second_cl[which(second_cl$variable %in% tmp_cir_df$variable),]$minmax
    
    
    weight = length(equal_enriched$variable) / (length(total_features[[first_cancer]]) + 
                                         length(total_features[[second_cancer]]) - 
                                         length(intersect(total_features[[first_cancer]], total_features[[second_cancer]])))
    
    # tmp_cir_df_filt = tmp_cir_df %>% 
    #   select(variable , cancertype, relative_importance,minmax,pval,classification) %>%
    #   mutate(weight_filt = weight * minmax * 10^4)

    tmp_cir_df_filt = tmp_cir_df %>% 
      select(variable , cancertype, relative_importance,minmax,pval,classification) %>%
      mutate(weight_filt =1)

  } else {
    next
  }
  
  total_network_sl_df = rbind(total_network_sl_df, tmp_cir_df_filt)
  
}

### filtering featrures at inner each cancer 

total_network_filt_df = data.frame()

for (type in unique(c(total_comb[,1],total_comb[,2]))) {
  tmp_network_df = total_network_sl_df %>% 
    filter( cancertype == type) %>%
    group_by(variable) %>%
    mutate(cancer_num = n()) %>%
    filter(row_number() == 1) %>%
    ungroup() 
  total_network_filt_df = rbind(total_network_filt_df, tmp_network_df)
}

### filtering features by duplicated features  

total_network_filt_fi_df = data.frame()
for (vari in unique(total_network_sl_df$variable)) {
  tmp_vari_ch = total_network_filt_df  %>%
    filter( variable == vari) %>%
    mutate(feature_num = n()) 
  
  total_network_filt_fi_df = rbind(total_network_filt_fi_df, tmp_vari_ch)
  
}

circlize_tbl <- total_network_filt_fi_df %>% 
  arrange(desc(feature_num)) 

### ordering
# sender = cancer type
# reciver = features
sender_order = circlize_tbl$cancertype %>% unique() # %>% sort()
sender_order = names(sort(table(circlize_tbl$cancertype)))

receiver_order = unique(circlize_tbl$variable)

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

receiver_color <- set_names(rep("black", nrow(circlize_tbl)), circlize_tbl$variable)
receiver_color = receiver_color[receiver_order]

grid_col = c( receiver_color,sender_color)

svg("figure3G.svg")

circos.clear()
# circos.par(start.degree = 100,
#            gap.degree=0.5,
#            gap.after = c("UCEC" = 10, "P15" = 10))

circos.par(start.degree = 90,
           gap.degree= 0.5,
           gap.after = c( rep(2, length(receiver_order)-1), 10,rep(2, length(sender_order)-1), 10))

chordDiagram(circlize_tbl %>% 
               select(cancertype  , variable  , weight_filt  ), directional = 1, 
             order=order,
             grid.col = grid_col,
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
