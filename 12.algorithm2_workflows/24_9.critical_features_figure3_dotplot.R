library(tidyverse)
library(RColorBrewer)
library(magrittr)
library(stringr)
library(data.table)
library(pheatmap)
library(readxl)
library(dplyr)
library(readxl)
library(tidyverse)
library(UpSetR)
library(ComplexHeatmap)
library(svglite)

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

######

total_cal_critical_features = data.frame()

for (num_CancerType in Cancerlist) {
  
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  Cancernum = gsub('[.-]','',gsub('[a-zA-Z]','', num_CancerType))
  # call input
  
  cancer_sl = read.csv(paste0(filepath,"04.Results/short_long/ttest_common/",CancerType, "_critical_features_short_long_common.csv"))
  tmp_cancer = cancer_sl %>% dplyr::select(-X.1, -X) %>%
    mutate(cancertype = gsub("TCGA-","",CancerType)) %>% 
    dplyr::select(cancertype , everything())
  
  # tmp_cancer = tmp_cancer %>% filter(classification != "common")
  
  total_cal_critical_features = rbind(total_cal_critical_features, tmp_cancer)
}

count_table = total_cal_critical_features %>%
  count(variable) %>%
  arrange(desc(n))

filt_count_table = count_table %>% filter(n > 3)

total_cf = total_cal_critical_features %>% 
  mutate(group = case_when( variable %in% filt_count_table$variable ~ "shared",
                                                          !variable %in% filt_count_table$variable ~ "notspe",
                                                          .default = NA)) %>%
  filter(!is.na(group)) 

tmp_spe %>% filter(variable == "P10P12")  

tmp_shared = total_cf %>% filter(group == "shared")

tmp_spe = total_cf %>% filter(group != "shared")

cancer_with_more_than_10_variables <- tmp_spe %>%
  group_by(cancertype) %>%
  summarize(num_variables = n()) %>%
  filter(num_variables > 10) %>%
  pull(cancertype)

# variable의 개수가 10개를 넘는 cancertype 추출
cancer_with_more_than_10_variables <- count_by_cancer %>%
  filter(num_variables > 10) %>%
  pull(cancertype)

test = tmp_spe %>%
  group_by(cancertype) %>%
  mutate(rank = row_number(desc(minmax))) %>%
  filter(!(cancertype %in% cancer_with_more_than_10_variables & rank > 10)) %>%
  select(-rank) %>%
  ungroup() 


# minmax 값을 기준으로 내림차순 정렬 후 상위 10개만 남기기
tmp_spe_filtered <- tmp_spe %>%
  group_by(cancertype) %>%
  mutate(rank = row_number(desc(minmax))) %>%
  filter(!(cancertype %in% cancer_with_more_than_10_variables & rank > 10)) %>%
  select(-rank) %>%
  ungroup() %>% 
  mutate(group = cancertype)

under_three = count_table %>% filter(n <= 3 & n !=1)
# ttt = "P44"
for (ttt in under_three$variable) {
  
  tmp_tt = tmp_spe_filtered %>% filter(variable == ttt)
  if (nrow(tmp_tt) >=2) {
    first_cancer = tmp_tt[1,]
    second_cancer = tmp_tt[2,]
    first_cancer$group = second_cancer$cancertype
    second_cancer$group = first_cancer$cancertype
   
    tmp_spe_filtered = rbind(tmp_spe_filtered,first_cancer,second_cancer)
  }
  
}

tmp_shared = tmp_shared %>% select(cancertype , variable,relative_importance,minmax,pval,classification,group)

tmp_total_cf = rbind(tmp_shared , tmp_spe_filtered)

name_pathway = read_xlsx(paste0(ref_path,"pathway_name.xlsx"))

name_pathway$filtered_pathway_name = gsub("_signaling.*$|_pathway$","",name_pathway$Kegg_pathway_name)

# test_tmp = tmp_total_cf
# tmp_total_cf$pathway_name = ""
# 
# for (ev in tmp_total_cf$variable) {
#   count <- str_count(ev, "P")
#   
#   if (count == 1) {
#     tmp_name =  paste0( name_pathway[which(name_pathway$num_pathway == ev),]$filtered_pathway_name, ":",ev )
#     
#   } else {
#     first_path = paste0("P",str_split(ev, "P")[[1]][2])
#     second_path = paste0("P",str_split(ev, "P")[[1]][3])
#     
#     tmp_name =  paste(paste0(name_pathway[which(name_pathway$num_pathway == first_path),]$filtered_pathway_name,":", ev),
#                       paste0(name_pathway[which(name_pathway$num_pathway == second_path),]$filtered_pathway_name,":", ev),
#                       sep = "\n")
#   }
#   
#   tmp_total_cf[which(tmp_total_cf$variable == ev),"pathway_name"] = tmp_name
#   
# }

tmp_total_cf$pathway_name = ""

for (ev in tmp_total_cf$variable) {
  count <- str_count(ev, "P")
  
  if (count == 1) {
    tmp_name =  paste0( name_pathway[which(name_pathway$num_pathway == ev),]$Shorter_name , " : ",ev )
    
  } else {
    first_path = paste0("P",str_split(ev, "P")[[1]][2])
    second_path = paste0("P",str_split(ev, "P")[[1]][3])
    
    tmp_name =  paste0(paste(name_pathway[which(name_pathway$num_pathway == first_path),]$Shorter_name , " / " ,
                             name_pathway[which(name_pathway$num_pathway == second_path),]$Shorter_name ), " : ",ev)
  }
  
  tmp_total_cf[which(tmp_total_cf$variable == ev),"pathway_name"] = tmp_name
  
}

tmp_total_cf = tmp_total_cf %>% filter(minmax !=0)


group_order = c("shared","BLCA","OV","UCEC","LIHC","STAD","LUAD","CESC","LUSC","BRCA","LGG")
tmp_total_cf$group <- factor(tmp_total_cf$group, levels = group_order)
tmp_total_cf$cancertype <- factor(tmp_total_cf$cancertype, levels = c("BLCA","OV","UCEC","LIHC","STAD","LUAD","CESC","LUSC","BRCA","LGG"))

tmp_total_cf <- tmp_total_cf %>%
  group_by(pathway_name) %>%
  mutate(cancertype_freq = n()) %>%
  ungroup() %>%
  arrange(desc(cancertype_freq))

tmp_total_cf = tmp_total_cf %>%
  group_by(pathway_name) %>%
  mutate(total_minmax = sum(minmax)) %>%
  ungroup()

five_order = tmp_total_cf %>% filter(cancertype_freq == 5) %>% arrange(desc(total_minmax))
four_order = tmp_total_cf %>% filter(cancertype_freq == 4) %>% arrange(desc(total_minmax))
three_order = tmp_total_cf %>% filter(cancertype_freq == 3) %>% arrange(desc(total_minmax))
two_order = tmp_total_cf %>% filter(cancertype_freq == 2) %>% arrange(desc(total_minmax))
one_order = tmp_total_cf %>% filter(cancertype_freq == 1) %>% arrange(desc(total_minmax))

path_order = c(unique(five_order$pathway_name), 
               unique(four_order$pathway_name),
               unique(three_order$pathway_name),
               unique(two_order$pathway_name),
               unique(one_order$pathway_name))

tmp_total_cf$pathway_name = factor(tmp_total_cf$pathway_name, levels = path_order)


dotplot = ggplot(tmp_total_cf, aes(x = pathway_name , y = cancertype , size = minmax, color = minmax)) +
  geom_point() +
  # scale_size_continuous(range = c(1, 10)) +
  scale_color_gradient(low = "white", high = "red") +
  facet_grid(. ~ group,
             scales = "free", 
             space = "free_x") +
  theme_classic() +
  
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face='bold', family = "arial" ),
        axis.text.y = element_text( face='bold', family = "arial"),
        panel.spacing = unit(0, "lines"),
        # panel.border = element_blank(),
        # panel.grid.major = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        strip.background = element_rect(color = "black", size = 1))
library(ggplot2)

ggplot2::ggsave(file="~/nas/04.Results/critical_features/figure3A_dotplot.svg", plot=dotplot, width=25, height=10)







#### 
# transfer_m_df = data.frame()
# 
# for (features in unique(filt_count_table$variable)) {
#   tmp_df = total_cal_critical_features %>% filter(variable == features)
#   
#   tmp2_df = as.data.frame(t(data.frame(as.numeric(gsub(".*\\.TCGA-", "", Cancerlist) %in% unique(tmp_df$cancertype)))))
#   colnames(tmp2_df) = gsub(".*\\.TCGA-", "", Cancerlist)
#   tmp2_df$features = features
#   rownames(tmp2_df) = NULL
#   tmp2_df = tmp2_df %>% dplyr::select(features, everything())
#   
#   transfer_m_df = rbind(transfer_m_df,tmp2_df)
# }
# 
# backup_m_df = transfer_m_df
# tmp_features = transfer_m_df$features
# transfer_m_df = transfer_m_df %>% dplyr::select(-features)
# 
# tmp_col = colnames(transfer_m_df)
# 
# colnames(transfer_m_df) = NULL
# 
# t_transfer_m_df = as.data.frame(t(transfer_m_df))
# 
# rownames(t_transfer_m_df) = tmp_col
# colnames(t_transfer_m_df) = tmp_features
# 
# m3= make_comb_mat(t_transfer_m_df)
# 
# cs_tmp = rowSums(t_transfer_m_df)
# cs_tmp = sort(cs_tmp , decreasing = T)
# name_cancer = cs_tmp
# names(cs_tmp) = names(cs)
# cs = cs_tmp
# 
# ss = set_size(m3)
# # cs = comb_size(m3)
# 
# col_top = c("#7E6148FF",
#             "#E64B35FF",
#             "#3C5488FF",
#             "#8491B4FF",
#             "#91D1C2FF",
#             "#B09C85FF",
#             "#F39B7FFF",
#             "#DC0000FF",
#             "#4DBBD5FF",
#             "#00A087FF")
# 
# names(col_top) = names(comb_degree(m3))
# 
# 
# upset_m3_df = data.frame()
# 
# for (upset_m3 in set_name(m3)) {
#   count <- str_count(upset_m3, "P")
#   
#   if (count == 1) {
#     tmp_pathway = data.frame( features = upset_m3,
#                               # pathway = upset_m3,
#                               pathway_name = name_pathway[which(name_pathway$num_pathway == upset_m3),]$Kegg_pathway_name
#     )
#     
#   } else {
#     first_path = paste0("P",str_split(upset_m3, "P")[[1]][2])
#     second_path = paste0("P",str_split(upset_m3, "P")[[1]][3])
#     
#     tmp_pathway_link = data.frame( features = upset_m3, 
#                                    pathway_name = paste(name_pathway[which(name_pathway$num_pathway == first_path),]$Kegg_pathway_name,
#                                                         name_pathway[which(name_pathway$num_pathway == second_path),]$Kegg_pathway_name,
#                                                         sep = "\n"))
#     
#     tmp_pathway = tmp_pathway_link
#   }
#   
#   upset_m3_df = rbind(upset_m3_df, tmp_pathway)
#   
# }
# 
# ht = UpSet(m3, 
#            pt_size = unit(10, "mm"), lwd = 7,
#            set_order = order(ss),
#            comb_order = order(-comb_degree(m3), -cs),
#            top_annotation = HeatmapAnnotation(
#              "Cancer counts" = anno_barplot(cs, 
#                                       ylim = c(0, max(cs)*1.1),
#                                       border = FALSE, 
#                                       gp = gpar(fill = col_top), 
#                                       height = unit(4, "cm")
#              ), 
#              annotation_name_side = "left", 
#              annotation_name_rot = 90),
#            right_annotation = rowAnnotation(
#              "counts" = anno_barplot(-ss, 
#                                      baseline = 0,
#                                      axis_param = list(
#                                        at = c(-5, -4, -3, -2, -1, 0),
#                                        labels = c(5,4,3,2,1,0),
#                                        direction = "reverse",
#                                        labels_rot = 0),
#                                      border = FALSE, 
#                                      gp = gpar(fill = "black"), 
#                                      width = unit(4, "cm")
#              ),
#              set_name = anno_text(upset_m3_df$pathway_name,
#                                   location = 1,
#                                   just = "right",
#                                   width = max_text_width(upset_m3_df$pathway_name) + unit(6, "mm"))
#            ), 
#            left_annotation = NULL,
#            show_row_names = T,
#            
#            )
# 
# svglite("~/nas/04.Results/critical_features/critical_features_upset.svg", width = 10, height = 10)
# 
# ht = draw(ht)
# od = column_order(ht)
# decorate_annotation("Cancer counts", {
#   grid.text(names(name_cancer[od]), x = seq_along(name_cancer), y = unit(name_cancer[od], "native") + unit(2, "pt"), 
#             default.units = "native", just = c( "bottom"), 
#             gp = gpar(fontsize = 10, col = "#404040"), rot = 0, )
# })
# 
# dev.off()
# 
# ######
# 
# 
# total_circos_sl_df = data.frame()
# for (i in 1:nrow(total_comb)) {
#   first_cancer = total_comb[i,1]
#   second_cancer = total_comb[i,2]
#   # first_cancer = "BLCA"
#   # second_cancer = "UCEC"
#   if (length(intersect(total_features[[first_cancer]], total_features[[second_cancer]])) != 0 ) {
#     
#     shared_features = intersect(total_features[[first_cancer]], total_features[[second_cancer]])
#     
#     first_cancer_sl = read.csv(paste0(filepath,"04.Results/short_long/ttest_common/",gsub('[.]','',gsub('\\d','',Cancerlist[grep(first_cancer, Cancerlist)])), "_critical_features_short_long_common.csv"))
#     second_cancer_sl = read.csv(paste0(filepath,"04.Results/short_long/ttest_common/",gsub('[.]','',gsub('\\d','',Cancerlist[grep(second_cancer, Cancerlist)])), "_critical_features_short_long_common.csv"))
#     
#     first_cl = first_cancer_sl %>% 
#       subset(., variable %in% shared_features) 
#     
#     second_cl = second_cancer_sl %>% 
#       subset(., variable %in% shared_features)
#     
#     # second_cl = second_cl %>% arrange(factor(variable, levels = first_cl$variable))
#     
#     tmp_equal_enriched = merge(first_cl , second_cl , by = "variable", all = FALSE )
#     
#     equal_enriched = tmp_equal_enriched %>% 
#       dplyr::select(-X.1.x, -X.x, -X.1.y,-X.y) %>% 
#       mutate(total_minmax = minmax.x * minmax.y) %>% 
#       filter(classification.x == classification.y) %>%
#       arrange(desc(total_minmax)) 
#     
#     equal_enriched = equal_enriched %>%
#       mutate( classification = classification.x) %>%
#       dplyr::select(-classification.x ,- classification.y)
#     
#     equal_enriched = equal_enriched %>% filter(classification != "common")
#     
#     if (nrow(equal_enriched) != 0) {
#       tmp_first_cir_df = equal_enriched %>% mutate(cancertype = first_cancer) %>%
#         dplyr::select(cancertype, variable , total_minmax , classification)
#       tmp_second_cir_df = equal_enriched %>% mutate(cancertype = second_cancer) %>%
#         dplyr::select(cancertype, variable , total_minmax , classification)
#     } else {
#       next
#     }
#     tmp_cir_df = rbind(tmp_first_cir_df, tmp_second_cir_df)
#     
#     # tmp_cir_df_filt = tmp_cir_df %>% 
#     #   select(variable , cancertype, relative_importance,minmax,pval,classification) %>%
#     #   mutate(weight_filt = weight * minmax * 10^4)
#     
#     tmp_cir_df_filt = tmp_cir_df %>% 
#       dplyr::select(variable , cancertype,classification)
#     
#   } else {
#     next
#   }
#   total_circos_sl_df = rbind(total_circos_sl_df, tmp_cir_df_filt)
#   
# }
# 
# top_enriched = head(sort(table(total_circos_sl_df$variable) / 2 , decreasing = T), n = 12)
# 
# top_enriched_df = total_circos_sl_df %>% filter(variable %in% names(top_enriched))
# 
# transfer_df = data.frame()
# 
# for (features in unique(top3_enriched_df$variable)) {
#   tmp_df = top3_enriched_df %>% filter(variable == features)
#   
#   tmp2_df = as.data.frame(t(data.frame(as.numeric(gsub(".*\\.TCGA-", "", Cancerlist) %in% unique(tmp_df$cancertype)))))
#   colnames(tmp2_df) = gsub(".*\\.TCGA-", "", Cancerlist)
#   tmp2_df$features = features
#   rownames(tmp2_df) = NULL
#   tmp2_df = tmp2_df %>% dplyr::select(features, everything())
#   
# 
#   transfer_df = rbind(transfer_df,tmp2_df)
# }
# backup_df = transfer_df
# tmp_features = transfer_df$features
# transfer_df = transfer_df %>% dplyr::select(-features)
# 
# tmp_col = colnames(transfer_df)
# 
# colnames(transfer_df) = NULL
# 
# t_transfer_df = as.data.frame(t(transfer_df))
# 
# rownames(t_transfer_df) = tmp_col
# colnames(t_transfer_df) = tmp_features
# 
# m2= make_comb_mat(t_transfer_df)
# 
# # total_shared_features_filt %>% filter(shared_features == "P54")
# 
# cs_slc_tmp = rowSums(t_transfer_df)
# cs_slc_tmp = sort(cs_slc_tmp , decreasing = T)
# name_cancer = cs_slc_tmp
# 
# cs = comb_size(m2)
# 
# names(cs_slc_tmp) = names(cs)
# cs_slc_tmp = cs_slc_tmp[!is.na(names(cs_slc_tmp))]
# 
# cs_slc = cs_slc_tmp
# ss = set_size(m2)
# 
# col_slc_top = c("#72BC6C",
#                 "#72BC6C",
#                 "#72BC6C",
#                 "#D3DFE5",
#                 "#72BC6C",
#                 "#C0392B",
#                 "#C0392B")
#                       
# 
# names(col_top) = names(comb_degree(m2))
# 
# name_pathway = read_xlsx(paste0(ref_path,"pathway_name.xlsx"))
# 
# upset_m2_df = data.frame()
# 
# for (upset_m2 in set_name(m2)) {
#   count <- str_count(upset_m2, "P")
#   
#   if (count == 1) {
#     tmp_pathway = data.frame( features = upset_m2,
#                               # pathway = upset_m3,
#                               pathway_name = name_pathway[which(name_pathway$num_pathway == upset_m2),]$Kegg_pathway_name
#     )
#     
#   } else {
#     first_path = paste0("P",str_split(upset_m2, "P")[[1]][2])
#     second_path = paste0("P",str_split(upset_m2, "P")[[1]][3])
#     
#     tmp_pathway_link = data.frame( features = upset_m2, 
#                                    pathway_name = paste(name_pathway[which(name_pathway$num_pathway == first_path),]$Kegg_pathway_name,
#                                                         name_pathway[which(name_pathway$num_pathway == second_path),]$Kegg_pathway_name,
#                                                         sep = "\n"))
#     
#     tmp_pathway = tmp_pathway_link
#   }
#   
#   upset_m2_df = rbind(upset_m2_df, tmp_pathway)
#   
# }
# 
# 
# ht_slc = UpSet(m2, 
#            pt_size = unit(10, "mm"), lwd = 7,
#            set_order = order(ss),
#            comb_order = order(-comb_degree(m2), -cs_slc),
#            top_annotation = HeatmapAnnotation(
#              "Cancer counts" = anno_barplot(cs_slc, 
#                                             ylim = c(0, max(cs_slc)*1.1),
#                                             border = FALSE, 
#                                             gp = gpar(fill = col_slc_top), 
#                                             height = unit(4, "cm")
#              ), 
#              annotation_name_side = "left", 
#              annotation_name_rot = 90),
#            right_annotation = rowAnnotation(
#              "counts" = anno_barplot(-ss, 
#                                      baseline = 0,
#                                      axis_param = list(
#                                        at = c(-5, -4, -3, -2, -1, 0),
#                                        labels = c(5,4,3,2,1,0),
#                                        direction = "reverse",
#                                        labels_rot = 0),
#                                      border = FALSE, 
#                                      gp = gpar(fill = "black"), 
#                                      width = unit(4, "cm")
#              ),
#              set_name = anno_text(upset_m2_df$pathway_name,
#                                   location = 1,
#                                   just = "right",
#                                   width = max_text_width(upset_m2_df$pathway_name) + unit(6, "mm"))
#            ), 
#            left_annotation = NULL,
#            show_row_names = T,
#            
# )
# 
# 
# svglite("~/nas/04.Results/critical_features/cf_enriched_slc_upset.svg", width = 10, height = 10)
# ht_slc = draw(ht_slc)
# od_slc = column_order(ht_slc)
# 
# decorate_annotation("Cancer counts", {
#   grid.text(names(name_cancer[od_slc]), x = seq_along(name_cancer[od_slc]), y = unit(name_cancer[od_slc], "native") + unit(2, "pt"), 
#             default.units = "native", just = c( "bottom"), 
#             gp = gpar(fontsize = 10, col = "#404040"), rot = 0, )
# })
# 
# dev.off()
# 
