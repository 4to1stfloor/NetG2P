# library(h2o)
filepath = "/home/seokwon/nas/"
ref_path = paste0(filepath, "99.reference/")
Cancerlist2 = dir(paste0(filepath, "/00.data/total/"))
Cancerlist = dir(paste0(filepath, "/00.data/filtered_TCGA/"))

Cancerlist2 = Cancerlist2[c(-1,-3,-6,-15,-16)]

# folder_name = "h2o_bias_pval_exp_dual_cut_50"
Cancerlist_total = sort(c(Cancerlist, Cancerlist2))
Cancerlist_total = Cancerlist_total[c(-1,-3,-7,-20,-21)]
# set.seed(13524)
# num_CancerType = "04.TCGA-CESC"
total_pat = data.frame()
# num_CancerType = "02.TCGA-UCS"

for (num_CancerType in Cancerlist) {

  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  
  data_tc_link = readRDS(file = paste0(main.path_tc,"/",CancerType,"_pathwaylink_all_0or1.rds"))
  
  if ("Not Reported" %in% data_tc_link$vitalstatus ) {
    if (length(which(data_tc_link$vitalstatus == "Not Reported")) == 1) {
      data_tc_link = data_tc_link[-which(data_tc_link$vitalstatus == "Not Reported"),]
    } else {
      while ("Not Reported" %in% data_tc_link$vitalstatus) {
        data_tc_link = data_tc_link[-which(data_tc_link$vitalstatus == "Not Reported")[1],]
      }
      
    }
    
  }
  tmp_summary = data.frame(cancer_type = CancerType , total_length = nrow(data_tc_link) , bias = sum(data_tc_link$vitalstatus == "Alive" ) / length(data_tc_link$vitalstatus))
  total_pat = rbind(total_pat , tmp_summary)
  
}

for (num_CancerType in Cancerlist2) {
  
  main.path_tc = paste0(filepath, "00.data/total/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  
  data_tc_link = readRDS(file = paste0(main.path_tc,"/",CancerType,"_pathwaylink_all_0or1.rds"))
  
  if ("Not Reported" %in% data_tc_link$vitalstatus ) {
    if (length(which(data_tc_link$vitalstatus == "Not Reported")) == 1) {
      data_tc_link = data_tc_link[-which(data_tc_link$vitalstatus == "Not Reported"),]
    } else {
      while ("Not Reported" %in% data_tc_link$vitalstatus) {
        data_tc_link = data_tc_link[-which(data_tc_link$vitalstatus == "Not Reported")[1],]
      }
      
    }
    
  }
  tmp_summary = data.frame(cancer_type = CancerType , total_length = nrow(data_tc_link) , bias = sum(data_tc_link$vitalstatus == "Alive" ) / length(data_tc_link$vitalstatus))
  total_pat = rbind(total_pat , tmp_summary)
  
}


total_pat = total_pat %>% filter(row_number() < n())
total_pat$cancer_type

total_pat = total_pat %>%
  mutate(cancer_type = str_remove(cancer_type, 'TCGA-'))

library(ggplot2)
library(ggrepel)

filtered_data <- total_pat %>%
  filter(total_length >= 200, bias > 0.125, bias < 0.875)

ggplot(data = total_pat, aes(x = total_length , y = bias, label = cancer_type))+
  geom_point(aes(color = ifelse(total_length >= 200 & bias > 0.125 & bias < 0.875, "Target", "Non-Target"))) +
  geom_text_repel(aes(label = cancer_type, color = ifelse(total_length >= 200 & bias > 0.125 & bias < 0.875, "Target", "Non-Target"),
                      fontface = ifelse(total_length >= 200 & bias > 0.125 & bias < 0.875, "bold", "plain")))+
  theme_bw() + 
  geom_vline(xintercept =200, color = "red", linetype = 2) + 
  geom_hline(yintercept = c(0.875, 0.125), color = "red", linetype = 2) +
  scale_color_manual(values = c("Target" = "#3C5488FF", "Non-Target" = "grey")) +
  theme(legend.position = "none")


