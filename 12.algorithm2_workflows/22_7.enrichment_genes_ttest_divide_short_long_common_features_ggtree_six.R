library(stringr)
library(data.table)
library(pheatmap)
library(readxl)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(dplyr)
library(readxl)
library(tidyverse)
library(rrvgo)
library(treemapify)

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
sce_path = "/mnt/gluster_server/data/raw/TCGA_data/00.data/"

surv_total_results = read_xlsx("~/nas/04.Results/Total_results_survpval2.xlsx")

# for fic
fig_path = paste0(filepath,"04.Results/GOenrichment_test/treemap_edit/five_category/")
if(!dir.exists(fig_path)){
  dir.create(fig_path)
  print(paste0("Created folder: ", fig_path))
} else {
  print(paste0("Folder already exists: ", fig_path))
}
setwd(fig_path)

# cf.) cell-cell adhesion, protein complex assembly
metabolism = c("metabolic ", " ion ","cadmium ion", "glycolytic ", "tricarboxylic acid",
               "glucose ", "catabolic ", "ATP ", 
               "proton transmembrane ", "acetyl-CoA", "electron transport","biosynthetic process","hexose","pentose-phosphate",
               "NADH regeneration","urea cycle","hydroxylation","NADH regeneration")

immune = c("cytokine", "leukocyte", "immune", "T cell", "lymphocyte", "Fc receptor", "cell activation","mononuclear cell",
           " lipopolysaccharide", " defense response", " antigen", 
           "phagocytosis", "response to virus","hemopoiesis", "chemotaxis", "cell killing", "autophagy", "MHC","repair","immunity", "response to biotic stimulus",
           "interferon","viral","inflammatory","bacterial origin","repair","complement activation")

signaling = c("signaling", "kinase", "cascade", "autophosphorylation","phosphorylation", "dephosphorylation","cellular response", "protein binding",
              "secretion","signal","signal transduction", "transcription factor","GTPase","peptide transport",
              "SMAD protein", "protein transport", "intracellular transport","zymogen","chemical synaptic transmission","immunoglobulin")

cell_cycle_develop = c("cell cycle","nuclear division", "cell division", "proliferation", "cell growth", "gland development", "gilogenesis", "cell development"
                       ,"renal system development", "protein localization", "angiogenesis","DNA replication", "differentiation", "epithelial tube",
                       "cell projection", "cellular component","pericardium dyevelopment",
                       "anterior/posterior pattern specification","ureteric bud development"," development","developmental growth","cell fate",
                       "telomere")

migration = c("migration", "actin", "lamellipodium", "epithelial to mesenchymal", "fiber assembly","mesenchyme development")  

apoptosis = c("cell death", "apoptosis","neuron death")


color_category = c("#7E6148FF","#8491B4FF","#00A087FF","#4DBBD5FF","#E64B35FF","#3C5488FF")
num_CancerType = "18.TCGA-LUAD"
for (num_CancerType in Cancerlist) {
  
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  
  # call input
  
  short_long_features = read_csv(paste0(filepath, "04.Results/short_long/ttest_common/",CancerType, "_critical_features_short_long_common.csv"))

  for (class_feature in unique(short_long_features$classification)) {
     assign(paste0("reducedTerms_",class_feature), read_rds(paste0(CancerType,"_reducedTerms",class_feature,".rds")))

    if (length(unique(get(paste0("reducedTerms_",class_feature))$parentTerm)) > 25) {
      assign(paste0("reducedTerms_",class_feature,"_edit"), 
             get(paste0("reducedTerms_",class_feature))[which(get(paste0("reducedTerms_",class_feature))$parentTerm %in%
                                                                unique(get(paste0("reducedTerms_",class_feature))[order(get(paste0("reducedTerms_",class_feature))$score, decreasing = T),]$parentTerm)[1:25]),] )
    } else {
      assign(paste0("reducedTerms_",class_feature,"_edit"), get(paste0("reducedTerms_",class_feature)))
    }

    assign(paste0("reducedTerms_",class_feature,"_edit"),get(paste0("reducedTerms_",class_feature,"_edit")) %>% 
             mutate(total_imp = case_when(parentTerm %in%
                                            get(paste0("reducedTerms_",class_feature,"_edit"))$parentTerm[grep(paste(metabolism, collapse = "|"), get(paste0("reducedTerms_",class_feature,"_edit"))$parentTerm)] ~ "meta",
                                          parentTerm %in%
                                            get(paste0("reducedTerms_",class_feature,"_edit"))$parentTerm[grep(paste(immune, collapse = "|"), get(paste0("reducedTerms_",class_feature,"_edit"))$parentTerm)] ~ "immune",
                                          parentTerm %in%
                                            get(paste0("reducedTerms_",class_feature,"_edit"))$parentTerm[grep(paste(signaling, collapse = "|"), get(paste0("reducedTerms_",class_feature,"_edit"))$parentTerm)] ~ "signaling",
                                          parentTerm %in%
                                            get(paste0("reducedTerms_",class_feature,"_edit"))$parentTerm[grep(paste(cell_cycle_develop, collapse = "|"), get(paste0("reducedTerms_",class_feature,"_edit"))$parentTerm)] ~ "cell_cycle_develop",
                                          parentTerm %in%
                                            get(paste0("reducedTerms_",class_feature,"_edit"))$parentTerm[grep(paste(migration, collapse = "|"), get(paste0("reducedTerms_",class_feature,"_edit"))$parentTerm)] ~ "migration",
                                          parentTerm %in%
                                            get(paste0("reducedTerms_",class_feature,"_edit"))$parentTerm[grep(paste(apoptosis, collapse = "|"), get(paste0("reducedTerms_",class_feature,"_edit"))$parentTerm)] ~ "apoptosis",
                                          .default = "not_include"))) 
   
    # assign(paste0("total_",class_feature,"_imp"),unique(get(paste0("reducedTerms_",class_feature,"_edit"))$total_imp))

    assign(paste0("reducedTerms_",class_feature,"_add"),get(paste0("reducedTerms_",class_feature,"_edit")) %>%
             mutate(color = case_when( total_imp == "meta" ~ "#7E6148FF",
                                       total_imp == "immune" ~ "#8491B4FF",
                                       total_imp == "signaling" ~ "#00A087FF",
                                       total_imp == "cell_cycle_develop" ~ "#4DBBD5FF",
                                       total_imp == "migration" ~ "#E64B35FF",
                                       total_imp == "apoptosis" ~ "#3C5488FF",
                                       
                                       .default = "#999999")) %>%
             mutate(subgroup_font_size = case_when( total_imp == "meta" ~ 25,
                                                    total_imp == "immune" ~ 25,
                                                    total_imp == "signaling" ~ 25,
                                                    total_imp == "cell_cycle_develop" ~ 25,
                                                    total_imp == "migration" ~ 25,
                                                    total_imp == "apoptosis" ~ 25,
                                                    
                                                    .default = 10))%>%
             
             mutate(subgroup_font_fontface = case_when(total_imp == "meta" ~ "bold",
                                                       total_imp == "immune" ~ "bold",
                                                       total_imp == "signaling" ~ "bold",
                                                       total_imp == "cell_cycle_develop" ~ "bold",
                                                       total_imp == "migration" ~ "bold",
                                                       total_imp == "apoptosis" ~ "bold", 
                                                       
                                                       .default = "plain")))
    
                             
                                           
    png(filename = paste0(CancerType,"_tree_",class_feature,"_cut_50_0.7_ttest_six_category.png"),
        width = 35, height = 35,  units = "cm" ,pointsize = 12,
        bg = "white", res = 1200, family = "")
    
    tree_mat =  ggplot(get(paste0("reducedTerms_",class_feature,"_add")), aes(area = score, fill = color,
                                                                              label = term, # country
                                                                              subgroup = parentTerm)) + # continent
      scale_fill_identity() +
      # 1. Draw country borders and fill colors
      geom_treemap(colour = "#FFFFFFDD", start = "topleft") +
      # 2. Draw continent borders
      geom_treemap_subgroup_border(
        # colour = "#00000080",
        start = "topleft",
        colour = "gray100",
        size = 5) +
      # 3. Print continent text
      
      geom_treemap_subgroup_text(place = "centre", colour = "black",
                                 reflow = T,
                                 fontface = get(paste0("reducedTerms_",class_feature,"_add"))$subgroup_font_fontface,
                                 # min.size = 0
                                 start = "topleft",
                                 grow = F,
                                 size= get(paste0("reducedTerms_",class_feature,"_add"))$subgroup_font_size) +
      # 4. Print country text
      geom_treemap_text(colour = "#FFFFFFDD",
                        place = "centre",
                        alpha = 0.5,
                        grow = F,
                        reflow = T,
                        start = "topleft") +
      # option
      theme(legend.position = 0)
    
    print(tree_mat)
    dev.off()
    
  }
}