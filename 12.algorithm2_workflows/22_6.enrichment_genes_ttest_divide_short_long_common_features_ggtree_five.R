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

wordcloudPlot_edit = function (reducedTerms, onlyParents = TRUE, ...) 
{
  if (!all(sapply(c("wordcloud", "tm", "slam"), requireNamespace, 
                  quietly = TRUE))) {
    stop("Package wordcloud and/or its dependencies (tm, slam) not available. ", 
         "Consider installing it before using this function.", 
         call. = FALSE)
  }
  if (onlyParents) {
    x <- tm::Corpus(tm::VectorSource(reducedTerms$term[reducedTerms$parent == 
                                                         reducedTerms$go]))
  }
  else {
    x <- tm::Corpus(tm::VectorSource(reducedTerms$term))
  }
  tdm <- tm::TermDocumentMatrix(x, control = list(removePunctuation = TRUE, 
                                                  stopwords = TRUE))
  m <- as.matrix(tdm)
  m <- m[which(!rownames(m) %in% c("regulation", "positive", "response",
                                   "development","process", "cell", "pathway",
                                   "activity","negative","involved","system","protein","via","levels","acid"
                                   ,"cellular","gene","eye","animal","precursor","class")),]
  v <- sort(rowSums(m), decreasing = TRUE)
  d <- data.frame(word = names(v), freq = v)
  wordcloud::wordcloud(d$word, d$freq, ...)
}

metabolism = c("metabolic ", " ion ","cadmium ion", "glycolytic ", "tricarboxylic ",
               "glucose ", "catabolic ", "ATP ", 
               "proton transmembrane ", "acetyl-CoA biosysthetic")

immune = c("cytokine", "leukocyte", "immune", "T cell", "lymphocyte", "Fc receptor", "cell activation","mononuclear cell",
           " lipopolysaccharide", " defense response", " antigen", "positive regulation of cell adhesion", " epithelial cell proliferation",
           "phagocytosis", "response to virus","hemopoiesis", "chemotaxis", "cell killing", "autophagy")

signaling = c("signaling", "kinase", "cascade", "autophosphorylation","phosphorylation", "dephosphorylation","cellular response", "protein binding",
              "secretion","signal","signal transduction")

cell_cycle_develop = c("cell cycle","nuclear division", "cell division", "proliferation", "cell growth", "gland development", "gilogenesis", "cell development"
                       ,"renal system development", "protein localization", "angiogenesis","DNA replication")
migration = c("migration", "actin", "lamellipodium")  

apoptosis = c("cell death", "apoptosis","neuron death")

# num_CancerType = "04.TCGA-CESC"
# for all
for (num_CancerType in Cancerlist) {
  
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  
  # call input
  
  short_long_features = read_csv(paste0(filepath, "04.Results/short_long/ttest_common/",CancerType, "_critical_features_short_long_common.csv"))
  
  for (class_feature in unique(short_long_features$classification)) {
    assign(paste0(class_feature,"_features"),short_long_features[which(short_long_features$classification == class_feature),]$variable )
    
    assign(paste0(class_feature,"_genes"), c())
    
    for (cf in get(paste0(class_feature,"_features"))) {
      count <- str_count(cf, "P")
      if (count == 1) {
        assign(paste0("tmp_",class_feature,"_genes"),single_genes[which(single_genes$Pathway == cf),]$Genes)
        
      } else {
        assign(paste0("tmp_",class_feature,"_genes"),link_genes_filtered_df[which(link_genes_filtered_df$Pathway == cf),]$Genes)
      }
      assign(paste0(class_feature,"_genes"),c(get(paste0(class_feature,"_genes")), get(paste0("tmp_",class_feature,"_genes")))) 
    }
    
    assign(paste0(class_feature,"_genes_en"),data.frame(AnnotationDbi::select(org.Hs.eg.db, get(paste0(class_feature,"_genes")), 'ENTREZID', 'SYMBOL')[
      which(!is.na(AnnotationDbi::select(org.Hs.eg.db, get(paste0(class_feature,"_genes")), 'ENTREZID', 'SYMBOL')$ENTREZID)),],row.names = NULL))
    
    assign(paste0(class_feature,"_enrichGO"),enrichGO(get(paste0(class_feature,"_genes_en"))$ENTREZID ,
                                                      OrgDb = org.Hs.eg.db ,
                                                      keyType = "ENTREZID" , 
                                                      ont = "BP" ,
                                                      pvalueCutoff  = 0.05,
                                                      qvalueCutoff = 0.01))
    
    saveRDS(get(paste0(class_feature,"_enrichGO")),paste0(CancerType,"_",class_feature,"_enrichGO.rds"))
    
    assign(paste0(class_feature,"_enrichGO_df"), as.data.frame(get(paste0(class_feature,"_enrichGO"))))
    
    assign(paste0("scores_",class_feature), setNames(-log10(get(paste0(class_feature,"_enrichGO_df"))$qvalue), get(paste0(class_feature,"_enrichGO_df"))$ID))
    
    assign(paste0("simMatrix_",class_feature), calculateSimMatrix(get(paste0(class_feature,"_enrichGO_df"))$ID,
                                                                  orgdb=org.Hs.eg.db,
                                                                  ont="BP",
                                                                  method="Rel"))
    
    saveRDS(get(paste0("simMatrix_",class_feature)),paste0(CancerType,"_","simMatrix_",class_feature,".rds"))
    
    if (sum(is.na(get(paste0("scores_",class_feature)))) != 0) {
      get(paste0("scores_",class_feature))[which(is.na(get(paste0("scores_",class_feature))))] = 0
    }
    
    if (sum(get(paste0("scores_",class_feature))) == 0) {
      assign(paste0("reducedTerms_",class_feature), reduceSimMatrix(get(paste0("simMatrix_",class_feature)),
                                                                    threshold=0.7,
                                                                    orgdb="org.Hs.eg.db"))
      
    } else {
      assign(paste0("reducedTerms_",class_feature), reduceSimMatrix(get(paste0("simMatrix_",class_feature)),
                                                                    threshold=0.7,
                                                                    get(paste0("scores_",class_feature)),
                                                                    orgdb="org.Hs.eg.db"))
      
    }
    
    saveRDS(get(paste0("reducedTerms_",class_feature)),paste0(CancerType,"_","reducedTerms",class_feature,".rds"))
    
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
                                          .default = "not_include"))) 
    
    assign(paste0("total_",class_feature,"_imp"),unique(get(paste0("reducedTerms_",class_feature,"_edit"))$total_short_imp))
    
    color_category = c("#7E6148FF","#8491B4FF","#00A087FF","#4DBBD5FF","#E64B35FF")
    n=0
    # category = "meta"
    for (category in c("meta", "immune","signaling","cell_cycle_develop","migration")) {
      n = n+1
      
      assign(paste0("reducedTerms_",class_feature,"_add"),get(paste0("reducedTerms_",class_feature,"_edit")) %>%
               mutate(color = case_when( total_imp == category ~ color_category[n],
                                         .default = "#999999")) %>%
               mutate(subgroup_font_size = case_when( total_imp == category ~ 25,
                                                      .default = 10))%>%
               mutate(subgroup_font_fontface = case_when( total_imp == category ~ "bold",
                                                          .default = "plain")))
      
      # fiq
      
      png(filename = paste0(CancerType,"_tree_",class_feature,"_cut_50_0.7_ttest_",category,".png"),
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
}
  
  
  
  
  
  ###################################################
  
  
  