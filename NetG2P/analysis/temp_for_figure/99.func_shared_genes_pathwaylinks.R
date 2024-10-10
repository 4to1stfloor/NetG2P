## 

filepath = "/home/seokwon/nas/"
ref_path = paste0(filepath, "99.reference/")

num_links_genes = readRDS(file = paste0(ref_path,"/KEGG_pathway_links_intersect_genes.rds"))
genes_in_pathway = readRDS(file = paste0(ref_path,"/Kegg_pathway_genes.rds"))

link_genes = num_links_genes
link_genes$name = gsub("[[:punct:]]", "", link_genes$name)
rownames(link_genes) = link_genes$name
link_genes$name = NULL
link_genes$shared_genes = NA

# get pathway links
pathway_combinations <- combn(unique(genes_in_pathway$pathway), 2, simplify = FALSE)
pathwaylinks <- lapply(pathway_combinations, function(x) {
  paste0(x, collapse = "")
})
pathwaylinks <- unlist(pathwaylinks)

for (num in 1:length(pathway_combinations)) {
  first_pathway = pathway_combinations[[num]][1]
  second_pathway = pathway_combinations[[num]][2]
  tmp_shared_genes = genes_in_pathway[which(genes_in_pathway$pathway==first_pathway),]$Symbol[which(genes_in_pathway[which(genes_in_pathway$pathway==first_pathway),]$Symbol %in% 
                                                                           genes_in_pathway[which(genes_in_pathway$pathway==second_pathway),]$Symbol)]
  tmp_shared_genes = paste(tmp_shared_genes, collapse = ",")
  link_genes[num,]$shared_genes = tmp_shared_genes
   
}

link_genes$pathway_name = rownames(link_genes)
write.csv(link_genes , paste0(ref_path, "/KEGG_pathway_shared_genes.csv"))
saveRDS(link_genes , paste0(ref_path, "/KEGG_pathway_shared_genes.rds"))

### pathway gene for other format wo 0 shared genes

onco_for_pathway_gene = data.frame(Genes = NA , Pathway = NA )
onco_for_pathway_gene$Genes = NA
links_genes_wo = link_genes[which(link_genes$n_genes != 0),]
n = 0
for ( total_pathway in 1:length(links_genes_wo$shared_genes) ) {
  for (sh_genes in strsplit(links_genes_wo$shared_genes[total_pathway], ",")[[1]]) {
    n = n +1
    onco_for_pathway_gene[n,]$Genes = sh_genes
    onco_for_pathway_gene[n,]$Pathway = links_genes_wo$pathway_name[total_pathway]
  }
 
}
onco_for_pathway_gene <- data.frame(onco_for_pathway_gene, row.names = NULL)

write.csv(onco_for_pathway_gene , paste0(ref_path, "/KEGG_pathway_shared_each_gene.csv"))
saveRDS(onco_for_pathway_gene , paste0(ref_path, "/KEGG_pathway_shared_each_gene.rds"))