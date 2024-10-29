GC_projection = function(result_rdpn,
                         output_path = "./final_res/",
                         random.num = 150){
  
  ###########################################################################################################################################################################
  ## This script is for hypergeometric test (both pathway and pathway link levels)
  ###########################################################################################################################################################################
  ### Direct to phenotype (Net.prop & RDPN --> GC --> 54 KEGG pathway links)
  ###########################################################################################################################################################################
  ## Function define
  ###########################################################################################################################################################################
  ### 1.function of pathway-level (e.g., P1, P2, P3, ..., P54) 
  
  packages <- c("dplyr") 
  
  install_if_missing <- function(pkg) {
    if (!require(pkg, character.only = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
      library(pkg, character.only = TRUE)
    }
  }
  
  sapply(packages, install_if_missing)
  
  library(dplyr)
  
  GC.hyper.pathway.direct  <- function(path.dat) {
    
    for (s in 1:ncol(pval.mat)) {
      pval.thres = 0.05 
      pval.sample = names(pval.mat[,s]) 
      group2 = names(pval.mat[,s][which(pval.mat[,s] <= pval.thres)]) 
      g.ppi.spec = induced.subgraph(g.ppi.conn, group2)
      c.ppi.spec = components(g.ppi.spec) 
      idx.max = which(c.ppi.spec$csize == max(c.ppi.spec$csize))
      #g.ppi.gc = induced.subgraph(g.ppi.spec, names(which(c.ppi.spec$membership == idx.max)))
      genes.gc = names(which(c.ppi.spec$membership == idx.max))
      group2 = genes.gc
      
      pathway_length = path %>%
        filter(nchar(Pathway) <= 3) %>%
        pull(Pathway) %>%
        unique() %>%
        length()
      
      for ( p in 1:pathway_length) {
        
        pathway.name = paste0("P", p)
        group1 = intersect(path$Genes[which(path$Pathway == pathway.name)], pval.sample)
        total.genes = unique(c(pval.sample, group1))
        overlapped = intersect(group1, group2) 
        enr.value = phyper(length(overlapped)-1, length(group2), (length(total.genes) - length(group2)), length(group1), lower.tail = FALSE) 
        if ( p == 1) {
          enr.accum = enr.value
          name.accum = pathway.name  
        } else {
          enr.accum = c(enr.accum, enr.value) 
          name.accum = c(name.accum, pathway.name)
        }
      }
      if (s == 1) {
        enr.mat = as.data.frame(enr.accum) 
      } else { 
        enr.mat = cbind(enr.mat, as.data.frame(enr.accum))
      }
    }
    
    rownames(enr.mat) = name.accum
    colnames(enr.mat) = colnames(pval.mat)
    
    return(enr.mat)
    
  }
  ###########################################################################################################################################################################
  
  ###########################################################################################################################################################################
  ### 2.function of pathway-link level (e.g., P1-P2, P3-P10, ..., )
  
  GC.hyper.path.net.direct  <- function(path.dat) {
    
    tmp_pathway = path %>%
      filter(nchar(Pathway) <= 3)  %>%
      pull(Pathway) %>%
      unique()
    
    pathway_links = apply(combn(tmp_pathway, 2), 2, paste, collapse = "")
    
    # giant clustering 
    for (s in 1:ncol(pval.mat)) {
      pval.thres = 0.05 
      pval.sample = names(pval.mat[,s]) 
      group2 = names(pval.mat[,s][which(pval.mat[,s] <= pval.thres)]) 
      g.ppi.spec = induced.subgraph(g.ppi.conn, group2)
      c.ppi.spec = components(g.ppi.spec) 
      idx.max = which(c.ppi.spec$csize == max(c.ppi.spec$csize))
      #g.ppi.gc = induced.subgraph(g.ppi.spec, names(which(c.ppi.spec$membership == idx.max)))
      genes.gc = names(which(c.ppi.spec$membership == idx.max))
      group2 = genes.gc
      
      enr.comb = c()
      names = c()
      for (pathway_link in pathway_links){
        
        if (path %>% filter(Pathway == pathway_link ) %>% nrow() > 0) {
          
          group1 = intersect(path %>% filter(Pathway == pathway_link) %>% pull(Genes), pval.sample) 
          total.genes = unique(c(pval.sample, group1))
          int = intersect(group1, group2) 
          enr = phyper(length(int)-1, length(group2), length(total.genes)-length(group2), length(group1), lower.tail = FALSE) 
          
        } else {
          enr = 1 # because there are not intersect genes between two other pathway
        }
        
        # if ( i == 1) {
        #   enr.comb = enr
        #   names = pathway_link[i]
        # } else {
        #   enr.comb = c(enr.comb, enr) 
        #   names = c(names, pathway_link[i])
        # }
        enr.comb = c(enr.comb, enr) 
        names = c(names, pathway_link)
      }
      
      if ( s == 1) {
        dat = as.data.frame(enr.comb) 
      } else {
        dat = cbind(dat, as.data.frame(enr.comb))
      }
    }
    
    rownames(dat) = names 
    colnames(dat) = colnames(pval.mat) 
    
    return(dat)
  }
  ###########################################################################################################################################################################
  
  #loading R packages
  require(data.table)
  require(igraph)
  require(data.table)
  
  if(!dir.exists(output_path)){
    dir.create(output_path)
    print(paste0("Created folder: ", output_path))
  } else {
    print(paste0("Folder already exists: ", output_path))
  }
  
  path.dat = "KEGG" #"DKShin"
  path = readRDS("../data/reference/KEGG_dual_total_genes.rds")
  g.ppi.conn = readRDS("../data/reference/ppi_backbone_20220117.rds")
  pval.mat = result_rdpn # output of RDPN.R 
  
  # Execute function 1
  path_res = GC.hyper.pathway.direct(path.dat)
  
  # Execute function 2
  link_res = GC.hyper.path.net.direct(path.dat) ## then "phyper_enrichment_54_KEGG_net_GC_20220117.rds" will be saved
  
  path_res_t = as.data.frame(t(path_res))
  link_res_t = as.data.frame(t(link_res))
  
  if (all.equal(rownames(path_res_t), rownames(link_res_t))) {
    total_res = cbind(path_res_t, link_res_t)
    total_res = -log(total_res)
    saveRDS(total_res , paste0(output_path,"/GC_projection_final.rds"))
    write.csv(total_res, paste0(output_path,"/GC_projection_final.csv"))
    
  } else {
    print("pathway and pathwaylink data are not matched rownames. Will be provided seperated file.")
    path_res_t = -log(path_res_t)
    link_res_t = -log(link_res_t)
    
    write.csv(path_res_t , paste0(output_path,"/GC_projection_final_pathway.csv"))
    write.csv(link_res_t, paste0(output_path,"/GC_projection_final_pathway_link.csv"))
  }
  
  return(total_res)
  
}







