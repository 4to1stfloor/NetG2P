GC_projection = function(result_rdpn,
                         output_path = "final_res/",
                         random.num = 150){
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
  path = readRDS("reference/KEGG_pathway_genes.rds")
  g.ppi.conn = readRDS("reference/ppi_backbone_20220117.rds")
  
  ###############################################################################################################################################
  ## This script is for hypergeometric test (both pathway and pathway link levels)
  ###########################################################################################################################################################################
  ### Direct to phenotype (Net.prop & RDPN --> GC --> 54 KEGG pathway links)
  ###########################################################################################################################################################################
  
  ## Function define
  
  ###########################################################################################################################################################################
  ### 1.function of pathway-level (e.g., P1, P2, P3, ..., P54) 
  
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
      for ( p in 1:length(unique(path$pathway))) {
        pathway.name = paste0("P", p)
        group1 = intersect(path$Symbol[which(path$pathway == pathway.name)], pval.sample)
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
    
    saveRDS(enr.mat, paste0(output_path, "phyper_enrichment_", length(unique(path$pathway)), "_", path.dat, "_GC.rds"))
    
  }
  ###########################################################################################################################################################################
  
  ###########################################################################################################################################################################
  ### 2.function of pathway-link level (e.g., P1-P2, P3-P10, ..., )
  
  GC.hyper.path.net.direct  <- function(path.dat) {
    
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
      for ( p in 1:length(unique(path$pathway))) {
        pathway.name = paste0("P", p)
        if (p == 1) {
          pathway.names = pathway.name 
        } else {
          pathway.names = c(pathway.names, pathway.name)
        }
      }
      com = combn(pathway.names, 2)
      for (i in 1:ncol(com)){
        pair = com[,i]
        first = unique(path$Genes[which(path$Pathway == pair[1])])
        second = unique(path$Genes[which(path$Pathway == pair[2])])
        pair.name = paste0(pair[1], "-", pair[2])
        group1 = intersect(intersect(first, second), pval.sample) 
        total.genes = unique(c(pval.sample, group1))
        int = intersect(group1, group2) 
        enr = phyper(length(int)-1, length(group2), length(total.genes)-length(group2), length(group1), lower.tail = FALSE) 
        if ( i == 1) {
          enr.comb = enr
          names = pair.name 
        } else {
          enr.comb = c(enr.comb, enr) 
          names = c(names, pair.name)
        }
      }
      if ( s == 1) {
        dat = as.data.frame(enr.comb) 
      } else {
        dat = cbind(dat, as.data.frame(enr.comb))
      }
    }
    rownames(dat) = names 
    colnames(dat) = colnames(pval.mat) 
    
    saveRDS(dat, paste0(output_path, "phyper_enrichment_", length(unique(path$pathway)), "_", path.dat, "_net_GC.rds"))
  }
  ###########################################################################################################################################################################

  pval.mat = result_rdpn # output of RDPN.R 
  
  # Execute function 1
  GC.hyper.pathway.direct(path.dat)
  
  # Execute function 2
  GC.hyper.path.net.direct(path.dat) ## then "phyper_enrichment_54_KEGG_net_GC_20220117.rds" will be saved
  
}





  

