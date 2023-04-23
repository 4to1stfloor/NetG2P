#loading R packages
library(data.table)
library(igraph)
library(survival)
library(survminer)

# global setting
filepath = paste0("/home/seokwon/nas/")
ref_path = paste0(filepath, "/99.reference/")

num.edge.multi = 10 #or 100 (for now)
random.num = 150
sourcepath = paste0(ref_path, "RDPN/",c(num.edge.multi),"X" ) 

# 1,2 setting
#cancer.type = "12.PRAD" 
ppi.date = "20220117" 
path.dat = "KEGG" #"Iorio2018" #"DKShin"
path.file = paste0(ref_path ,"Kegg_pathway_genes.rds") 
#loading PPI Backbone
backbone.file = paste0(ref_path,"ppi_backbone_" ,ppi.date,".rds")
g.ppi.conn = readRDS(backbone.file)
#loading CHG data 
path = readRDS(path.file)

# 2 setting
level = "_net" # ""
filt.genes.n = 5 # at pathway-link level 

# 3 setting
## extracting pathway links with having more than 5 genes in a pathway link 
dat = readRDS(paste0(ref_path,"KEGG_pathway_links_intersect_genes.rds"))

# zscore normalization
tcga.calc.zscore = function(sce, target.genes){
  message("Calculating z-score with respect to all diploid cells. Version 2022.11.07")
  common.genes = intersect(rownames(sce), target.genes)
  if (length(common.genes) == 0) {
    stop("None of the target genes found in sce. Please check your nomenclature")
  } else if (all(common.genes %in% target.genes)) {
    print("Check!")
  }
  sce.sub = subset(sce, rownames(sce) %in% target.genes,)
  #i am not checking assay names.
  count.mat = assay(sce.sub, 1)
  cnv.mat = assay(sce.sub, 3)
  z.mat = matrix(data = NA, nrow = length(target.genes), ncol = ncol(count.mat))
  colnames(z.mat) = colnames(count.mat)
  rownames(z.mat) = target.genes
  for (i in 1:nrow(count.mat)) {
    idx.di = which(cnv.mat[i,] == 2)
    query.mean = mean(count.mat[i, idx.di], na.rm = T)
    query.sd = sd(count.mat[i, idx.di], na.rm = T)
    idx.symb = which(rownames(z.mat) == rownames(count.mat)[i])
    z.mat[i,] = (count.mat[i,] - query.mean)/query.sd
  }
  return(z.mat)
}

###############################################################################################################################################
## This script is for hypergeometric test (both pathway and pathway link levels)
###########################################################################################################################################################################
### Direct to phenotype (Net.prop & RDPN --> GC --> 54 KEGG pathway links)
###########################################################################################################################################################################

## Function define

###########################################################################################################################################################################
### 1.function of pathway-level (e.g., P1, P2, P3, ..., P54) 

GC.hyper.pathway.direct  <- function(CancerType, ppi.date, path.dat, path.file) {
  
  for (s in 1:ncol(pval.mat)) {
    pval.sample = names(pval.mat[,s]) 
    group2 = names(pval.mat[,s][which(abs(pval.mat[,s]) >= 1)]) 
    g.ppi.spec = induced.subgraph(g.ppi.conn, group2)
    c.ppi.spec = components(g.ppi.spec) 
    idx.max = which(c.ppi.spec$csize == max(c.ppi.spec$csize))
    #g.ppi.gc = induced.subgraph(g.ppi.spec, names(which(c.ppi.spec$membership == idx.max)))
    genes.gc = names(which(c.ppi.spec$membership == idx.max))
    group2 = genes.gc
    for ( p in 1:length(unique(path$Pathway))) {
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
  
  saveRDS(enr.mat, paste0(main.path_tc, "/phyper_enrichment_exp_", length(unique(path$Pathway)), "_", path.dat, "_GC_",ppi.date, ".rds"))
  
  # transforming the p-values obtained from the hypergeometric tests with log10(pvalue)
  #enr = readRDS(paste0("/home/jisushin/project/", cancer.type, "/phyper_enrichment_", chg.dat, "_10chg_GC_", ppi.date, "_grn.rds"))
  #enr.trans = -log10(enr) 
  #saveRDS(enr.trans, paste0("/home/jisushin/project/", cancer.type, "/phyper_enrichment_",chg.dat, "_10chg_GC_", ppi.date, "_log_trans_grn.rds"))
}
###########################################################################################################################################################################

###########################################################################################################################################################################
### 2.function of pathway-link level (e.g., P1-P2, P3-P10, ..., )

GC.hyper.path.net.direct  <- function(CancerType, ppi.date, path.dat, path.file) {
  
  # giant clustering 
  for (s in 1:ncol(pval.mat)) {
    
    pval.sample = names(pval.mat[,s]) 
    group2 = names(pval.mat[,s][which(abs(pval.mat[,s]) >= 1)]) 
    g.ppi.spec = induced.subgraph(g.ppi.conn, group2)
    c.ppi.spec = components(g.ppi.spec) 
    idx.max = which(c.ppi.spec$csize == max(c.ppi.spec$csize))
    #g.ppi.gc = induced.subgraph(g.ppi.spec, names(which(c.ppi.spec$membership == idx.max)))
    genes.gc = names(which(c.ppi.spec$membership == idx.max))
    group2 = genes.gc
    for ( p in 1:length(unique(path$Pathway))) {
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
  
  saveRDS(dat, paste0(main.path_tc, "/phyper_enrichment_exp_", length(unique(path$Pathway)), "_", path.dat, "_net_GC_", ppi.date, ".rds"))
}
###########################################################################################################################################################################

###########################################################################################################################################################################
### 3. function of obtaining p values from the survival analysis (pathway level)

cal.surv.each.p <- function(CancerType, ppi.date, path.dat, filt.genes.n) { 
  
  hyper.mat = readRDS(paste0(main.path_tc, "/phyper_enrichment_exp_", length(unique(path$Pathway)), "_", path.dat, "_GC_",ppi.date, ".rds"))
  #print(hyper.mat)
  # if (level != ""){
  #   sig.pair = dat$name[which(dat$n_genes > as.numeric(filt.genes.n))] 
  #   idx = match(sig.pair, rownames(hyper.mat))
  #   hyper.filt = hyper.mat[idx, ]
  #   } else {
  hyper.filt = hyper.mat
  # }
  
  #exp.log.mat = readRDS(paste0(CancerType, '_exp_TPM_mat_filt_log_gene_symbol.rds')) #56404 973 
  mut.mat = readRDS(paste0(main.path_tc,"/", CancerType ,'_mut_count_filt_data.rds'))#16131 985
  mut.count = as.data.frame(colSums(mut.mat))
  colnames(mut.count) = "counts" 
  
  if (CancerType %in% c("TCGA-COADREAD","TCGA-KIDNEY")) {
    mut.count$patient_id = substr(rownames(mut.count), 1, 12) 
  } else {
    mut.count$patient_id = substr(rownames(mut.count), 1, 16) 
  }
  
  mut.filt = mut.count[which(mut.count$counts < 1000), ]
  
  
  if (CancerType %in% c("TCGA-COADREAD","TCGA-KIDNEY")) {
    int.id = intersect(substr(colnames(hyper.filt),1,12), mut.filt$patient_id) 
    idx.id = match(int.id, substr(colnames(hyper.filt),1,12))
  } else {
    int.id = intersect(colnames(hyper.filt), mut.filt$patient_id) 
    idx.id = match(int.id, colnames(hyper.filt))
  }
  
  hyper.filt = hyper.filt[, idx.id ]
  
  #loading clinical information for survival analysis 
  # Patient ID -> submitter_id, Overall_survival -> vital_status , Overall_survival(Month) -> days_to_death or days_to_last_follow_up
  # time -> overall_survival  
  
  cli = fread(paste0(ref_path,"all_clin_indexed.csv")) 
  
  if (CancerType == "TCGA-COADREAD") {
    cli_surv = cli[cli$project %in% c("TCGA-COAD","TCGA-READ"),
                   c("submitter_id",
                     "vital_status",
                     "days_to_death",
                     "days_to_last_follow_up")]
  } else if (CancerType == "TCGA-KIDNEY") {
    cli_surv = cli[cli$project %in% c("TCGA-KIRP","TCGA-KIRC","TCGA-KICH"),
                   c("submitter_id",
                     "vital_status",
                     "days_to_death",
                     "days_to_last_follow_up")]
  } else {
    cli_surv = cli[cli$project == CancerType,
                   c("submitter_id",
                     "vital_status",
                     "days_to_death",
                     "days_to_last_follow_up")]
  }
  
  cli_surv$deceased = cli_surv$vital_status == "Dead"
  cli_surv$overall_survival = ifelse(cli_surv$deceased,
                                     cli_surv$days_to_death,
                                     cli_surv$days_to_last_follow_up)
  
  cli_surv = cli_surv[!is.na(cli_surv$overall_survival),]
  cli_surv = cli_surv[!is.na(cli_surv$vital_status),]
  
  cli_surv$status = NA
  cli_surv$status[which(cli_surv$vital_status == "Dead")] = 1
  cli_surv$status[which(cli_surv$vital_status == "Alive")] = 0
  
  cli_surv$status = as.numeric(cli_surv$status)

  wrong_cluster = vector()
  pval.com =c()
  for (p in 1:as.numeric(dim(hyper.filt)[1])) {
    hyper.each = hyper.filt[p,]
    hyper.dat = as.data.frame(cbind(id = substr(names(hyper.each), 1, 12), pvals= as.numeric(hyper.each)))
    #hyper.dat$group = "" 
    hyper.dat['group'] = ""
    hyper.dat$group[which(hyper.dat$pvals <= 0.05)] = "enriched"
    hyper.dat$group[which(hyper.dat$pvals > 0.05)] = "not_enriched"
    int = intersect(hyper.dat$id, cli_surv$submitter_id)
    idx = match(int, cli_surv$submitter_id) 
    idx2 = match(int, hyper.dat$id) 
    # print("hear")
    cli_surv_filt = cli_surv[idx,]
    
    hyper.dat.filt = hyper.dat[idx2,]
    cli_surv_filt$cluster = hyper.dat.filt$group
    
    if (length(unique(cli_surv_filt$cluster)) == 1) {
      wrong_cluster = c(wrong_cluster, p) 
      next
    }
    
    fit = survminer::surv_fit(Surv(overall_survival,status) ~ cluster, data = cli_surv_filt) 

    pval.each = surv_pvalue(fit)[2]
    # print(p)
    if( length(pval.com) == 0) { 
      pval.com = as.numeric(pval.each)
    } else {
      pval.com = c(pval.com, as.numeric(pval.each))
      #print(pval.com)
    } 
  }
  
  pval.com = as.data.frame(pval.com) 
  if (length(wrong_cluster) == 0) {
    rownames(pval.com) = rownames(hyper.filt)
  } else if (length(wrong_cluster) != 0 ) {
    rownames(pval.com) = rownames(hyper.filt)[-wrong_cluster] 
  } else {
    print("I dont know")
  }
  
  colnames(pval.com) = gsub(".*\\.", "", CancerType)
  # print(pval.com)
  
  saveRDS(pval.com, paste0(main.path_tc, "/surv_pvals_54_exp_", path.dat, "_GC_", ppi.date, "_filt", filt.genes.n, ".rds"))
  
  wrong_cluster = vector()
}
###########################################################################################################################################################################

###########################################################################################################################################################################
### 4. function of obtaining p values from the survival analysis (pathway-link level)

cal.surv.p <- function(CancerType, ppi.date, path.dat, level, filt.genes.n) { 
  
  hyper.mat = readRDS(paste0(main.path_tc, "/phyper_enrichment_exp_", length(unique(path$Pathway)), "_", path.dat, "_net_GC_", ppi.date, ".rds"))
  #print(hyper.mat)
  if (level != ""){
    sig.pair = dat$name[which(dat$n_genes > as.numeric(filt.genes.n))] 
    idx = match(sig.pair, rownames(hyper.mat))
    hyper.filt = hyper.mat[idx, ]
  } else {
    hyper.filt = hyper.mat
  }
  
  

  
  #exp.log.mat = readRDS(paste0(CancerType, '_exp_TPM_mat_filt_log_gene_symbol.rds')) #56404 973 
  mut.mat = readRDS(paste0(main.path_tc,"/", CancerType ,'_mut_count_filt_data.rds'))#16131 985
  mut.count = as.data.frame(colSums(mut.mat))
  colnames(mut.count) = "counts" 
  
  if (CancerType %in% c("TCGA-COADREAD","TCGA-KIDNEY")) {
    mut.count$patient_id = substr(rownames(mut.count), 1, 12) 
  } else {
    mut.count$patient_id = substr(rownames(mut.count), 1, 16) 
  }
  
  mut.filt = mut.count[which(mut.count$counts < 1000), ]
  
  if (CancerType %in% c("TCGA-COADREAD","TCGA-KIDNEY")) {
    int.id = intersect(substr(colnames(hyper.filt),1,12), mut.filt$patient_id) 
    idx.id = match(int.id, substr(colnames(hyper.filt),1,12))
  } else {
    int.id = intersect(colnames(hyper.filt), mut.filt$patient_id) 
    idx.id = match(int.id, colnames(hyper.filt))
  }
  
  hyper.filt = hyper.filt[, idx.id ]
  
  #loading clinical information for survival analysis 
  # Patient ID -> submitter_id, Overall_survival -> vital_status , Overall_survival(Month) -> days_to_death or days_to_last_follow_up
  # time -> overall_survival  
  
  cli = fread(paste0(ref_path,"all_clin_indexed.csv")) 
  
  if (CancerType == "TCGA-COADREAD") {
    cli_surv = cli[cli$project %in% c("TCGA-COAD","TCGA-READ"),
                   c("submitter_id",
                     "vital_status",
                     "days_to_death",
                     "days_to_last_follow_up")]
  } else if (CancerType == "TCGA-KIDNEY") {
    cli_surv = cli[cli$project %in% c("TCGA-KIRP","TCGA-KIRC","TCGA-KICH"),
                   c("submitter_id",
                     "vital_status",
                     "days_to_death",
                     "days_to_last_follow_up")]
  } else {
    cli_surv = cli[cli$project == CancerType,
                   c("submitter_id",
                     "vital_status",
                     "days_to_death",
                     "days_to_last_follow_up")]
  }
  
  cli_surv$deceased = cli_surv$vital_status == "Dead"
  cli_surv$overall_survival = ifelse(cli_surv$deceased,
                                     cli_surv$days_to_death,
                                     cli_surv$days_to_last_follow_up)
  
  cli_surv = cli_surv[!is.na(cli_surv$overall_survival),]
  cli_surv = cli_surv[!is.na(cli_surv$vital_status),]
  
  cli_surv$status = NA
  cli_surv$status[which(cli_surv$vital_status == "Dead")] = 1
  cli_surv$status[which(cli_surv$vital_status == "Alive")] = 0
  
  cli_surv$status = as.numeric(cli_surv$status)
  
  wrong_cluster = vector()
  pval.com = c()
  for (p in 1:as.numeric(dim(hyper.filt)[1])) {
    hyper.each = hyper.filt[p,]
    hyper.dat = as.data.frame(cbind(id = substr(names(hyper.each), 1, 12), pvals= as.numeric(hyper.each)))
    #hyper.dat$group = "" 
    hyper.dat['group'] = ""
    hyper.dat$group[which(hyper.dat$pvals <= 0.05)] = "enriched"
    hyper.dat$group[which(hyper.dat$pvals > 0.05)] = "not_enriched"
    int = intersect(hyper.dat$id, cli_surv$submitter_id)
    idx = match(int, cli_surv$submitter_id) 
    idx2 = match(int, hyper.dat$id) 
    # print("hear")
    cli_surv_filt = cli_surv[idx,]
    
    # print("hear")
    hyper.dat.filt = hyper.dat[idx2,]
    cli_surv_filt$cluster = hyper.dat.filt$group

    if (length(unique(cli_surv_filt$cluster)) == 1) {
      wrong_cluster = c(wrong_cluster, p) 
      next
    }
    
    fit = survminer::surv_fit(Surv(overall_survival,status) ~ cluster, data = cli_surv_filt) 
    #fit = survfit(Surv(overall_survival,status) ~ cluster, data = cli_surv_filt)
    # print(p)
    pval.each = surv_pvalue(fit)[2]
    # print(p)
    if(  length(pval.com) == 0) { 
      pval.com = as.numeric(pval.each)
    } else {
      pval.com = c(pval.com, as.numeric(pval.each))
      #print(pval.com)
    } 
  }
  
  pval.com = as.data.frame(pval.com) 
  if (length(wrong_cluster) == 0) {
    rownames(pval.com) = rownames(hyper.filt)
  } else if (length(wrong_cluster) != 0 ) {
    rownames(pval.com) = rownames(hyper.filt)[-wrong_cluster] 
  } else {
    print("I dont know")
  }
  
  colnames(pval.com) = gsub(".*\\.", "", CancerType)
  print(pval.com)
  if (level != "") {
    saveRDS(pval.com, paste0(main.path_tc, "/surv_pvals_54_exp_", path.dat, level, "_GC_", ppi.date, "_filt", filt.genes.n, ".rds"))
  } else {
    saveRDS(pval.com, paste0(main.path_tc,"/surv_pvals_54_exp_", path.dat, level, "_GC_", ppi.date, ".rds"))
  }
  
  wrong_cluster = vector()
}
###########################################################################################################################################################################

Cancerlist = dir(paste0(filepath, "/00.data/filtered_TCGA/"))
sce_path = "/mnt/gluster_server/data/raw/TCGA_data/00.data/"
# Cancerlist = Cancerlist[-7]

# pval.mat 이거를 expression matrix로 변경해야함 hvg는 위에 뽑음 

for (num_CancerType in Cancerlist) {
  
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  
  # call input
  
  sce = readRDS(paste0(sce_path,num_CancerType,"/", CancerType,".sce_raw_snv_cnv.rds"))
  sce_for_hvgs = sce 
  
  sce_for_hvgs = logNormCounts(sce_for_hvgs)
  sce_for_hvgs = subset(sce_for_hvgs, rownames(sce_for_hvgs) %in% V(g.ppi.conn)$name,)
  
  assays_to_remove <- c("counts", "snv", "cnv")
  for (assay_name in assays_to_remove) {
    assay(sce_for_hvgs, assay_name) <- NULL
  }
  
  dec.sce <- modelGeneVar(sce_for_hvgs)
  hvgs = getTopHVGs(dec.sce, n=3000)
  sce_exp = tcga.calc.zscore(sce = sce, hvgs)
  #loading RDPN output and function 2 output 
  sce_exp = as.data.frame(sce_exp)
  sce_exp = sce_exp[complete.cases(sce_exp),]
  sce_exp_filt = sce_exp
  
  pval.mat = as.matrix(sce_exp_filt) # expression of hvg!

  if (file.exists(paste0(main.path_tc, "surv_pvals_54_", path.dat, level, "_GC_", ppi.date, "_filt", filt.genes.n, ".rds"))) {
    next
  }
  # Execute function 1
  GC.hyper.pathway.direct(cancer.type, ppi.date, path.dat, path.file)
  
  # Execute function 2
  GC.hyper.path.net.direct(CancerType, ppi.date, path.dat, path.file) ## then "phyper_enrichment_54_KEGG_net_GC_20220117.rds" will be saved
  
  # Execute function 3
  cal.surv.each.p(CancerType, ppi.date, path.dat, filt.genes.n ) # this will generate "surv_pvals_54_KEGG_GC_20220117_filt5.rds" 
  
  # Execute function 4
  cal.surv.p(CancerType, ppi.date, path.dat, level, filt.genes.n ) # this will generate "surv_pvals_54_KEGG_net_GC_20220117_filt5.rds" 
  
}



