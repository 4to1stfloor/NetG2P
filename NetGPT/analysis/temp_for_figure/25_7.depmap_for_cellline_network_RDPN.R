
## This script is for running RDPN (random degree preserving network)  with new PPI backbone 

#loading tools 
require(igraph) 
require(data.table)
require(foreach) 
require(doParallel) 
library(tictoc)
library(tidyverse)
library(dplyr)
library(igraph)
library(tictoc) 
library(pcg)
library(SingleCellExperiment)

# setting
dep.snv.caller = function(depmap.ver = "23Q4", sample.name){
  library(maftools)
  depmap.path = "/mnt/gluster_server/data/reference/DepMap/"
  maf.file = "OmicsSomaticMutations.rds"
  if (file.exists(paste0(depmap.path, depmap.ver, "/", maf.file))) {
    message("Processed cnv maf file found.")
    maf.read = readRDS(paste0(depmap.path, depmap.ver, "/", maf.file))
  } else {
    message("Processing the .csv file provided from the DepMap. May require manual fixing...")
    #test.read = read.table(paste0(depmap.path, depmap.ver, "/", maf.file), sep = ",", fill = NA, header = T)
    maf.csv.file = "OmicsSomaticMutations.csv"
    library(data.table)
    test.fread = fread(paste0(depmap.path, depmap.ver, "/", maf.csv.file))
    #testing if it could be saved as maf
    colnames(test.fread)[colnames(test.fread) == "HugoSymbol"] = "Hugo_Symbol"
    colnames(test.fread)[colnames(test.fread) == "Chrom"] = "Chromosome"
    colnames(test.fread)[colnames(test.fread) == "Pos"] = "Start_Position"
    #calculating end position
    test.fread$End_Position = test.fread$Start_Position + nchar(test.fread$Ref) -1
    colnames(test.fread)[colnames(test.fread) == "Pos"] = "Start_Position"
    colnames(test.fread)[colnames(test.fread) == "Ref"] = "Reference_Allele"
    test.fread$Tumor_Seq_Allele1 = test.fread$Reference_Allele
    colnames(test.fread)[colnames(test.fread) == "Alt"] = "Tumor_Seq_Allele2"
    colnames(test.fread)[colnames(test.fread) == "VariantInfo"] = "Variant_Classification"
    #fixing names
    test.fread$Variant_Classification[test.fread$Variant_Classification == "MISSENSE"] = "Missense_Mutation"
    test.fread$Variant_Classification[test.fread$Variant_Classification == "SILENT"] = "Silent"
    test.fread$Variant_Classification[test.fread$Variant_Classification == "IN_FRAME_INS"] = "In_Frame_Ins"
    test.fread$Variant_Classification[test.fread$Variant_Classification == "SPLICE_SITE"] = "Splice_Site"
    test.fread$Variant_Classification[test.fread$Variant_Classification == "NONSENSE"] = "Nonsense_Mutation"
    test.fread$Variant_Classification[test.fread$Variant_Classification == "FRAME_SHIFT_DEL"] = "Frame_Shift_Del"
    test.fread$Variant_Classification[test.fread$Variant_Classification == "NONSTOP"] = "Nonstop_Mutation"
    test.fread$Variant_Classification[test.fread$Variant_Classification == "START_CODON_SNP"] = "Translation_Start_Site"
    test.fread$Variant_Classification[test.fread$Variant_Classification == "IN_FRAME_DEL"] = "In_Frame_Del"
    test.fread$Variant_Classification[test.fread$Variant_Classification == "FRAME_SHIFT_INS"] = "Frame_Shift_Ins"
    test.fread$Variant_Classification[test.fread$Variant_Classification == "START_CODON_INS"] = "Translation_Start_Site"
    test.fread$Variant_Classification[test.fread$Variant_Classification == "FIVE_PRIME_FLANK"] = "5'Flank"
    test.fread$Variant_Classification[test.fread$Variant_Classification == "INTRON"] = "Intron"
    test.fread$Variant_Classification[test.fread$Variant_Classification == "THREE_PRIME_UTR"] = "3'UTR"
    colnames(test.fread)[colnames(test.fread) == "VariantType"] = "Variant_Type"
    colnames(test.fread)[colnames(test.fread) == "ModelID"] = "Tumor_Sample_Barcode"
    #test number
    maf.read = read.maf(test.fread)
    saveRDS(maf.read, paste0(depmap.path, depmap.ver, "/", maf.file))
  }
  tumor.barc = as.character(maf.read@variants.per.sample$Tumor_Sample_Barcode)
  sample.check = intersect(tumor.barc, sample.name)
  if (length(sample.check) < 2) {
    stop("Less than 2 samples found. Cannot proceed")
  } else {
    message(length(sample.check), " samples found. Processing...")
  }
  maf.sub = subsetMaf(maf.read, tsb = sample.check)
  output = mutCountMatrix(maf.sub)
  return(output)
}

filepath = "/home/seokwon/nas/"
ref_path = paste0(filepath, "99.reference/")

sce_depmap = readRDS("/mnt/gluster_server/data/reference/DepMap/23Q4/sce_depmap.rds")
exp_count_mat = assay(sce_depmap, 1)
# mut_count_mat = mutCountMatrix(mut)
exp_count_mat_tpm = (2^exp_count_mat)-1

meta_for_TCGA = read.csv(paste0(ref_path, "/Depmap_meta_filt_to_TCGA.csv"))
meta_tbl <- read.csv(file = '/mnt/gluster_server/data/reference/DepMap/23Q4/Model.csv', 
                     sep=',', header=T, fill=T)


Cancerlist = dir(paste0(filepath, "/00.data/"))
rdpn_path = "/mnt/gluster_server/06.seokwon/"

ppi.date = "20220117"

num.edge.multi = 10 #or 100 (for now)
random.num = 150
sourcepath = paste0(ref_path, "RDPN/",c(num.edge.multi),"X" ) 

Cancerlist = dir(paste0(filepath, "/00.data/filtered_TCGA/"))
# num_CancerType = "04.TCGA-CESC"
######################################################################################################################
## Network propagation using 'ppi_backbone_*_RDPN_#.rds' 
######################################################################################################################

for (num_CancerType in Cancerlist) {
  
  main.path = paste0(filepath, "/00.data/filtered_TCGA/", num_CancerType, "/")
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  
  # rdpn_main.path = paste0(rdpn_path, "00.data/", num_CancerType, "/")
  result_savepath = paste0(main.path, "03.RDPN_depmap/RDPN_",c(num.edge.multi),"X/")
  
  ######################################################################################################################
  ## check the exsited folder ##########################################################################################
  ######################################################################################################################
  
  ck.file = list.files(path =  result_savepath)
  ck.file = grep('depmap_cellline', ck.file)
  ck.rdpn.file = list.files(path =  paste0(main.path, "03.RDPN_depmap/RDPN_",c(num.edge.multi),"X/"))
  
  if (length(ck.file) == 1) {
    print("RDPN finish")
    next
    # print("1111")
    
  } else if (dir.exists(paste0(main.path, "03.RDPN_depmap/RDPN_",c(num.edge.multi),"X/")) && length(ck.rdpn.file) != 150) {
    
    print("it stopped")
    
    cell_oi <- meta_tbl %>% 
      subset(GrowthPattern != 'Organoid' & 
               PrimaryOrMetastasis == 'Primary' &
               OncotreeLineage == unique(meta_for_TCGA[which(meta_for_TCGA$TCGA == CancerType),]$oncotree) &
               OncotreeCode %in% unique(meta_for_TCGA[which(meta_for_TCGA$TCGA == CancerType),]$oncotreecode)) %>%
      dplyr::select(ModelID) %>% unlist(use.names=F)
    
    if (length(cell_oi) != 0) {
      exP_count_mat_tpm_filt = exp_count_mat_tpm %>%
        as.data.frame() %>%
        subset(select = names(.) %in% cell_oi)
      
      mut_count_mat = dep.snv.caller(sample.name = cell_oi)
      
      mut_count_mat_filt = mut_count_mat %>%
        as.data.frame() %>%
        subset(select = names(.) %in% cell_oi)
    } else {
      next
    }
    
    # input data (exp data and mut data)
    mut.mat = mut_count_mat_filt
    exp.log.mat = exP_count_mat_tpm_filt
    raw.prop.score = readRDS(paste0(main.path, dir(path = main.path , pattern = "net_prop_depmap_total")))
    query.genes = rownames(raw.prop.score) #12902 
    savepath = paste0(main.path, "03.RDPN_depmap/RDPN_",c(num.edge.multi),"X/") #mine RDPN ppi backbone based on the new ppi backbone 
    result_savepath = paste0(main.path, "03.RDPN_depmap/RDPN_",c(num.edge.multi),"X/")
    
    for (i in length(ck.rdpn.file):random.num) {  # when there are 50 iterations 
      if (length(ck.rdpn.file) == 150) {
        next
      }
      g.rdpn = readRDS(paste0(sourcepath,"/ppi_backbone_",ppi.date,"_RDPN",c(i),".rds")) #13064 614246
      g.rdpn.sub = induced_subgraph(g.rdpn, query.genes) #12974 613673 (node, edge)
      c.rdpn = components(g.rdpn.sub) #10253 (maximum) # c.rdpn$membership[10130] differ from the guery.names
      
      query.names = names(which(c.rdpn$membership == which(c.rdpn$csize == max(c.rdpn$csize)))) #10252
      g.rdpn.gs = induced_subgraph(g.rdpn, query.names) #12973 613673
      
      #fixing exp.mat 
      idx.mat = match(query.names, rownames(exp.log.mat)) #10252
      exp.rdpn = exp.log.mat[idx.mat, ] #10252 152 
      
      #fixing mut.mat 
      common.mut = intersect(query.names, rownames(mut.mat)) 
      idx.mut = match(common.mut, rownames(mut.mat))
      mut.rdpn = as.matrix(mut.mat[idx.mut,]) # 7244 152
      
      #fixing the samples for mut.mat and exp.mat 
      #mut.col = do.call("rbind", strsplit(colnames(mut.rdpn),"-"))
      #n.mut.col = paste(mut.col[,1],mut.col[,2],mut.col[,3],sep="-")
      colnames(mut.rdpn) = substr(colnames(mut.rdpn), 1, 16)
      
      #exp.col = do.call("rbind", strsplit(colnames(exp.rdpn),"-")) 
      #n.exp.col = paste(exp.col[,1],exp.col[,2],exp.col[,3],sep="-")
      #colnames(exp.rdpn) = n.exp.col 
      colnames(exp.rdpn) = substr(colnames(exp.rdpn), 1, 16)
      
      int.col.names = intersect(colnames(mut.rdpn), colnames(exp.rdpn)) #151
      
      if (length(int.col.names)==length(colnames(exp.rdpn))){
        idx.col.mut = match(int.col.names, colnames(mut.rdpn)) 
      } else {
        print("IDs are not exactly matched.....please check ! ")
      }
      mut.rdpn = mut.rdpn[,idx.col.mut] #11012 973
      
      #added 2021-Dec-21
      idx.col.exp = match(int.col.names, colnames(exp.rdpn))
      exp.rdpn = exp.rdpn[,idx.col.exp]
      
      #stretching 
      idx.stretch = match(rownames(mut.rdpn), rownames(exp.rdpn)) 
      n.mut.stretch = matrix(data = 0, nrow = nrow(exp.rdpn), ncol(exp.rdpn)) #12973 973
      colnames(n.mut.stretch) = colnames(exp.rdpn) 
      rownames(n.mut.stretch) = rownames(exp.rdpn) 
      
      #binarizing and stretching the mutation matrix 
      for (j in 1:ncol(mut.rdpn)) {
        query = mut.rdpn[,j]
        query[query > 0] = 1 
        n.mut.stretch[idx.stretch,j] = query 
      }
      
      idx.empty = colSums(n.mut.stretch) > 1  #check length(whihc(idx.empty == TRUE))
      n.mut.stretch = n.mut.stretch[,idx.empty] # 972 1 sample was dropped out 
      
      exp.rdpn = exp.rdpn[,idx.empty]
      sn = 1e-14 
      exp.rdpn[exp.rdpn == 0] = sn 
      
      #sample fixing for a single sample testing 
      #idx.sample = match(n.col.net, colnames(exp.rdpn)) 
      #exp.rdpn = as.matrix(exp.rdpn[,idx.sample])
      #n.mut.stretch = as.matrix(n.mut.stretch[,idx.sample])
      
      #testing stuff for compressed matrix 
      adj.rdpn = as.matrix(as_adjacency_matrix(g.rdpn.gs)) #12973 12973 
      
      N.pat <- dim(exp.rdpn)[2] #972  
      N.net <- dim(adj.rdpn)[1] #12973 
      alphav <- 0.7 
      
      # when not using parallel processing 
      #source("/home/jisushin/project/99.scripts/func_netprop.r")
      #for(j in 1:N.pat) {
      #    prop.res.new = net.propagation(exp.rdpn[,j], adj.rdpn, n.mut.stretch[,j], alphav)
      #    if(j == 1){
      #        prop.res.all = prop.res.new 
      #    } else {
      #        prop.res.all = cbind(prop.res.all, prop.res.new)
      #    }
      #}
      
      # for doing parallel processing
      cl <- makeCluster(20) 
      registerDoParallel(cl) 
      
      tic() 
      # .combine = > combine mode of results (it must be r function)
      prop.res.all = foreach(j=1:N.pat, .combine =cbind) %dopar% {
        source(paste0(ref_path,"func_netprop.r"))
        prop.res.new <- net.propagation(exp.rdpn[,j], adj.rdpn, n.mut.stretch[,j], alphav) 
        prop.res.new       
      }
      toc() 
      stopCluster(cl) 
      
      rownames(prop.res.all) = rownames(exp.rdpn) 
      colnames(prop.res.all) = colnames(exp.rdpn) 
      
      
      saveRDS(prop.res.all, file = paste0(savepath,"net_prop_RDPN",c(i),"_",ppi.date,".rds"))
      
    }
    
    ######################################################################################################################
    ## calculating the p-values by comparing the original and generated with the random degree preserving network 
    ######################################################################################################################
    
    rdpn.path = savepath
    
    #raw.net.file = list.files(ppi.path, pattern = ppi.date)
    raw.net.file = list.files(main.path, pattern ='net_prop_total')
    raw.prop.score = readRDS(paste0(main.path, raw.net.file))
    
    # mut.path = paste0("/home/jisushin/project/", cancer.type, "/00.raw/mutations/")
    # mut.file = list.files(mut.path, pattern = "mut_count_filt_data.rds") 
    # mut.mat = readRDS(paste0(mut.path,mut.file))
    
    pat.ids.test = colnames(raw.prop.score)
    colnames(raw.prop.score) = substr(pat.ids.test,1,16)
    
    colnames(mut.mat) = substr(colnames(mut.mat),1,16)
    
    count.mat = matrix(data = 0, nrow = nrow(raw.prop.score), ncol = ncol(raw.prop.score))
    rownames(count.mat) = rownames(raw.prop.score) 
    
    rand.mat = matrix(data = 0, nrow = nrow(raw.prop.score), ncol = ncol(raw.prop.score)) 
    rownames(rand.mat) = rownames(raw.prop.score) 
    
    tic() 
    
    for (i in 1:random.num) {
      rdpn.file = paste0(rdpn.path, "net_prop_RDPN", c(i), "_", ppi.date, ".rds")
      rand.prop.score = readRDS(rdpn.file)
      
      idx = match(rownames(rand.prop.score), rownames(raw.prop.score)) 
      idx.ID = match(colnames(rand.prop.score), colnames(raw.prop.score)) 
      
      count.mat[idx, idx.ID] = count.mat[idx, idx.ID]+1  #count up when gene is present 
      norm.score = rand.prop.score / raw.prop.score[idx, idx.ID]
      norm.logi = norm.score >= 1
      rand.mat[idx, idx.ID] = rand.mat[idx, idx.ID] + norm.logi 
    }
    toc() 
    
    pval.mat = rand.mat / count.mat 
    colnames(pval.mat) = colnames(raw.prop.score) 
    
    idx.mut.id = match(colnames(pval.mat), colnames(mut.mat)) 
    mut.mat = mut.mat[,idx.mut.id] 
    
    for (j in 1:ncol(pval.mat)) {
      query = pval.mat[,j] 
      mut.col.idx = match(colnames(pval.mat)[j], colnames(mut.mat)) 
      query.mut = as.matrix(mut.mat)[,mut.col.idx] 
      
      mut.genes = names(query.mut[which(query.mut > 0)])
      common = intersect(mut.genes, names(query)) 
      idx = match(common, names(query)) 
      pval.mat[idx, j] = 0 # to treat as very significant when the mutation exists 
    }
    
    saveRDS(pval.mat, file = paste0(result_savepath, "net_prop_pval_", random.num, "_depmap_cellline_total_samples_", ppi.date, ".rds"))
    
  } else {
    
    # create folder
    dir.create(paste0(main.path,"03.RDPN_depmap/"))
    dir.create(paste0(main.path,"03.RDPN_depmap/RDPN_",c(num.edge.multi),"X/"))
    # dir.create(paste0(rdpn_main.path,"03.RDPN_depmap/"))
    # dir.create(paste0(rdpn_main.path,"03.RDPN_depmap/RDPN_",c(num.edge.multi),"X/"))
    
    cell_oi <- meta_tbl %>% 
      subset(GrowthPattern != 'Organoid' & 
               PrimaryOrMetastasis == 'Primary' &
               OncotreeLineage == unique(meta_for_TCGA[which(meta_for_TCGA$TCGA == CancerType),]$oncotree) &
               OncotreeCode %in% unique(meta_for_TCGA[which(meta_for_TCGA$TCGA == CancerType),]$oncotreecode)) %>%
      dplyr::select(ModelID) %>% unlist(use.names=F)
    
    if (length(cell_oi) != 0) {
      exP_count_mat_tpm_filt = exp_count_mat_tpm %>%
        as.data.frame() %>%
        subset(select = names(.) %in% cell_oi)
      
      mut_count_mat = dep.snv.caller(sample.name = cell_oi)
      
      mut_count_mat_filt = mut_count_mat %>%
        as.data.frame() %>%
        subset(select = names(.) %in% cell_oi)
    } else {
      next
    }
    
    # input data (exp data and mut data)
    mut.mat = mut_count_mat_filt
    exp.log.mat = exP_count_mat_tpm_filt
    raw.prop.score = readRDS(paste0(main.path, dir(path = main.path , pattern = "net_prop_depmap_total")))
    query.genes = rownames(raw.prop.score) #12902 
    savepath = paste0(main.path, "03.RDPN_depmap/RDPN_",c(num.edge.multi),"X/") #mine RDPN ppi backbone based on the new ppi backbone 
    result_savepath = paste0(main.path, "03.RDPN_depmap/RDPN_",c(num.edge.multi),"X/")
    i = 1
    for (i in 1:random.num) {  # when there are 50 iterations 
      
      g.rdpn = readRDS(paste0(sourcepath,"/ppi_backbone_",ppi.date,"_RDPN",c(i),".rds")) #13064 614246
      g.rdpn.sub = induced_subgraph(g.rdpn, query.genes) #12974 613673 (node, edge)
      c.rdpn = components(g.rdpn.sub) #10253 (maximum) # c.rdpn$membership[10130] differ from the guery.names
      
      query.names = names(which(c.rdpn$membership == which(c.rdpn$csize == max(c.rdpn$csize)))) #10252
      g.rdpn.gs = induced_subgraph(g.rdpn, query.names) #12973 613673
      
      #fixing exp.mat 
      idx.mat = match(query.names, rownames(exp.log.mat)) #10252
      exp.rdpn = exp.log.mat[idx.mat, ] #10252 152 
      
      #fixing mut.mat 
      common.mut = intersect(query.names, rownames(mut.mat)) 
      idx.mut = match(common.mut, rownames(mut.mat))
      mut.rdpn = as.matrix(mut.mat[idx.mut,]) # 7244 152
      
      #fixing the samples for mut.mat and exp.mat 
      #mut.col = do.call("rbind", strsplit(colnames(mut.rdpn),"-"))
      #n.mut.col = paste(mut.col[,1],mut.col[,2],mut.col[,3],sep="-")
      colnames(mut.rdpn) = substr(colnames(mut.rdpn), 1, 16)
      
      #exp.col = do.call("rbind", strsplit(colnames(exp.rdpn),"-")) 
      #n.exp.col = paste(exp.col[,1],exp.col[,2],exp.col[,3],sep="-")
      #colnames(exp.rdpn) = n.exp.col 
      colnames(exp.rdpn) = substr(colnames(exp.rdpn), 1, 16)
      
      int.col.names = intersect(colnames(mut.rdpn), colnames(exp.rdpn)) #151
      
      if (length(int.col.names)==length(colnames(exp.rdpn))){
        idx.col.mut = match(int.col.names, colnames(mut.rdpn)) 
      } else {
        print("IDs are not exactly matched.....please check ! ")
      }
      mut.rdpn = mut.rdpn[,idx.col.mut] #11012 973
      
      #added 2021-Dec-21
      idx.col.exp = match(int.col.names, colnames(exp.rdpn))
      exp.rdpn = exp.rdpn[,idx.col.exp]
      
      #stretching 
      idx.stretch = match(rownames(mut.rdpn), rownames(exp.rdpn)) 
      n.mut.stretch = matrix(data = 0, nrow = nrow(exp.rdpn), ncol(exp.rdpn)) #12973 973
      colnames(n.mut.stretch) = colnames(exp.rdpn) 
      rownames(n.mut.stretch) = rownames(exp.rdpn) 
      
      #binarizing and stretching the mutation matrix 
      for (j in 1:ncol(mut.rdpn)) {
        query = mut.rdpn[,j]
        query[query > 0] = 1 
        n.mut.stretch[idx.stretch,j] = query 
      }
      
      idx.empty = colSums(n.mut.stretch) > 1  #check length(whihc(idx.empty == TRUE))
      n.mut.stretch = n.mut.stretch[,idx.empty] # 972 1 sample was dropped out 
      
      exp.rdpn = exp.rdpn[,idx.empty]
      sn = 1e-14 
      exp.rdpn[exp.rdpn == 0] = sn 
      
      #sample fixing for a single sample testing 
      #idx.sample = match(n.col.net, colnames(exp.rdpn)) 
      #exp.rdpn = as.matrix(exp.rdpn[,idx.sample])
      #n.mut.stretch = as.matrix(n.mut.stretch[,idx.sample])
      
      #testing stuff for compressed matrix 
      adj.rdpn = as.matrix(as_adjacency_matrix(g.rdpn.gs)) #12973 12973 
      
      N.pat <- dim(exp.rdpn)[2] #972  
      N.net <- dim(adj.rdpn)[1] #12973 
      alphav <- 0.7 
      
      # when not using parallel processing 
      #source("/home/jisushin/project/99.scripts/func_netprop.r")
      #for(j in 1:N.pat) {
      #    prop.res.new = net.propagation(exp.rdpn[,j], adj.rdpn, n.mut.stretch[,j], alphav)
      #    if(j == 1){
      #        prop.res.all = prop.res.new 
      #    } else {
      #        prop.res.all = cbind(prop.res.all, prop.res.new)
      #    }
      #}
      
      # for doing parallel processing
      cl <- makeCluster(20) 
      registerDoParallel(cl) 
      
      tic() 
      # .combine = > combine mode of results (it must be r function)
      prop.res.all = foreach(j=1:N.pat, .combine =cbind) %dopar% {
        source(paste0(ref_path,"func_netprop.r"))
        prop.res.new <- net.propagation(exp.rdpn[,j], adj.rdpn, n.mut.stretch[,j], alphav) 
        prop.res.new       
      }
      toc() 
      stopCluster(cl) 
      
      rownames(prop.res.all) = rownames(exp.rdpn) 
      colnames(prop.res.all) = colnames(exp.rdpn) 
      
      
      saveRDS(prop.res.all, file = paste0(savepath,"net_prop_RDPN",c(i),"_",ppi.date,".rds"))
      
    }
    
    
    ######################################################################################################################
    ## calculating the p-values by comparing the original and generated with the random degree preserving network 
    ######################################################################################################################
    
    rdpn.path = savepath
    
    #raw.net.file = list.files(ppi.path, pattern = ppi.date)
    raw.net.file = list.files(main.path, pattern ='net_prop_depmap_total')
    raw.prop.score = readRDS(paste0(main.path, raw.net.file))
    
    # mut.path = paste0("/home/jisushin/project/", cancer.type, "/00.raw/mutations/")
    # mut.file = list.files(mut.path, pattern = "mut_count_filt_data.rds") 
    # mut.mat = readRDS(paste0(mut.path,mut.file))
    
    pat.ids.test = colnames(raw.prop.score)
    colnames(raw.prop.score) = substr(pat.ids.test,1,16)
    
    colnames(mut.mat) = substr(colnames(mut.mat),1,16)
    
    count.mat = matrix(data = 0, nrow = nrow(raw.prop.score), ncol = ncol(raw.prop.score))
    rownames(count.mat) = rownames(raw.prop.score) 
    
    rand.mat = matrix(data = 0, nrow = nrow(raw.prop.score), ncol = ncol(raw.prop.score)) 
    rownames(rand.mat) = rownames(raw.prop.score) 
    
    tic() 
    
    for (i in 1:random.num) {
      print (i)
      rdpn.file = paste0(rdpn.path, "net_prop_RDPN", c(i), "_", ppi.date, ".rds")
      rand.prop.score = readRDS(rdpn.file)
      
      idx = match(rownames(rand.prop.score), rownames(raw.prop.score)) 
      idx.ID = match(colnames(rand.prop.score), colnames(raw.prop.score)) 
      
      count.mat[idx, idx.ID] = count.mat[idx, idx.ID]+1  #count up when gene is present 
      norm.score = rand.prop.score / raw.prop.score[idx, idx.ID]
      norm.logi = norm.score >= 1
      rand.mat[idx, idx.ID] = rand.mat[idx, idx.ID] + norm.logi 
    }
    toc() 
    
    pval.mat = rand.mat / count.mat 
    colnames(pval.mat) = colnames(raw.prop.score) 
    
    idx.mut.id = match(colnames(pval.mat), colnames(mut.mat)) 
    mut.mat = mut.mat[,idx.mut.id] 
    
    for (j in 1:ncol(pval.mat)) {
      query = pval.mat[,j] 
      mut.col.idx = match(colnames(pval.mat)[j], colnames(mut.mat)) 
      query.mut = as.matrix(mut.mat)[,mut.col.idx] 
      
      mut.genes = names(query.mut[which(query.mut > 0)])
      common = intersect(mut.genes, names(query)) 
      idx = match(common, names(query)) 
      pval.mat[idx, j] = 0 # to treat as very significant when the mutation exists 
    }
    
    saveRDS(pval.mat, file = paste0(result_savepath, "net_prop_pval_", random.num, "_depmap_cellline_total_samples_", ppi.date, ".rds"))
    
  }
  
} 





