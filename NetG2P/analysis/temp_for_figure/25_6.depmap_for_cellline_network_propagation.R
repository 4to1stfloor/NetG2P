library(tidyverse)
library(dplyr)
library(igraph)
library(tictoc) 
library(pcg)
library(SingleCellExperiment)

filepath = "/home/seokwon/nas/"
ref_path = paste0(filepath, "99.reference/")

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

sce_depmap = readRDS("/mnt/gluster_server/data/reference/DepMap/23Q4/sce_depmap.rds")
# mut = readRDS("/mnt/gluster_server/data/reference/DepMap/23Q4/OmicsSomaticMutations.rds")
meta_for_TCGA = read.csv(paste0(ref_path, "/Depmap_meta_filt_to_TCGA.csv"))

exp_count_mat = assay(sce_depmap, 1)
# mut_count_mat = mutCountMatrix(mut)
exp_count_mat_tpm = (2^exp_count_mat)-1

meta_tbl <- read.csv(file = '/mnt/gluster_server/data/reference/DepMap/23Q4/Model.csv', 
                     sep=',', header=T, fill=T)

##########################################################################################################################
## network propagation 
##########################################################################################################################

Cancerlist = dir(paste0(filepath, "/00.data/filtered_TCGA/"))

for (num_CancerType in Cancerlist) {
  
  main.path = paste0(filepath, "/00.data/filtered_TCGA/", num_CancerType, "/")
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  
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
    
    # input data (exp data and mut data)
    exp.log.mat = exP_count_mat_tpm_filt
    mut.mat = mut_count_mat_filt
    
    # protein to protein 
    g.ppi.conn = readRDS(file = paste0(ref_path, "ppi_backbone_20220117.rds"))
    g.ppi.conn.simp = simplify(g.ppi.conn) # symobolic simplification of an expression or function 
    gene.backbone = names(V(g.ppi.conn.simp))
    
    gene.query = rownames(exp.log.mat) #doing with exp.log.mat. taking max stuff #60466 in TCGA-COAD 
    common_ppi_exp = intersect(gene.backbone, gene.query) # probably okay without filtering from exp.log. if error appears, do that. #13024 in TCGA-BRCA
    
    g.ppi.spec = induced.subgraph(g.ppi.conn, common_ppi_exp) #10257 
    c.ppi.spec = components(g.ppi.spec)
    idx.max = which(c.ppi.spec$csize == max(c.ppi.spec$csize)) ## clustered network groups, and then select the biggest network based on the cluster size 
    
    g.ppi.gc = induced.subgraph(g.ppi.spec, names(which(c.ppi.spec$membership == idx.max))) #12902 138948 
    n.adj.mat = as.matrix(as_adjacency_matrix(g.ppi.gc)) #as_adjacency_matrix function in igraph package, to create the adjacency matrix for undirected graphs 
    gene.adj = rownames(n.adj.mat) #less than common. because we are isolating the largest component
    
    common.exp = intersect(gene.adj, rownames(exp.log.mat)) #12902
    idx.exp = match(common.exp, rownames(exp.log.mat))
    n.log.mat = exp.log.mat[idx.exp,] #12902 973 
    
    #fixing the issue of abs 0s (0 value may sometimes cause errors while doing network propagation, insert the smallest value instead 0)
    sn =1e-14
    n.log.mat[n.log.mat == 0] = sn #length(n.log.mat[n.log.mat ==0]) 1015081 
    
    #fixing it 
    # first shirinking the mut.mat to contain only the common genes found in n.log.mat 
    common = intersect(rownames(mut.mat),rownames(n.log.mat)) # 9782 
    idx.shrink = match(common, rownames(mut.mat)) 
    mut.mat.shrink = mut.mat[idx.shrink, ] #9782 436 (#genes #patients)
    
    ################################################################################################################################
    ## add the process for matching the number of samples (if the samples are not matched yet in the previous preparation steps)  
    ################################################################################################################################
    #colnames(mut.mat.shrink) = substr(colnames(mut.mat.shrink),1,16) #973
    #colnames(n.log.mat) = substr(colnames(n.log.mat),1,16) #1222 
    
    #common_smp = intersect(colnames(mut.mat.shrink), colnames(n.log.mat)) #973
    
    #idx.smp = match(common_smp, colnames(n.log.mat)) 
    #n.exp.mat = n.log.mat[,idx.smp] #12902 973 
    
    #if this process is not need (when the number of samples are matched already)
    n.exp.mat = n.log.mat 
    
    ###################################################################################################
    # then streching it to the size of n.exp.mat (making mut.mat to have same genes with the expression mat) 
    idx.stretch = match(rownames(mut.mat.shrink), rownames(n.exp.mat))
    
    n.mut.stretch = matrix (data = 0 , nrow=nrow(n.exp.mat), ncol = ncol(n.exp.mat))  #12902 973 
    colnames(n.mut.stretch) = colnames(n.exp.mat) 
    rownames(n.mut.stretch) = rownames(n.exp.mat) 
    
    for ( i in 1:ncol(n.mut.stretch)) {
      query = mut.mat.shrink[,i]
      query[which(query > 0)] = 1 
      n.mut.stretch[idx.stretch,i]=query 
    }
    ##################################################################################################
    
    idx.muts = which(colSums(n.mut.stretch) > 1) # taking stuff with at least 2 muts #971 
    n.mut.stretch = n.mut.stretch[,idx.muts] #10253 436 
    n.exp.mat = n.exp.mat[,idx.muts] #10253 436 
    
    if (all.equal(rownames(n.exp.mat), rownames(n.mut.stretch)) && all.equal(rownames(n.exp.mat), rownames(n.adj.mat)) == TRUE) {
      ###################################################################################################
      ## Network propagation using parallel, depends on the sample size  
      ###################################################################################################
      N.pat <- dim(n.exp.mat)[2] #971 
      N.net <- dim(n.adj.mat)[1] #12902 
      alphav <- 0.7 # propagation parameter (0.3 is assumed as returning back to the starting point, other 0.7 is assumed as going to the end network, not returning)
      
      require(foreach)
      require(doParallel)
      
      # parallel::detectCores()-1 #detect the number of available cores 
      binlimit = 192 # 24 * 8. using 24 cores is 1.5 times faster than 36 cores due to memory issues (= 24 cores for 8 samples each time) 
      
      ## test
      if (N.pat < binlimit) {
        cl <- makeCluster(24, outfile = "log.txt") 
        registerDoParallel(cl) 
        tic() 
        prop.res.all = foreach (i=1:N.pat, .combine = cbind) %dopar% {
          source(paste0(ref_path,"func_netprop.r"))
          prop.res.new <- net.propagation(n.exp.mat[,i], n.adj.mat, n.mut.stretch[,i], alphav) 
          prop.res.new 
        }
        toc() 
        stopCluster(cl) 
      } else { 
        # subsetting data
        bin = ceiling(N.pat / binlimit) 
        message(paste("data is too large. Dividing into", c(bin), "bins"))
        cl <- makeCluster(24, outfile = "log.txt") 
        registerDoParallel(cl) 
        tic() 
        prop.res.all = foreach(i=1:binlimit, .combine = cbind) %dopar% { 
          source(paste0(ref_path,"func_netprop.r"))
          prop.res.new <- net.propagation(n.exp.mat[,i], n.adj.mat, n.mut.stretch[,i], alphav) 
          prop.res.new
        }
        toc() 
        stopCluster(cl) 
        for (j in 2:bin) {
          max.query = binlimit*j 
          if (max.query > N.pat) { 
            max.query = N.pat 
          }
          cl <- makeCluster(24, outfile = "log.txt")
          registerDoParallel(cl) 
          tic() 
          prop.res.add = foreach(i=((binlimit*(j-1))+1): max.query, .combine = cbind) %dopar% {
            source(paste0(ref_path,"func_netprop.r")) 
            prop.res.new <- net.propagation(n.exp.mat[,i], n.adj.mat, n.mut.stretch[,i], alphav) 
            prop.res.new 
          }
          toc() 
          stopCluster(cl) 
          prop.res.all = cbind(prop.res.all, prop.res.add) 
        }
      }
      
      rownames(prop.res.all) = rownames(n.exp.mat)
      colnames(prop.res.all) = colnames(n.exp.mat) 
      
      
      saveRDS(prop.res.all, file = paste0(main.path, CancerType,"_net_prop_depmap_total_", N.pat, ".rds"))
      
      ######################################################################################################################
      
    } 
  } else {
    next
  }
  
}











