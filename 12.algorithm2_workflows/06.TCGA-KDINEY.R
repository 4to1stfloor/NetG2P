
## Date: 2021.12.01 / 2022.07.20
## Writer: Jisu Shin 
## edit : Seok-Won Jang

##########################################################################################################################
## network propagation 
##########################################################################################################################

library(igraph)
library(tictoc) 
library(pcg)

filepath = setwd("/home/seokwon")
ref_path = paste0(filepath, "/99.reference/")
num_CancerType = "35.TCGA-KIDNEY"
main.path = paste0(filepath, "/00.data/", num_CancerType, "/")
CancerType = "KIDNEY"

# input data (exp data and mut data)
mut.mat = readRDS(paste0(main.path,CancerType,"_mut_count_filt_data.rds")) 
exp.log.mat = readRDS(paste0(main.path,CancerType,"_exp_TPM_mat_filt_log_gene_symbol.rds"))

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
  
  
  saveRDS(prop.res.all, file = paste0(main.path,"net_prop_total_", N.pat, ".rds"))
  
  ######################################################################################################################
  
} 

## 03.RDPN
## data: 2022.01.06 
## Writer: Jisu Shin 
## edit : Seok-Won Jang

## /home/seokwon -> there is not a empty space anymore
## mv RDPN to /mnt/gluster_server/06.seokwon

## This script is for running RDPN (random degree preserving network)  with new PPI backbone 

#loading tools 
require(igraph) 
require(data.table)
require(foreach) 
require(doParallel) 
library(tictoc)

# setting

filepath = setwd("/home/seokwon")
ref_path = paste0(filepath, "/99.reference/")
Cancerlist = dir(paste0(filepath, "/00.data/"))
rdpn_path = "/mnt/gluster_server/06.seokwon/"

ppi.date = "20220117"

num.edge.multi = 10 #or 100 (for now)
random.num = 150
sourcepath = paste0(ref_path, "RDPN/",c(num.edge.multi),"X" ) 



######################################################################################################################
## Network propagation using 'ppi_backbone_*_RDPN_#.rds' 
######################################################################################################################



main.path = paste0(filepath, "/00.data/", num_CancerType, "/")
rdpn_main.path = paste0(rdpn_path, "00.data/", num_CancerType, "/")
result_savepath = paste0(main.path, "03.RDPN/RDPN_",c(num.edge.multi),"X/")

CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))

######################################################################################################################
## check the exsited folder ##########################################################################################
######################################################################################################################

ck.file = list.files(path =  result_savepath)
ck.rdpn.file = list.files(path =  paste0(rdpn_main.path, "03.RDPN/RDPN_",c(num.edge.multi),"X/"))

if (length(ck.file) == 1) {
  print("RDPN finish")
  next
  # print("1111")
  
  } else if (dir.exists(paste0(rdpn_main.path, "03.RDPN/RDPN_",c(num.edge.multi),"X/"))&& length(ck.rdpn.file) !=0 && length(ck.rdpn.file) != 150) {
    print("it stopped")
  
    # input data (exp data and mut data)
    mut.mat = readRDS(paste0(main.path,CancerType,"_mut_count_filt_data.rds")) 
    exp.log.mat = readRDS(paste0(main.path,CancerType,"_exp_TPM_mat_filt_log_gene_symbol.rds"))
    raw.prop.score = readRDS(paste0(main.path,dir(path = main.path , pattern = "net_prop_total")))
    query.genes = rownames(raw.prop.score) #12902 
    savepath = paste0(rdpn_main.path, "03.RDPN/RDPN_",c(num.edge.multi),"X/") #mine RDPN ppi backbone based on the new ppi backbone 
    
    for (i in length(ck.file):random.num) {  # when there are 50 iterations 
      if (length(ck.file) == 150) {
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
    
    saveRDS(pval.mat, file = paste0(result_savepath, "net_prop_pval_", random.num, "rep_total_samples_", ppi.date, ".rds"))
    
    } else {
    
    # input data (exp data and mut data)
    mut.mat = readRDS(paste0(main.path,CancerType,"_mut_count_filt_data.rds")) 
    exp.log.mat = readRDS(paste0(main.path,CancerType,"_exp_TPM_mat_filt_log_gene_symbol.rds"))
    raw.prop.score = readRDS(paste0(main.path,dir(path = main.path , pattern = "net_prop_total")))
    query.genes = rownames(raw.prop.score) #12902 
    savepath = paste0(rdpn_main.path, "03.RDPN/RDPN_",c(num.edge.multi),"X/") #mine RDPN ppi backbone based on the new ppi backbone 
    result_savepath = paste0(main.path, "03.RDPN/RDPN_",c(num.edge.multi),"X/")
    
    
    # create folder
    dir.create(paste0(main.path,"03.RDPN/"))
    dir.create(paste0(main.path,"03.RDPN/RDPN_",c(num.edge.multi),"X/"))
    dir.create(paste0(rdpn_main.path,"03.RDPN/"))
    dir.create(paste0(rdpn_main.path,"03.RDPN/RDPN_",c(num.edge.multi),"X/"))
    
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
    
    saveRDS(pval.mat, file = paste0(result_savepath, "net_prop_pval_", random.num, "rep_total_samples_", ppi.date, ".rds"))
    
    }

## Date: 2022-03-28
## Writer :Jisu Shin 
## edit : Seok-Won Jang

#loading R packages
library(data.table)
library(igraph)
library(data.table)
library(survival)
library(survminer)

# global setting
filepath = setwd("/home/seokwon")
ref_path = paste0(filepath, "/99.reference/")
CancerType
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

###############################################################################################################################################
## This script is for hypergeometric test (both pathway and pathway link levels)
###########################################################################################################################################################################
### Direct to phenotype (Net.prop & RDPN --> GC --> 54 KEGG pathway links)
###########################################################################################################################################################################

## Function define

###########################################################################################################################################################################
### 1.function of pathway-level (e.g., P1, P2, P3, ..., P54) 
###########################################################################################################################################################################
GC.hyper.pathway.direct  <- function(CancerType, ppi.date, path.dat, path.file) {
  
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
  
  saveRDS(enr.mat, paste0(main.path, "phyper_enrichment_", length(unique(path$pathway)), "_", path.dat, "_GC_",ppi.date, ".rds"))
  
  # transforming the p-values obtained from the hypergeometric tests with log10(pvalue)
  #enr = readRDS(paste0("/home/jisushin/project/", cancer.type, "/phyper_enrichment_", chg.dat, "_10chg_GC_", ppi.date, "_grn.rds"))
  #enr.trans = -log10(enr) 
  #saveRDS(enr.trans, paste0("/home/jisushin/project/", cancer.type, "/phyper_enrichment_",chg.dat, "_10chg_GC_", ppi.date, "_log_trans_grn.rds"))
}
###########################################################################################################################################################################
### 2.function of pathway-link level (e.g., P1-P2, P3-P10, ..., )
###########################################################################################################################################################################
GC.hyper.path.net.direct  <- function(CancerType, ppi.date, path.dat, path.file) {
  
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
      first = unique(path$Symbol[which(path$pathway == pair[1])])
      second = unique(path$Symbol[which(path$pathway == pair[2])])
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
  
  saveRDS(dat, paste0(main.path, "phyper_enrichment_", length(unique(path$pathway)), "_", path.dat, "_net_GC_", ppi.date, ".rds"))
}
###########################################################################################################################################################################
### 3. function of obtaining p values from the survival analysis
###########################################################################################################################################################################
cal.surv.p <- function(CancerType, ppi.date, path.dat, level, filt.genes.n) { 
  
  hyper.mat = readRDS(paste0(main.path, "phyper_enrichment_", length(unique(path$pathway)), "_", path.dat, "_net_GC_", ppi.date, ".rds"))
  #print(hyper.mat)
  if (level != ""){
    sig.pair = dat$name[which(dat$n_genes > as.numeric(filt.genes.n))] 
    idx = match(sig.pair, rownames(hyper.mat))
    hyper.filt = hyper.mat[idx, ]
  } else {
    hyper.filt = hyper.mat
  }
  
  #exp.log.mat = readRDS(paste0(CancerType, '_exp_TPM_mat_filt_log_gene_symbol.rds')) #56404 973 
  mut.mat = readRDS(paste0(main.path, CancerType ,'_mut_count_filt_data.rds'))#16131 985
  mut.count = as.data.frame(colSums(mut.mat))
  colnames(mut.count) = "counts" 
  mut.count$patient_id = substr(rownames(mut.count), 1, 16) 
  mut.filt = mut.count[which(mut.count$counts < 1000), ]
  int.id = intersect(colnames(hyper.filt), mut.filt$patient_id) 
  idx.id = match(int.id, colnames(hyper.filt))
  hyper.filt = hyper.filt[, idx.id ]
  
  #loading clinical information for survival analysis 
  # Patient ID -> submitter_id, Overall_survival -> vital_status , Overall_survival(Month) -> days_to_death or days_to_last_follow_up
  # time -> overall_survival  
  
  cli = fread(paste0(ref_path,"all_clin_indexed.csv")) 
  
  cli_surv1 = cli[cli$project == "TCGA-KICH",
                  c("submitter_id",
                    "vital_status",
                    "days_to_death",
                    "days_to_last_follow_up")]
  
  cli_surv2 = cli[cli$project == "TCGA-KIRC",
                  c("submitter_id",
                    "vital_status",
                    "days_to_death",
                    "days_to_last_follow_up")]
  cli_surv3 = cli[cli$project == "TCGA-KIRP",
                  c("submitter_id",
                    "vital_status",
                    "days_to_death",
                    "days_to_last_follow_up")]
  
  cli_surv = rbind(cli_surv1,cli_surv2,cli_surv3)
  
  cli_surv$deceased = cli_surv$vital_status == "Dead"
  cli_surv$overall_survival = ifelse(cli_surv$deceased,
                                     cli_surv$days_to_death,
                                     cli_surv$days_to_last_follow_up)
  
  cli_surv = cli_surv[!is.na(cli_surv$vital_status),]
  
  cli_surv$status = NA
  cli_surv$status[which(cli_surv$vital_status == "Dead")] = 1
  cli_surv$status[which(cli_surv$vital_status == "Alive")] = 0
  
  cli_surv$status = as.numeric(cli_surv$status)
  cli_surv
  
  wrong_cluster = vector()
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
    # print (cli_surv_filt)
    
    
    if (length(unique(cli_surv_filt$cluster)) == 1) {
      wrong_cluster = c(wrong_cluster, p) 
      next
    }
    
    fit = survminer::surv_fit(Surv(overall_survival,status) ~ cluster, data = cli_surv_filt) 
    #fit = survfit(Surv(overall_survival,status) ~ cluster, data = cli_surv_filt)
    # print(p)
    pval.each = surv_pvalue(fit)[2]
    # print(p)
    if( p == 1) { 
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
    saveRDS(pval.com, paste0(main.path, "/surv_pvals_54_", path.dat, level, "_GC_", ppi.date, "_filt", filt.genes.n, ".rds"))
  } else {
    saveRDS(pval.com, paste0(main.path,"/surv_pvals_54_", path.dat, level, "_GC_", ppi.date, ".rds"))
  }
  
  wrong_cluster = vector()
}
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# 4. function of obtaining p values from the survival analysis (for each pathway )
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

cal.surv.each.p <- function(CancerType, ppi.date, path.dat) { 
  
  hyper.mat = readRDS(paste0(main.path, "phyper_enrichment_", length(unique(path$pathway)), "_", path.dat, "_GC_", ppi.date, ".rds"))
  hyper.filt = hyper.mat
  
  #exp.log.mat = readRDS(paste0(CancerType, '_exp_TPM_mat_filt_log_gene_symbol.rds')) #56404 973 
  mut.mat = readRDS(paste0(main.path, CancerType ,'_mut_count_filt_data.rds'))#16131 985
  mut.count = as.data.frame(colSums(mut.mat))
  colnames(mut.count) = "counts" 
  mut.count$patient_id = substr(rownames(mut.count) , 1 ,16)
  mut.filt = mut.count[which(mut.count$counts < 1000), ]
  int.id = intersect(colnames(hyper.filt), mut.filt$patient_id) 
  idx.id = match(int.id, colnames(hyper.filt))
  hyper.filt = hyper.filt[, idx.id ]
  
  #loading clinical information for survival analysis 
  # Patient ID -> submitter_id, Overall_survival -> vital_status , Overall_survival(Month) -> days_to_death or days_to_last_follow_up
  # time -> overall_survival  
  
  cli = fread(paste0(ref_path,"all_clin_indexed.csv")) 
  cli
  cli_surv1 = cli[cli$project == "TCGA-KICH",
                  c("submitter_id",
                    "vital_status",
                    "days_to_death",
                    "days_to_last_follow_up")]
  
  cli_surv2 = cli[cli$project == "TCGA-KIRC",
                  c("submitter_id",
                    "vital_status",
                    "days_to_death",
                    "days_to_last_follow_up")]
  cli_surv3 = cli[cli$project == "TCGA-KIRP",
                  c("submitter_id",
                    "vital_status",
                    "days_to_death",
                    "days_to_last_follow_up")]
  
  cli_surv = rbind(cli_surv1,cli_surv2,cli_surv3)
  cli_surv$deceased = cli_surv$vital_status == "Dead"
  cli_surv$overall_survival = ifelse(cli_surv$deceased,
                                     cli_surv$days_to_death,
                                     cli_surv$days_to_last_follow_up)
  
  cli_surv = cli_surv[!is.na(cli_surv$vital_status),]
  
  cli_surv$status = NA
  cli_surv$status[which(cli_surv$vital_status == "Dead")] = 1
  cli_surv$status[which(cli_surv$vital_status == "Alive")] = 0
  
  cli_surv$status = as.numeric(cli_surv$status)
  cli_surv
  
  
  wrong_cluster = vector()
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
    # print (cli_surv_filt)
    
    
    if (length(unique(cli_surv_filt$cluster)) == 1) {
      wrong_cluster = c(wrong_cluster, p) 
      next
    }
    
    fit = survminer::surv_fit(Surv(overall_survival,status) ~ cluster, data = cli_surv_filt) 
    #fit = survfit(Surv(overall_survival,status) ~ cluster, data = cli_surv_filt)
    # print(p)
    pval.each = surv_pvalue(fit)[2]
    # print(p)
    if( p == 1) { 
      pval.com = as.numeric(pval.each)
    } else {
      pval.com = c(pval.com, as.numeric(pval.each))
      #print(pval.com)
    } 
  }
  
  pval.com = as.data.frame(pval.com)
  print(wrong_cluster)
  if (length(wrong_cluster) == 0) {
    rownames(pval.com) = rownames(hyper.filt)
  } else if (length(wrong_cluster) != 0 ) {
    rownames(pval.com) = rownames(hyper.filt)[-wrong_cluster] 
  } else {
    print("I dont know")
  }
  
  colnames(pval.com) = gsub(".*\\.", "", CancerType)
  # print(pval.com)
  
  saveRDS(pval.com, paste0(main.path,"/surv_pvals_54_", path.dat, "_GC_", ppi.date, ".rds"))
  
  wrong_cluster = vector()
}



#loading RDPN output and function 2 output 
pval.path = paste0(main.path, "/03.RDPN/RDPN_",c(num.edge.multi),"X/")
pval.mat = readRDS(paste0(pval.path,"net_prop_pval_",random.num ,"rep_total_samples_",ppi.date,".rds")) # output of 03.RDPN.r 

# Execute function 1
GC.hyper.pathway.direct(CancerType, ppi.date, path.dat, path.file)

# Execute function 2
GC.hyper.path.net.direct(CancerType, ppi.date, path.dat, path.file) ## then "phyper_enrichment_54_KEGG_net_GC_20220117.rds" will be saved

# Execute function 3
cal.surv.p(CancerType, ppi.date, path.dat, level, filt.genes.n ) # this will generate "surv_pvals_54_KEGG_net_GC_20220117_filt5.rds" 

# Execute function 4
cal.surv.each.p(CancerType, ppi.date, path.dat ) 




















