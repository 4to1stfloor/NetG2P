require(biomaRt)
require(TCGAbiolinks)
require(SummarizedExperiment)
require(SingleCellExperiment)
require(maftools)
require(dplyr)
require(DT)
library(dplyr)

##########################################################################################################################
## 01. data prepare 
##########################################################################################################################
getwd()
kidney_path ="/home/seokwon/01.externaldata/kidney/"
CancerType_in = "rt_target_2018_pub"
kidney_indi_path = paste0(kidney_path, CancerType_in,"/")

# expretion matrix and convert form RPKM to TPM
kidney_rpkm = read.table(paste0(kidney_indi_path, "data_mrna_seq_rpkm.txt" ), header = T,fill = TRUE)
kidney_rpkm_filt = kidney_rpkm
while (sum(duplicated(kidney_rpkm_filt$Hugo_Symbol)) != 0 ) {
  print(sum(duplicated(kidney_rpkm_filt$Hugo_Symbol)))
  dupl = which(duplicated(kidney_rpkm_filt$Hugo_Symbol))[1]
  dupl_gene = kidney_rpkm_filt[dupl,]$Hugo_Symbol
  tmp.df = kidney_rpkm_filt[which(kidney_rpkm_filt$Hugo_Symbol == dupl_gene),]
  tmp.df.2 = tmp.df[,-1]
  tmp.df.2 = t(as.data.frame(sapply(tmp.df.2, mean)))
  rownames(tmp.df.2) = dupl
  tmp.df.2 = as.data.frame(tmp.df.2)
  tmp.df.2$Hugo_Symbol = tmp.df$Hugo_Symbol[1]
  # tmp.df.2$Entrez_Gene_Id = tmp.df$Entrez_Gene_Id[1]
  tmp.df.2= tmp.df.2 %>% relocate("Hugo_Symbol", .before = colnames(kidney_rpkm_filt)[1])
  # tmp.df.2= tmp.df.2 %>% relocate("Entrez_Gene_Id", .before = "CH79T")
  kidney_rpkm_filt = kidney_rpkm_filt[-which(kidney_rpkm_filt$Hugo_Symbol == dupl_gene),]
  kidney_rpkm_filt = rbind(kidney_rpkm_filt,tmp.df.2)
  tmp.df.2 = NULL
  tmp.df =NULL
}


# if rna-seq was done by paired-end, below code should execute.
sum(is.na(kidney_rpkm_filt))

kidney_rpkm_filt_tmp = kidney_rpkm_filt[,c(-1,-2)]
kidney_rpkm_filt_tmp = kidney_rpkm_filt_tmp/2
kidney_rpkm_filt_tmp[is.na(kidney_rpkm_filt_tmp)] = 0
fpkmToTpm <- function(fpkm)
{
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}

kidney_tpm_filt =sapply(kidney_rpkm_filt_tmp, fpkmToTpm)  
kidney_tpm_filt = as.data.frame(kidney_tpm_filt)
kidney_tpm_filt$Hugo_Symbol = kidney_rpkm_filt$Hugo_Symbol
# luad_tpm_filt$Entrez_Gene_Id = luad_rpkm_filt$Entrez_Gene_Id
kidney_tpm_filt= kidney_tpm_filt %>% relocate("Hugo_Symbol", .before =  colnames(kidney_tpm_filt)[1])
# eurprad_tpm_filt= eurprad_tpm_filt %>% relocate("Entrez_Gene_Id", .before = "CH79T")
rownames(kidney_tpm_filt) = kidney_tpm_filt$Hugo_Symbol
kidney_tpm_filt = kidney_tpm_filt[,-1]

# mutation count matrix
kidney_maf.data = read.maf(paste0(kidney_indi_path, "data_mutations.txt"))
mut.count.mat_kidney = mutCountMatrix(kidney_maf.data) 

colnames(kidney_tpm_filt) = gsub(pattern = "\\.", replacement = "-", x = colnames(kidney_tpm_filt))
kidney_tpm_filt = kidney_tpm_filt[,colnames(mut.count.mat_kidney)]

saveRDS(kidney_tpm_filt, file=paste0(kidney_indi_path, CancerType_in ,"_exp_TPM_mat_filt_log_gene_symbol.rds")) #60466 436 
saveRDS(mut.count.mat_kidney, file=paste0(kidney_indi_path , CancerType_in ,"_mut_count_filt_data.rds"))


##########################################################################################################################
## 02. network propagation 
##########################################################################################################################

library(igraph)
library(tictoc) 
library(pcg)

filepath = "/home/seokwon/"
ref_path = paste0(filepath, "99.reference/")

main.path = paste0(kidney_path,CancerType_in, "/")

# input data (exp data and mut data)
mut.mat = readRDS(paste0(main.path,CancerType_in,"_mut_count_filt_data.rds")) 
exp.log.mat = readRDS(paste0(main.path,CancerType_in,"_exp_TPM_mat_filt_log_gene_symbol.rds"))

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
mut.mat.shrink = mut.mat[idx.shrink, ] #6143  106 (#genes #patients)

#
# add the process for matching the number of samples (if the samples are not matched yet in the previous preparation steps)  
#
#colnames(mut.mat.shrink) = substr(colnames(mut.mat.shrink),1,16) #973
#colnames(n.log.mat) = substr(colnames(n.log.mat),1,16) #1222 

#common_smp = intersect(colnames(mut.mat.shrink), colnames(n.log.mat)) #973

#idx.smp = match(common_smp, colnames(n.log.mat)) 
#n.exp.mat = n.log.mat[,idx.smp] #12902 973 

#if this process is not need (when the number of samples are matched already)
n.exp.mat = n.log.mat 

#
# then streching it to the size of n.exp.mat (making mut.mat to have same genes with the expression mat) 
idx.stretch = match(rownames(mut.mat.shrink), rownames(n.exp.mat))

n.mut.stretch = matrix (data = 0 , nrow=nrow(n.exp.mat), ncol = ncol(n.exp.mat))  #7887  106
colnames(n.mut.stretch) = colnames(n.exp.mat) 
rownames(n.mut.stretch) = rownames(n.exp.mat) 

for ( i in 1:ncol(n.mut.stretch)) {
  query = mut.mat.shrink[,i]
  query[which(query > 0)] = 1 
  n.mut.stretch[idx.stretch,i]=query 
}

idx.muts = which(colSums(n.mut.stretch) > 1) # taking stuff with at least 2 muts #971 
n.mut.stretch = n.mut.stretch[,idx.muts] #7887  106 
n.exp.mat = n.exp.mat[,idx.muts] #7887  106  

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
##########################################################################################################################
## 03. RDPN 
########################################################################################################################## 

#loading tools 
require(igraph) 
require(data.table)
require(foreach) 
require(doParallel) 
library(tictoc)

# setting

rdpn_main.path = "/mnt/gluster_server/06.seokwon/01.externaldata/"
num.edge.multi = 10
sourcepath = paste0(ref_path, "RDPN/",c(num.edge.multi),"X" ) 
random.num = 150
ppi.date = 20220117
result_savepath = paste0(main.path,"03.RDPN/RDPN_",c(num.edge.multi),"X/")
#
# Network propagation using 'ppi_backbone_*_RDPN_#.rds' 
#

# input data (exp data and mut data)
mut.mat = readRDS(paste0(main.path,CancerType_in,"_mut_count_filt_data.rds")) 
exp.log.mat = readRDS(paste0(main.path,CancerType_in,"_exp_TPM_mat_filt_log_gene_symbol.rds"))
raw.prop.score = readRDS(paste0(main.path,dir(path = main.path , pattern = "net_prop_total")))
query.genes = rownames(raw.prop.score) #12902 
savepath = paste0(rdpn_main.path,"kidney/", CancerType_in , "/03.RDPN/RDPN_",c(num.edge.multi),"X/") #mine RDPN ppi backbone based on the new ppi backbone 

# create folder
dir.create(paste0(main.path,"03.RDPN/"))
dir.create(paste0(main.path,"03.RDPN/RDPN_",c(num.edge.multi),"X/"))
dir.create(paste0(rdpn_main.path,"kidney/", CancerType_in,"/03.RDPN"))
dir.create(paste0(rdpn_main.path,"kidney/", CancerType_in, "/03.RDPN/RDPN_",c(num.edge.multi),"X/"))


#
# check the exsited folder 
#

ck.file = list.files(path =  result_savepath)
ck.rdpn.file = list.files(path =  paste0(rdpn_main.path, "kidney/", CancerType_in, "/03.RDPN/RDPN_",c(num.edge.multi),"X/"))

if (length(ck.file) == 1) {
  print("RDPN finish")
  next
  # print("1111")
} else if (dir.exists(paste0(rdpn_main.path, "kidney/", CancerType_in, "/03.RDPN/RDPN_",c(num.edge.multi),"X/")) && length(ck.rdpn.file) != 0 && length(ck.rdpn.file) != 150) {
  print("it stopped")
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
    # colnames(mut.rdpn) = substr(colnames(mut.rdpn), 1, 16)
    
    #exp.col = do.call("rbind", strsplit(colnames(exp.rdpn),"-")) 
    #n.exp.col = paste(exp.col[,1],exp.col[,2],exp.col[,3],sep="-")
    #colnames(exp.rdpn) = n.exp.col 
    # colnames(exp.rdpn) = substr(colnames(exp.rdpn), 1, 16)
    
    dim(mut.rdpn)
    dim(exp.rdpn)
    
    int.col.names = intersect(colnames(mut.rdpn), colnames(exp.rdpn)) #151
    
    idx.col.mut = match(int.col.names, colnames(mut.rdpn))
    idx.col.exp = match(int.col.names, colnames(exp.rdpn))
    
    mut.rdpn = mut.rdpn[,idx.col.mut] #11012 973
    exp.rdpn = exp.rdpn[,idx.col.exp] #11012 973
    
    if (length(colnames(mut.rdpn)) == length(colnames(exp.rdpn))) {
      print("it matched")
    } else {
      stop("it did not match")
      
    }
    
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
  # colnames(raw.prop.score) = substr(pat.ids.test,1,16)
  # 
  # colnames(mut.mat) = substr(colnames(mut.mat),1,16)
  
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
  
  # idx.mut.id = match(colnames(pval.mat), colnames(mut.mat)) 
  # mut.mat = mut.mat[,idx.mut.id] 
  
  mut.mat = mut.mat[,intersect(colnames(pval.mat), colnames(mut.mat))]
  
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
  for (i in 1:random.num) {  # when there are 50 iterations 
    print(i)
    if (length(ck.rdpn.file) == 150) {
      next
    }
    print("here")
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
    # colnames(mut.rdpn) = substr(colnames(mut.rdpn), 1, 16)
    
    #exp.col = do.call("rbind", strsplit(colnames(exp.rdpn),"-")) 
    #n.exp.col = paste(exp.col[,1],exp.col[,2],exp.col[,3],sep="-")
    #colnames(exp.rdpn) = n.exp.col 
    # colnames(exp.rdpn) = substr(colnames(exp.rdpn), 1, 16)
    dim(mut.rdpn)
    dim(exp.rdpn)
    
    int.col.names = intersect(colnames(mut.rdpn), colnames(exp.rdpn)) #151
    
    idx.col.mut = match(int.col.names, colnames(mut.rdpn))
    idx.col.exp = match(int.col.names, colnames(exp.rdpn))
    
    mut.rdpn = mut.rdpn[,idx.col.mut] #11012 973
    exp.rdpn = exp.rdpn[,idx.col.exp] #11012 973
    
    if (length(colnames(mut.rdpn)) == length(colnames(exp.rdpn))) {
      print("it matched")
    } else {
      stop("it did not match")
      
    }
    
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
  
  # pat.ids.test = colnames(raw.prop.score)
  # colnames(raw.prop.score) = substr(pat.ids.test,1,16)
  # 
  # colnames(mut.mat) = substr(colnames(mut.mat),1,16)
  
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
    norm.logi2 = ifelse(is.na(norm.logi) , FALSE , norm.logi) 
    rand.mat[idx, idx.ID] = rand.mat[idx, idx.ID] + norm.logi2 
  }
  toc() 
  
  pval.mat = rand.mat / count.mat 
  colnames(pval.mat) = colnames(raw.prop.score) 
  
  # idx.mut.id = match(colnames(pval.mat), colnames(mut.mat))
  # mut.mat = mut.mat[,idx.mut.id] 
  
  mut.mat = mut.mat[,intersect(colnames(pval.mat), colnames(mut.mat))]
  
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


##########################################################################################################################
## 04. func_GC_hypergeometirc_54_KEGG
##########################################################################################################################

#loading R packages
library(data.table)
library(igraph)
library(data.table)
library(survival)
library(survminer)

# @@@@@@@@@@@@@@@@@
# global setting
# @@@@@@@@@@@@@@@@@

filepath = "/home/seokwon/"
ref_path = paste0(filepath, "99.reference/")
num.edge.multi = 10 #or 100 (for now)
random.num = 150
sourcepath = paste0(ref_path, "RDPN/",c(num.edge.multi),"X" ) 

# @@@@@@@@@@@@@@@@@
# 1,2 setting
# @@@@@@@@@@@@@@@@@

ppi.date = "20220117" 
path.dat = "KEGG" #"Iorio2018" #"DKShin"
path.file = paste0(ref_path ,"Kegg_pathway_genes.rds") 
#loading PPI Backbone
backbone.file = paste0(ref_path,"ppi_backbone_" ,ppi.date,".rds")
g.ppi.conn = readRDS(backbone.file)
#loading CHG data 
path = readRDS(path.file)

# @@@@@@@@@@@@@@@@@
# 2 setting
# @@@@@@@@@@@@@@@@@
level = "_net" # ""
filt.genes.n = 5 # at pathway-link level 

# @@@@@@@@@@@@@@@@@
# 3 setting
# @@@@@@@@@@@@@@@@@

# extracting pathway links with having more than 5 genes in a pathway link 
dat = readRDS(paste0(ref_path,"KEGG_pathway_links_intersect_genes.rds"))

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# This script is for hypergeometric test (both pathway and pathway link levels) $$$$$$$
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Direct to phenotype (Net.prop & RDPN --> GC --> 54 KEGG pathway links) $$$$$$$$$$$$$$
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Function define $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# 1.function of pathway-level (e.g., P1, P2, P3, ..., P54) $$$$$$$$$$$$$$$$$$$$$$$$$$$$
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

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

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# 2.function of pathway-link level (e.g., P1-P2, P3-P10, ..., )$$$$$$$$$$$$$$$$$$$$$$$$
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

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

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# 3. function of obtaining p values from the survival analysis $$$$$$$$$$$$$$$$$$$$$$$$
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

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
  mut.count$patient_id = rownames(mut.count) 
  mut.filt = mut.count[which(mut.count$counts < 1000), ]
  int.id = intersect(colnames(hyper.filt), mut.filt$patient_id) 
  idx.id = match(int.id, colnames(hyper.filt))
  hyper.filt = hyper.filt[, idx.id ]
  
  #loading clinical information for survival analysis 
  # Patient ID -> submitter_id, Overall_survival -> vital_status , Overall_survival(Month) -> days_to_death or days_to_last_follow_up
  # time -> overall_survival  
  
  cli = fread(paste0(main.path,CancerType_in,"_clinical_data.tsv")) 
  cli
  cli_surv = cli[,
                 c("Patient ID",
                   "Overall Survival Days",
                   "Overall Survival Status")]
  
  # cli_surv$deceased = cli_surv$vital_status == "Dead"
  # cli_surv$overall_survival = ifelse(cli_surv$deceased,
  #                                    cli_surv$days_to_death,
  #                                    cli_surv$days_to_last_follow_up)
  # 
  # cli_surv = cli_surv[!is.na(cli_surv$vital_status),]
  # 
  cli_surv$status = NA
  cli_surv$status[which(cli_surv$`Overall Survival Status` == "1:DECEASED")] = 1
  cli_surv$status[which(cli_surv$`Overall Survival Status` == "0:LIVING")] = 0
  
  # 
  # cli_surv$overall_survival = cli_surv$`Overall Survival (Months)`*30
  cli_surv$overall_survival = cli_surv$`Overall Survival Days`
  cli_surv = cli_surv[!is.na(cli_surv$overall_survival),]
  cli_surv$status = as.numeric(cli_surv$status)
  # cli_surv
  
  wrong_cluster = vector()
  for (p in 1:as.numeric(dim(hyper.filt)[1])) {
    hyper.each = hyper.filt[p,]
    hyper.dat = as.data.frame(cbind(id = substr(names(hyper.each),1,16), pvals= as.numeric(hyper.each)))
    #hyper.dat$group = "" 
    hyper.dat['group'] = ""
    hyper.dat$group[which(hyper.dat$pvals <= 0.05)] = "enriched"
    hyper.dat$group[which(hyper.dat$pvals > 0.05)] = "not_enriched"
    int = intersect(hyper.dat$id, cli_surv$`Patient ID`)
    idx = match(int, cli_surv$`Patient ID`) 
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
    #fit = survfit(Surv(cli_surv_filt$overall_survaival,cli_surv_filt$status) ~ cluster, data = cli_surv_filt)
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
  
  rownames(pval.com) = rownames(hyper.filt)[-wrong_cluster]
  colnames(pval.com) = gsub(".*\\.", "", CancerType)
  print(pval.com)
  
  saveRDS(pval.com, paste0(main.path, "/surv_pvals_54_", path.dat, level, "_GC_", ppi.date, "_filt", filt.genes.n, ".rds"))
  
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
  mut.count$patient_id = rownames(mut.count)
  mut.filt = mut.count[which(mut.count$counts < 1000), ]
  int.id = intersect(colnames(hyper.filt), mut.filt$patient_id) 
  idx.id = match(int.id, colnames(hyper.filt))
  hyper.filt = hyper.filt[, idx.id ]
  
  #loading clinical information for survival analysis 
  # Patient ID -> submitter_id, Overall_survival -> vital_status , Overall_survival(Month) -> days_to_death or days_to_last_follow_up
  # time -> overall_survival  
  
  cli = fread(paste0(main.path,CancerType_in,"_clinical_data.tsv")) 
  cli
  cli_surv = cli[,
                 c("Patient ID",
                   "Overall Survival Days",
                   "Overall Survival Status")]
  
  # cli_surv$deceased = cli_surv$vital_status == "Dead"
  # cli_surv$overall_survival = ifelse(cli_surv$deceased,
  #                                    cli_surv$days_to_death,
  #                                    cli_surv$days_to_last_follow_up)
  # 
  # cli_surv = cli_surv[!is.na(cli_surv$vital_status),]
  # 
  cli_surv$status = NA
  cli_surv$status[which(cli_surv$`Overall Survival Status` == "1:DECEASED")] = 1
  cli_surv$status[which(cli_surv$`Overall Survival Status` == "0:LIVING")] = 0
  
  # 
  # cli_surv$overall_survival = cli_surv$`Overall Survival (Months)`*30
  cli_surv$overall_survival = cli_surv$`Overall Survival Days`
  cli_surv = cli_surv[!is.na(cli_surv$overall_survival),]
  cli_surv$status = as.numeric(cli_surv$status)
  
  wrong_cluster = vector()
  for (p in 1:as.numeric(dim(hyper.filt)[1])) {
    hyper.each = hyper.filt[p,]
    hyper.dat = as.data.frame(cbind(id = substr(names(hyper.each),1,16), pvals= as.numeric(hyper.each)))
    #hyper.dat$group = "" 
    hyper.dat['group'] = ""
    hyper.dat$group[which(hyper.dat$pvals <= 0.05)] = "enriched"
    hyper.dat$group[which(hyper.dat$pvals > 0.05)] = "not_enriched"
    int = intersect(hyper.dat$id, cli_surv$`Patient ID`)
    idx = match(int, cli_surv$`Patient ID`) 
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
  # print(wrong_cluster)
  if (is.logical(wrong_cluster) || wrong_cluster == 0) {
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
pval.path = paste0(main.path, "03.RDPN/RDPN_",c(num.edge.multi),"X/")
pval.mat = readRDS(paste0(pval.path,"net_prop_pval_",random.num ,"rep_total_samples_",ppi.date,".rds")) # output of 03.RDPN.r 

# Execute function 1
GC.hyper.pathway.direct(CancerType_in, ppi.date, path.dat, path.file)

# Execute function 2
GC.hyper.path.net.direct(CancerType_in, ppi.date, path.dat, path.file) ## then "phyper_enrichment_54_KEGG_net_GC_20220117.rds" will be saved

# Execute function 3
cal.surv.p(CancerType_in, ppi.date, path.dat, level, filt.genes.n) # this will generate "surv_pvals_54_KEGG_net_GC_20220117_filt5.rds" 

# Execute function 4
cal.surv.each.p(CancerType_in, ppi.date, path.dat )


##########################################################################################################################
## 05. Algorithm 2
##########################################################################################################################

library(survminer) 
library(ggplot2) 
library(pheatmap)
library(survival) 
library(data.table) 
library(RColorBrewer) 
library(NbClust)
library(cluster) 
library(factoextra) 
library(dplyr)
library(stringr)
library(nnet)
require(RSNNS)
require(clusterGeneration)

# 1) prepare for NN
# (1) for TCGA data upload

CancerType = "TCGA-LUAD"
main.path_tc = paste0(filepath, "01.externaldata/luad/", CancerType)

#   A. for each pathway
surv.pvals = readRDS(paste0(main.path_tc,"/surv_pvals_54_KEGG_GC_20220117.rds"))
phyper = readRDS(paste0(main.path_tc, "/phyper_enrichment_54_KEGG_GC_20220117.rds"))

#   B. for pathwaylink
surv.pvals = readRDS(paste0(main.path_tc,"/surv_pvals_54_KEGG_net_GC_20220117_filt5.rds"))
phyper = readRDS(paste0(main.path_tc, "/phyper_enrichment_54_KEGG_net_GC_20220117.rds"))

mut.mat = readRDS(paste0(main.path_tc,"/", CancerType,"_mut_count_filt_data.rds"))
mut.count = as.data.frame(colSums(mut.mat))
colnames(mut.count) = "counts"

# (1) - clinical data upload

cli = fread(paste0(ref_path , "all_clin_indexed.csv"))
cli_surv = cli[cli$project == CancerType,
               c("submitter_id",
                 "vital_status",
                 "days_to_death",
                 "days_to_last_follow_up")]

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


# (1) - filtering only 01A

mut.count$patient_id = substr(rownames(mut.count), 1, 16)
mut.count$forfilt =  rownames(mut.count)
mut.count.filt = mut.count %>% filter(str_detect(forfilt, "-01A-")) 
mut.filt = mut.count.filt[which(mut.count.filt$counts < 1000),]
mut.filt$forfilt = NULL

int.id = intersect(colnames(phyper), mut.filt$patient_id)
idx.id = match(int.id, colnames(phyper))
phyper.filt = phyper[,idx.id]

# (1) - transform -log p.val and transverse

phyper.filt_t = t(phyper.filt)
phyper.filt_t_num= phyper.filt_t[,-which(colSums(phyper.filt_t) == dim(phyper.filt_t)[1])]

phyper.filt_t_log = -log(phyper.filt_t)

phyper.filt_t_log = t(phyper.filt_t_log)
colnames(phyper.filt_t_log) = substr(colnames(phyper.filt_t_log), 1,12)
df.phyper.filt_t_log = as.data.frame(phyper.filt_t_log)
df.phyper.filt_t_log_num = data.frame(matrix(nrow = length(rownames(df.phyper.filt_t_log))))
for (all_sub in cli_surv$submitter_id) {
  if (all_sub %in% colnames(df.phyper.filt_t_log)) {
    # df.phyper.all.link_num = cbind(df.phyper.all.link_num,df.phyper.filt_num[,all_sub])
    tmp.phyer.filt_num = as.data.frame(df.phyper.filt_t_log[,all_sub])
    colnames(tmp.phyer.filt_num) = all_sub
    rownames(tmp.phyer.filt_num) = rownames(df.phyper.filt_t_log)
    # print(tmp.phyer.filt_num)
    df.phyper.filt_t_log_num = cbind(df.phyper.filt_t_log_num,tmp.phyer.filt_num)
    tmp.phyer.filt_num = NULL
  }
  
} 
df.phyper.filt_t_log_num
df.phyper.filt_t_log_num= df.phyper.filt_t_log_num[,-1]

df.phyper.filt_t_log_num = t(df.phyper.filt_t_log_num)

df.phyper.filt_t_log_num = as.data.frame(df.phyper.filt_t_log_num)
df.phyper.filt_t_log_num$submitter_id = rownames(df.phyper.filt_t_log_num)
df.phyper.filt_t_log_num

df.phyper.merge_log_all = as.data.frame(merge(cli_surv , df.phyper.filt_t_log_num, by = "submitter_id"))
rownames(df.phyper.merge_log_all) = df.phyper.merge_log_all$submitter_id

# (1) - Alive or Dead

df.phyper.merge_log_all

df.phyper.merge_log_all3 = df.phyper.merge_log_all[,c(-1,-3,-4,-5,-6,-7)]
df.phyper.merge_log_all3
df.phyper.merge_vital_log_all2= df.phyper.merge_log_all3

colnames(df.phyper.merge_vital_log_all2)= gsub("[[:punct:]]", "", colnames(df.phyper.merge_vital_log_all2))
df.phyper.merge_vital_log_all2= df.phyper.merge_vital_log_all2 %>% relocate("vitalstatus", .after = "P53P54")

# (2) - for external data

main.path
main.path_in = main.path

#   A. for each pathway
surv.pvals_in = readRDS(paste0(main.path_in,"/surv_pvals_54_KEGG_GC_20220117.rds"))
phyper_in = readRDS(paste0(main.path_in, "/phyper_enrichment_54_KEGG_GC_20220117.rds"))

#   B. for pathwaylink
surv.pvals_in = readRDS(paste0(main.path_in,"/surv_pvals_54_KEGG_net_GC_20220117_filt5.rds"))
phyper_in = readRDS(paste0(main.path_in, "/phyper_enrichment_54_KEGG_net_GC_20220117.rds"))

mut.mat_in = readRDS(paste0(main.path, CancerType_in,"_mut_count_filt_data.rds"))
mut.count_in = as.data.frame(colSums(mut.mat_in))
colnames(mut.count_in) = "counts"


# (2) - clinical data upload

cli_sclc_uco = fread(paste0(main.path_in,CancerType_in,"_clinical_data.tsv")) 
cli_sclc_uco
cli_sclc_uco_surv = cli_sclc_uco[,
                                 c("Patient ID",
                                   "Overall Survival (Months)",
                                   "Overall Survival Status")]

cli_sclc_uco_surv$vital_status = NA
cli_sclc_uco_surv$vital_status[which(cli_sclc_uco_surv$`Overall Survival Status` == "1:DECEASED")] = "Dead"
cli_sclc_uco_surv$vital_status[which(cli_sclc_uco_surv$`Overall Survival Status` == "0:LIVING")] = "Alive"

cli_sclc_uco_surv$overall_survival = cli_sclc_uco_surv$`Overall Survival (Months)`*30
cli_sclc_uco_surv = cli_sclc_uco_surv[!is.na(cli_sclc_uco_surv$`Overall Survival Status`),]
cli_sclc_uco_surv$status = as.numeric(cli_sclc_uco_surv$status)


# (2) - filtering only 01A 

mut.count_in$patient_id = rownames(mut.count_in)
mut.filt_in = mut.count_in[which(mut.count_in$counts < 1000),]

int.id_in = intersect(colnames(phyper_in), mut.filt_in$patient_id)
idx.id_in = match(int.id_in, colnames(phyper_in))
phyper.filt_in = phyper_in[,idx.id_in]

# (2) - transform -log p.val and transverse

phyper.filt_in_t = t(phyper.filt_in)
phyper.filt_in_t_num= phyper.filt_in_t[,-which(colSums(phyper.filt_in_t) == dim(phyper.filt_in_t)[1])]

phyper.filt_in_t_log = -log(phyper.filt_in_t)
phyper.filt_in_t_log = t(phyper.filt_in_t_log)

df.phyper.filt_in_t_log = as.data.frame(phyper.filt_in_t_log)
df.phyper.filt_in_t_log_num = data.frame(matrix(nrow = length(rownames(df.phyper.filt_in_t_log))))
for (all_sub in cli_sclc_uco_surv$`Patient ID`) {
  if (all_sub %in% colnames(df.phyper.filt_in_t_log)) {
    # df.phyper.all.link_num = cbind(df.phyper.all.link_num,df.phyper.filt_num[,all_sub])
    tmp.phyer.filt_num = as.data.frame(df.phyper.filt_in_t_log[,all_sub])
    colnames(tmp.phyer.filt_num) = all_sub
    rownames(tmp.phyer.filt_num) = rownames(df.phyper.filt_in_t_log)
    # print(tmp.phyer.filt_num)
    df.phyper.filt_in_t_log_num = cbind(df.phyper.filt_in_t_log_num,tmp.phyer.filt_num)
    tmp.phyer.filt_num = NULL
  }
  
} 
df.phyper.filt_in_t_log_num
df.phyper.filt_in_t_log_num= df.phyper.filt_in_t_log_num[,-1]

df.phyper.filt_in_t_log_num = t(df.phyper.filt_in_t_log_num)

df.phyper.filt_in_t_log_num = as.data.frame(df.phyper.filt_in_t_log_num)
df.phyper.filt_in_t_log_num$submitter_id = rownames(df.phyper.filt_in_t_log_num)
df.phyper.filt_in_t_log_num

cli_sclc_uco_surv$submitter_id = cli_sclc_uco_surv$`Patient ID`

df.phyper.merge_in_log_all = as.data.frame(merge(cli_sclc_uco_surv , df.phyper.filt_in_t_log_num, by = "submitter_id"))
rownames(df.phyper.merge_in_log_all) = df.phyper.merge_in_log_all$submitter_id

# (2) - Alive or Dead

df.phyper.merge_in_log_all_filt = df.phyper.merge_in_log_all

df.phyper.merge_log_in3 = df.phyper.merge_in_log_all_filt[,c(-1,-2,-3,-4,-6,-7)]
df.phyper.merge_log_in3
df.phyper.merge_vital_log_in2= df.phyper.merge_log_in3
test.df = df.phyper.merge_vital_log_in2[,-1]
test.filt.df = test.df[,-which(colSums(test.df) == 0)]
test.filt.df = scale(test.filt.df)
colnames(df.phyper.merge_vital_log_in2)= gsub("[[:punct:]]", "", colnames(df.phyper.merge_vital_log_in2)) # '-' does not allow
df.phyper.merge_vital_log_in2= df.phyper.merge_vital_log_in2 %>% relocate("vitalstatus", .after = "P53P54")
df.phyper.merge_vital_log_in2 # this will be input data for validate
df.phyper.merge_vital_log_in2[,-length(colnames(df.phyper.merge_vital_log_in2))] = scale(df.phyper.merge_vital_log_in2[,-length(colnames(df.phyper.merge_vital_log_in2))])

# 2) build a model

train_val = sample(1:nrow(df.phyper.merge_vital_log_all2), nrow(df.phyper.merge_vital_log_all2)*0.7) 
test_val = -train_val

# (1) train by TCGA data
library(neuralnet)
first_layer = 20
second_layer = 20
third_layer = 20
# hidden = c(first_layer,second_layer )
hidden = c(first_layer,second_layer, third_layer )

cli_lu_log_each_neural <- neuralnet( vitalstatus~., data=df.phyper.merge_vital_log_all2[train_val,], hidden=hidden, err.fct = "sse" ,algorithm = "rprop+",linear.output=TRUE, threshold = 0.08,stepmax=1e6)

# plot(cli_lu_log_each_neural)

#   A. TCGA train set
cli_lu_log_train_neural.results <- neuralnet::compute(cli_lu_log_each_neural, df.phyper.merge_vital_log_all2[train_val,])
predicted_train_result = cli_lu_log_train_neural.results$net.result

df.predicted_train_result = as.data.frame(predicted_train_result)
df.predicted_train_result$prediction = NA
for (predic in 1:length(rownames(df.predicted_train_result))) {
  if (df.predicted_train_result[predic,]$V1 > df.predicted_train_result[predic,]$V2) {
    df.predicted_train_result[predic,]$prediction = "Alive"
  } else {
    df.predicted_train_result[predic,]$prediction = "Dead"
  }
}


results.cli_train_alive_log_neural <- data.frame(actual = df.phyper.merge_vital_log_all2[train_val,]$vitalstatus, prediction = df.predicted_train_result$prediction)

# (1) - prediction
predicted_alive_train_neural_log = table(prediction = results.cli_train_alive_log_neural$prediction, actual = results.cli_train_alive_log_neural$actual)

acc_alive_train_neural = round((predicted_alive_train_neural_log[1,1] + predicted_alive_train_neural_log[2,2]) / sum(predicted_alive_train_neural_log) * 100, 2)
acc_alive_train_neural

#   B. TCGA test set
cli_lu_log_test_neural.results <- neuralnet::compute(cli_lu_log_each_neural, df.phyper.merge_vital_log_all2[test_val,])
predicted_test_result = cli_lu_log_test_neural.results$net.result

df.predicted_test_result = as.data.frame(predicted_test_result)
df.predicted_test_result$prediction = NA
for (predic in 1:length(rownames(df.predicted_test_result))) {
  if (df.predicted_test_result[predic,]$V1 > df.predicted_test_result[predic,]$V2) {
    df.predicted_test_result[predic,]$prediction = "Alive"
  } else {
    df.predicted_test_result[predic,]$prediction = "Dead"
  }
}

results.cli_test_alive_log_neural <- data.frame(actual = df.phyper.merge_vital_log_all2[test_val,]$vitalstatus, prediction = df.predicted_test_result$prediction)

# (1) - prediction
predicted_alive_test_neural_log = table(prediction = results.cli_test_alive_log_neural$prediction, actual = results.cli_test_alive_log_neural$actual)

acc_alive_test_neural = round((predicted_alive_test_neural_log[1,1] + predicted_alive_test_neural_log[2,2]) / sum(predicted_alive_test_neural_log) * 100, 2)
acc_alive_test_neural

# 3) validate with external data set
cli_ex_log_in_neural.results <- neuralnet::compute(cli_lu_log_each_neural, df.phyper.merge_vital_log_in2)

predicted_in_result = cli_ex_log_in_neural.results$net.result

df.predicted_in_result = as.data.frame(predicted_in_result)
df.predicted_in_result$prediction = NA
for (predic in 1:length(rownames(df.predicted_in_result))) {
  if (df.predicted_in_result[predic,]$V1 > df.predicted_in_result[predic,]$V2) {
    df.predicted_in_result[predic,]$prediction = "Alive"
  } else {
    df.predicted_in_result[predic,]$prediction = "Dead"
  }
}

results.cli_ex_alive_log_neural <- data.frame(actual = df.phyper.merge_vital_log_in2$vitalstatus, prediction = df.predicted_in_result$prediction)

# (1) - prediction
predicted_alive_ex_neural_log = table(prediction = results.cli_ex_alive_log_neural$prediction, actual = results.cli_ex_alive_log_neural$actual)

acc_alive_ex_neural = round((predicted_alive_ex_neural_log[1,1] + predicted_alive_ex_neural_log[2,2]) / sum(predicted_alive_ex_neural_log) * 100, 2)
acc_alive_ex_neural

neuralnet::gwplot(cli_lu_log_each_neural, selected.covariate = 3,selected.response = 2)

# 4) model validation with ROC curve
library(ROCR)
results.cli_train_alive_log_neural$actual_bi = NA
results.cli_train_alive_log_neural$predic_bi = NA
results.cli_train_alive_log_neural$actual_bi[which(results.cli_train_alive_log_neural$actual == "Alive")] = 1
results.cli_train_alive_log_neural$actual_bi[which(results.cli_train_alive_log_neural$actual == "Dead")] = 0
results.cli_train_alive_log_neural$predic_bi[which(results.cli_train_alive_log_neural$prediction == "Alive")] = 1
results.cli_train_alive_log_neural$predic_bi[which(results.cli_train_alive_log_neural$prediction == "Dead")] = 0

pr_neural = ROCR::prediction( results.cli_train_alive_log_neural$actual_bi, results.cli_train_alive_log_neural$predic_bi)
prf_neural <- performance(pr_neural, measure = "tpr", x.measure = "fpr")
plot(prf_neural, main = "ROC Curve")

# 5) shrink variables

# relative importance function
# relative importance of input variables for Alive or Dead. you can change the str

source(paste0(ref_path, "gar_fun_neuralnet.R"))

rel.imp_neural = gar.fun.neural('Alive',cli_lu_log_each_neural , bar.plot=FALSE)
# rel.imp_neural = gar.fun.neural('Dead',cli_lu_log_each_neural , bar.plot=FALSE)

x.names = cli_lu_log_each_neural$model.list$variables
abs.rel.imp_neural = rel.imp_neural^2

to_plo <- abs.rel.imp_neural[order(abs.rel.imp_neural$rel.imp),,drop = F]

to_plo$x.names = factor(x.names[order(abs.rel.imp_neural)], levels = x.names[order(abs.rel.imp_neural)])

ggplot(to_plo, aes(x = x.names, y = rel.imp, fill = rel.imp,
                   colour = rel.imp)) + 
  geom_bar(stat = 'identity') + 
  scale_x_discrete(element_blank()) +
  scale_y_continuous(c("Dead"))
tail(to_plo)
best_rel_vari = to_plo[which(to_plo$rel.imp>0.25),]

# 6) test for best
variables = rownames(best_rel_vari)

new_vari = c()
for (best in variables) {
  new_vari = c(new_vari, best)
  new_vari = c(new_vari, "+")
}
new_vari = new_vari[-length(new_vari)]

cat(new_vari)

# (1) train by TCGA data
cli_lu_log_best_each_neural <- neuralnet( vitalstatus~ P40P42 + P38P44 + P38P42 + P40P44 , data=df.phyper.merge_vital_log_all2[train_val,], hidden=hidden, err.fct = "sse" ,algorithm = "rprop+",linear.output=TRUE, threshold = 0.08,stepmax=1e6)

# (2-1) validation with TCGA train set
cli_lu_log_best_train_neural.results <- neuralnet::compute(cli_lu_log_best_each_neural, df.phyper.merge_vital_log_all2[train_val,])
predicted_best_train_result = cli_lu_log_best_train_neural.results$net.result

df.predicted_best_train_result = as.data.frame(predicted_best_train_result)
df.predicted_best_train_result$prediction = NA
for (predic in 1:length(rownames(df.predicted_best_train_result))) {
  if (df.predicted_best_train_result[predic,]$V1 > df.predicted_best_train_result[predic,]$V2) {
    df.predicted_best_train_result[predic,]$prediction = "Alive"
  } else {
    df.predicted_best_train_result[predic,]$prediction = "Dead"
  }
}

results.cli_train_best_alive_log_neural <- data.frame(actual = df.phyper.merge_vital_log_all2[train_val,]$vitalstatus, prediction = df.predicted_best_train_result$prediction)

#   A. acc results

predicted_best_alive_train_neural_log = table(prediction = results.cli_train_best_alive_log_neural$prediction, actual = results.cli_train_best_alive_log_neural$actual)

acc_best_alive_train_neural = round((predicted_best_alive_train_neural_log[1,1] + predicted_best_alive_train_neural_log[2,2]) / sum(predicted_best_alive_train_neural_log) * 100, 2)
acc_best_alive_train_neural

# (2-2) validation with TCGA test set
cli_lu_log_best_test_neural.results <- neuralnet::compute(cli_lu_log_best_each_neural, df.phyper.merge_vital_log_all2[test_val,])
predicted_best_test_result = cli_lu_log_best_test_neural.results$net.result

df.predicted_best_test_result = as.data.frame(predicted_best_test_result)
df.predicted_best_test_result$prediction = NA
for (predic in 1:length(rownames(df.predicted_best_test_result))) {
  if (df.predicted_best_test_result[predic,]$V1 > df.predicted_best_test_result[predic,]$V2) {
    df.predicted_best_test_result[predic,]$prediction = "Alive"
  } else {
    df.predicted_best_test_result[predic,]$prediction = "Dead"
  }
}

results.cli_test_best_alive_log_neural <- data.frame(actual = df.phyper.merge_vital_log_all2[test_val,]$vitalstatus, prediction = df.predicted_best_test_result$prediction)

#   A. acc results

predicted_best_alive_test_neural_log = table(prediction = results.cli_test_best_alive_log_neural$prediction, actual = results.cli_test_best_alive_log_neural$actual)

acc_best_alive_test_neural = round((predicted_best_alive_test_neural_log[1,1] + predicted_best_alive_test_neural_log[2,2]) / sum(predicted_best_alive_test_neural_log) * 100, 2)
acc_best_alive_test_neural

# (3) validation with external data by best variables
cli_lu_log_best_in_neural.results <- neuralnet::compute(cli_lu_log_best_each_neural, df.phyper.merge_vital_log_in2)
predicted_best_in_result = cli_lu_log_best_in_neural.results$net.result

df.predicted_best_in_result = as.data.frame(predicted_best_in_result)
df.predicted_best_in_result$prediction = NA
for (predic in 1:length(rownames(df.predicted_best_in_result))) {
  if (df.predicted_best_in_result[predic,]$V1 > df.predicted_best_in_result[predic,]$V2) {
    df.predicted_best_in_result[predic,]$prediction = "Alive"
  } else {
    df.predicted_best_in_result[predic,]$prediction = "Dead"
  }
}


results.cli_ex_best_alive_in_log_neural <- data.frame(actual = df.phyper.merge_vital_log_in2$vitalstatus, prediction = df.predicted_best_in_result$prediction)

#   A. acc results

predicted_best_alive_ex_neural_log = table(prediction = results.cli_ex_best_alive_in_log_neural$prediction, actual = results.cli_ex_best_alive_in_log_neural$actual)

acc_best_alive_ex_neural = round((predicted_best_alive_ex_neural_log[1,1] + predicted_best_alive_ex_neural_log[2,2]) / sum(predicted_best_alive_ex_neural_log) * 100, 2)
acc_best_alive_ex_neural


# 7) save output

hidden_char = ""

for (layer in 1:length(hidden)) {
  if (layer == 1) {
    hidden_char= paste0(hidden[layer])
    
  } else {
    hidden_char= paste0(hidden_char, ",", hidden[layer])
  }
  
}

variables_char = ""
variables_name = ""
for (b_vari in 1:length(variables)) {
  if (b_vari == 1) {
    variables_char= paste0(variables[b_vari])
    variables_name = paste0(variables[b_vari])
    
  } else {
    variables_char= paste0(variables_char, ",", variables[b_vari])
    variables_name = paste0(variables_name, "_", variables[b_vari])
  }
  
}

merge_results = data.frame(acc_train = acc_alive_train_neural, acc_test = acc_alive_test_neural, acc_ex = acc_alive_ex_neural,
                           acc_best_train = acc_best_alive_train_neural , acc_best_test = acc_best_alive_test_neural , 
                           acc_best_ex = acc_best_alive_ex_neural,
                           hidden = hidden_char, best_vari = variables_char)

write.csv(merge_results, paste0(main.path ,"results_of_",variables_name,"neuralnet_pathwaylink.csv"), row.names = CancerType_in)


####### fin









# rel.imp_neural = gar.fun.neural('Alive',cli_log_each_neural , bar.plot=FALSE)$rel.imp
rel.imp_best_neural = gar.fun.neural('Alive',cli_lu_log_best_each_neural , bar.plot=FALSE)
rel.imp_best_neural = gar.fun.neural('Dead',cli_lu_log_best_each_neural , bar.plot=FALSE)

x.best.names = cli_lu_log_best_each_neural$model.list$variables
abs.rel.imp_neural = rel.imp_neural^2

tobest_plo <- rel.imp_best_neural[order(rel.imp_best_neural$rel.imp),,drop = F]

tobest_plo$x.best.names = factor(x.best.names[order(rel.imp_best_neural)], levels = x.best.names[order(rel.imp_best_neural)])

ggplot(tobest_plo, aes(x = x.best.names, y = rel.imp, fill = rel.imp,
                       colour = rel.imp)) + 
  geom_bar(stat = 'identity') + 
  scale_x_discrete(element_blank()) +
  scale_y_continuous(c("Dead"))



##### plot 

plot.nnet(cli_log_neural,pos.col='darkgreen',neg.col='darkblue',alpha.val=0.7,rel.rsc=15,
          circle.cex=10,cex=1.4,
          circle.col='brown')

# find weight
plot.nnet(cli_nnv, wts.only=T)

# plot.nnet(cli_nn,pos.col='darkgreen',neg.col='darkblue',alpha.val=0.7,rel.rsc=15,
#           circle.cex=10,cex=1.4,circle.col='brown',all.in='Sepal W.',all.out='Alive')

plot.nnet(cli_nnv,pos.col='darkgreen',neg.col='darkblue',alpha.val=0.7,rel.rsc=15,
          circle.cex=10,cex=1.4,circle.col='brown',all.out='Alive')


##
source('https://gist.githubusercontent.com/fawda123/6860630/raw/b8bf4a6c88d6b392b1bfa6ef24759ae98f31877c/lek_fun.r')

library(NeuralNetTools)
garson(cli_log_each_nnv) # it will not work at two output mode
lekprofile(cli_log_each_nnv)



#color vector based on relative importance of input values
# the value of 22 meaning is the number of input column
cols<-colorRampPalette(c('green','red'))(length(rank(rel.imp)))[rank(rel.imp)]

# for nnet
toprank = as.data.frame(tail(sort(rank(rel.imp),  decreasing = TRUE), n=10))
colnames(toprank) = "bottom10"
toprank$bottom10 = seq(length(rel.imp)-9,length(rel.imp),1)
toprank2 = as.data.frame(head(sort(rank(rel.imp),  decreasing = TRUE), n=10))
colnames(toprank2) = "top10"
toprank2$top10 = seq(1,10,1)
write.csv(toprank, "weight_bottom10_each_alllink.csv")
write.csv(toprank2, "weight_top10_each_alllink.csv")
write.csv(rel.imp , "weight_for_each_neuralnet.csv")

# for neuralnet
out$variables = rownames(out)
out.order = out[c(order(out$rel.imp)),]
toprank = out.order[1:10,]
toprank$variables = NULL
colnames(toprank) = "bottom10"
toprank$weight = seq(length(rownames(out)),length(rownames(out))-9,-1)

toprank2 = out.order[c(length(rownames(out.order))-9):length(rownames(out.order)),]
toprank2$variables = NULL
colnames(toprank2) = "top10"
toprank2$weight = seq(10,1,-1)
write.csv(out , "weight_for_each_neuralnet.csv")
write.csv(toprank, "weight_bottom10_each_neuralnet.csv")
write.csv(toprank2, "weight_top10_each_neuralnet.csv")




# weight 뽑기 
neuralweights(cli_log_each_neural)
##
#plotting function
source_url('https://gist.githubusercontent.com/fawda123/7471137/raw/466c1474d0a505ff044412703516c34f1a4684a5/nnet_plot_update.r')

# plot model with new color vector
# separate colors for input vectors using a list for 'circle.col'
# all.out does not automatically matched with the results of rel.imp. 
# So you can manually change the all.out value correlate with rel.imp

plot(cli_log_each_neural,circle.col=list(cols,'lightblue'), all.out = 1)

plot.nnet(cli_log_each_neural,pos.col='darkgreen',neg.col='darkblue',alpha.val=0.7,rel.rsc=15,
          circle.cex=5,cex=1.4,circle.col=list(cols,'lightblue'))

