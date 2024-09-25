RDPN = function(expression_data, 
                mutation_data, 
                net_propa_data, 
                random_out_path = "RDPN/", 
                final_out_path = "RDPN_res/",
                num.edge.multi = 10, 
                random.num = 150){
  
  
  #loading tools 
  require(igraph) 
  require(data.table)
  require(foreach) 
  require(doParallel) 
  require(tictoc)
  
  ######################################################################################################################
  ## Network propagation using 'ppi_backbone_*_RDPN_#.rds' 
  ######################################################################################################################
  
  # input data (exp data and mut data)
  mut.mat = mutation_data
  exp.log.mat = expression_data
  raw.prop.score = net_propa_data
  query.genes = rownames(raw.prop.score) #12902
  
  # savepath == random_out_path
  # result_savepath == final_out_path
  
  # create folder
  if(!dir.exists(random_out_path)){
    dir.create(random_out_path)
    print(paste0("Created folder: ", random_out_path))
  } else {
    print(paste0("Folder already exists: ", random_out_path))
  }
  
  if(!dir.exists(final_out_path)){
    dir.create(final_out_path)
    print(paste0("Created folder: ", final_out_path))
  } else {
    print(paste0("Folder already exists: ", final_out_path))
  }
  
  for (i in 1:random.num) { 
    
    g.rdpn = readRDS(paste0("reference/RDPN/10X/ppi_backbone_20220117_RDPN",c(i),".rds")) 
    g.rdpn.sub = induced_subgraph(g.rdpn, query.genes) 
    c.rdpn = components(g.rdpn.sub) 
    
    query.names = names(which(c.rdpn$membership == which(c.rdpn$csize == max(c.rdpn$csize)))) 
    g.rdpn.gs = induced_subgraph(g.rdpn, query.names) 
    
    #fixing exp.mat 
    idx.mat = match(query.names, rownames(exp.log.mat)) 
    exp.rdpn = exp.log.mat[idx.mat, ] 
    
    #fixing mut.mat 
    common.mut = intersect(query.names, rownames(mut.mat)) 
    idx.mut = match(common.mut, rownames(mut.mat))
    mut.rdpn = as.matrix(mut.mat[idx.mut,]) 
    
    # colnames(mut.rdpn) = substr(colnames(mut.rdpn), 1, 16)
    # colnames(exp.rdpn) = substr(colnames(exp.rdpn), 1, 16)
    
    int.col.names = intersect(colnames(mut.rdpn), colnames(exp.rdpn))
    
    if (length(int.col.names)==length(colnames(exp.rdpn))){
      idx.col.mut = match(int.col.names, colnames(mut.rdpn)) 
    } else {
      print("IDs are not exactly matched.....please check ! ")
    }
    mut.rdpn = mut.rdpn[,idx.col.mut] 
    
    idx.col.exp = match(int.col.names, colnames(exp.rdpn))
    exp.rdpn = exp.rdpn[,idx.col.exp]
    
    #stretching 
    idx.stretch = match(rownames(mut.rdpn), rownames(exp.rdpn)) 
    n.mut.stretch = matrix(data = 0, nrow = nrow(exp.rdpn), ncol(exp.rdpn))
    colnames(n.mut.stretch) = colnames(exp.rdpn) 
    rownames(n.mut.stretch) = rownames(exp.rdpn) 
    
    #binarizing and stretching the mutation matrix 
    for (j in 1:ncol(mut.rdpn)) {
      query = mut.rdpn[,j]
      query[query > 0] = 1 
      n.mut.stretch[idx.stretch,j] = query 
    }
    
    idx.empty = colSums(n.mut.stretch) > 1 
    n.mut.stretch = n.mut.stretch[,idx.empty] 
    
    exp.rdpn = exp.rdpn[,idx.empty]
    sn = 1e-14 
    exp.rdpn[exp.rdpn == 0] = sn 
    
    #testing stuff for compressed matrix 
    adj.rdpn = as.matrix(as_adjacency_matrix(g.rdpn.gs))
    
    N.pat <- dim(exp.rdpn)[2] 
    N.net <- dim(adj.rdpn)[1] 
    alphav <- 0.7 
    
    # for doing parallel processing
    cl <- makeCluster(20) 
    registerDoParallel(cl) 
    
    tic() 
    # .combine = > combine mode of results (it must be r function)
    prop.res.all = foreach(j=1:N.pat, .combine =cbind) %dopar% {
      source("NetGPT/func_netprop.r")
      prop.res.new <- net.propagation(exp.rdpn[,j], adj.rdpn, n.mut.stretch[,j], alphav) 
      prop.res.new       
    }
    toc() 
    stopCluster(cl) 
    
    rownames(prop.res.all) = rownames(exp.rdpn) 
    colnames(prop.res.all) = colnames(exp.rdpn) 
    
    saveRDS(prop.res.all, file = paste0(random_out_path,"net_prop_RDPN",c(i),".rds"))
    
  }
  
  
  ######################################################################################################################
  ## calculating the p-values by comparing the original and generated with the random degree preserving network 
  ######################################################################################################################
  
  count.mat = matrix(data = 0, nrow = nrow(net_propa_data), ncol = ncol(net_propa_data))
  rownames(count.mat) = rownames(net_propa_data) 
  
  rand.mat = matrix(data = 0, nrow = nrow(net_propa_data), ncol = ncol(net_propa_data)) 
  rownames(rand.mat) = rownames(net_propa_data) 
  
  tic() 
  
  for (i in 1:random.num) {
    # print (i)
    rdpn.file = paste0(random_out_path, "net_prop_RDPN", c(i), "_.rds")
    rand.prop.score = readRDS(rdpn.file)
    
    idx = match(rownames(rand.prop.score), rownames(net_propa_data)) 
    idx.ID = match(colnames(rand.prop.score), colnames(net_propa_data)) 
    
    count.mat[idx, idx.ID] = count.mat[idx, idx.ID]+1  #count up when gene is present 
    norm.score = rand.prop.score / net_propa_data[idx, idx.ID]
    norm.logi = norm.score >= 1
    rand.mat[idx, idx.ID] = rand.mat[idx, idx.ID] + norm.logi 
  }
  toc() 
  
  pval.mat = rand.mat / count.mat 
  colnames(pval.mat) = colnames(net_propa_data) 
  
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
  
  saveRDS(pval.mat, file = paste0(final_out_path, "net_prop_pval_", random.num, "rep_total_samples.rds"))
  return(pval.mat)
  
}

  




  





