##########################################################################################################################
## network propagation 
##########################################################################################################################

Network_propagation = function(expression_data, mutation_data, binlimit = 192) {
  
  packages <- c("igraph","tictoc","pcg","foreach","doParallel") 
  
  install_if_missing <- function(pkg) {
    if (!require(pkg, character.only = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
      library(pkg, character.only = TRUE)
    }
  }
  
  sapply(packages, install_if_missing)
  
  library(igraph)
  library(tictoc) 
  library(pcg)
  
  # input data (exp data and mut data)
  mut.mat = mutation_data
  exp.log.mat = expression_data
  
  # protein to protein 
  g.ppi.conn = readRDS(file = "../data/reference/ppi_backbone_20220117.rds")
  
  g.ppi.conn.simp = simplify(g.ppi.conn) 
  gene.backbone = names(V(g.ppi.conn.simp))
  
  gene.query = rownames(exp.log.mat) 
  common_ppi_exp = intersect(gene.backbone, gene.query)
  
  g.ppi.spec = induced.subgraph(g.ppi.conn, common_ppi_exp) 
  c.ppi.spec = components(g.ppi.spec)
  idx.max = which(c.ppi.spec$csize == max(c.ppi.spec$csize)) ## clustered network groups, and then select the biggest network based on the cluster size 
  
  g.ppi.gc = induced.subgraph(g.ppi.spec, names(which(c.ppi.spec$membership == idx.max)))
  n.adj.mat = as.matrix(as_adjacency_matrix(g.ppi.gc)) #as_adjacency_matrix function in igraph package, to create the adjacency matrix for undirected graphs 
  gene.adj = rownames(n.adj.mat) #less than common. because we are isolating the largest component
  
  common.exp = intersect(gene.adj, rownames(exp.log.mat)) 
  idx.exp = match(common.exp, rownames(exp.log.mat))
  n.log.mat = exp.log.mat[idx.exp,] 
  
  sn =1e-14
  n.log.mat[n.log.mat == 0] = sn 
  
  # first shirinking the mut.mat to contain only the common genes found in n.log.mat 
  common = intersect(rownames(mut.mat),rownames(n.log.mat)) 
  idx.shrink = match(common, rownames(mut.mat)) 
  mut.mat.shrink = mut.mat[idx.shrink, ] 
  
  n.exp.mat = n.log.mat 
  
  ###################################################################################################
  # then streching it to the size of n.exp.mat (making mut.mat to have same genes with the expression mat) 
  idx.stretch = match(rownames(mut.mat.shrink), rownames(n.exp.mat))
  
  n.mut.stretch = matrix (data = 0 , nrow=nrow(n.exp.mat), ncol = ncol(n.exp.mat)) 
  colnames(n.mut.stretch) = colnames(n.exp.mat) 
  rownames(n.mut.stretch) = rownames(n.exp.mat) 
  
  for ( i in 1:ncol(n.mut.stretch)) {
    query = mut.mat.shrink[,i]
    query[which(query > 0)] = 1 
    n.mut.stretch[idx.stretch,i]=query 
  }
  ##################################################################################################
  
  idx.muts = which(colSums(n.mut.stretch) > 1) 
  n.mut.stretch = n.mut.stretch[,idx.muts]
  n.exp.mat = n.exp.mat[,idx.muts] 
  
  if (all.equal(rownames(n.exp.mat), rownames(n.mut.stretch)) && all.equal(rownames(n.exp.mat), rownames(n.adj.mat)) == TRUE) {
    ###################################################################################################
    ## Network propagation using parallel, depends on the sample size  
    ###################################################################################################
    
    N.pat <- dim(n.exp.mat)[2] #971 
    N.net <- dim(n.adj.mat)[1] #12902 
    alphav <- 0.7 # propagation parameter (0.3 is assumed as returning back to the starting point, other 0.7 is assumed as going to the end network, not returning)
    
    library(foreach)
    library(doParallel)
    
    # parallel::detectCores()-1 #detect the number of available cores 
    # binlimit = 24 * 8. using 24 cores is 1.5 times faster than 36 cores due to memory issues (= 24 cores for 8 samples each time) 

    if (N.pat < binlimit) {
      cl <- makeCluster(24) 
      registerDoParallel(cl) 
      tic() 
      prop.res.all = foreach (i=1:N.pat, .combine = cbind) %dopar% {
        source("./func_network.r")
        
        prop.res.new <- network_propagation(n.exp.mat[,i], n.adj.mat, n.mut.stretch[,i], alphav) 
        # prop.res.new 
      }
      toc() 
      stopCluster(cl) 
    } else { 
      # subsetting data
      bin = ceiling(N.pat / binlimit) 
      message(paste("data is too large. Dividing into", c(bin), "bins"))
      cl <- makeCluster(24) 
      registerDoParallel(cl) 
      tic() 
      prop.res.all = foreach(i=1:binlimit, .combine = cbind) %dopar% { 
        source("./func_network.r")
        
        prop.res.new <- network_propagation(n.exp.mat[,i], n.adj.mat, n.mut.stretch[,i], alphav) 
        prop.res.new
      }
      toc() 
      stopCluster(cl) 
      for (j in 2:bin) {
        max.query = binlimit*j 
        if (max.query > N.pat) { 
          max.query = N.pat 
        }
        cl <- makeCluster(24)
        registerDoParallel(cl) 
        tic() 
        prop.res.add = foreach(i=((binlimit*(j-1))+1): max.query, .combine = cbind) %dopar% {
          source("./func_network.r") 
          
          prop.res.new <- network_propagation(n.exp.mat[,i], n.adj.mat, n.mut.stretch[,i], alphav) 
          # prop.res.new 
        }
        toc() 
        stopCluster(cl) 
        prop.res.all = cbind(prop.res.all, prop.res.add) 
      }
    }
    
    rownames(prop.res.all) = rownames(n.exp.mat)
    colnames(prop.res.all) = colnames(n.exp.mat) 
    
    return(prop.res.all)
    
    ######################################################################################################################
    
  } 
  
}














