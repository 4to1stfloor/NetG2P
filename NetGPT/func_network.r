#########################################################################################
### Define network propagation function 
#########################################################################################

# exported from Dr.Shin source file. Keeping it here to do the sourcing easier (be sourced for one of the GPMapper steps)
network_propagation <- function (exp.vec, adj.mat, source.mat, alphav) {
    # Calculating network propagation for given perturbation sources. 

    #input: 
    # exp.vec: expression vector of a subject. log normalization recommended. 
    # adj.mat: adjacency matrix of a nominal PPI network. Symmetric matrix. 0 or 1 
    # source.mat: a colum vector with perturbed node to be 1 for a trial, a matrix for multiple trials. 

    # exp.vec = n.exp.mat[,ipat]
    # adj.mat = n.adj.mat 
    # source.mat = n.mut.abs[,ipat]

    #output: 
    # prop.vec: N x 1 

    library(pcg) 
    library(igraph)
    N.net = dim(adj.mat)[1] 

    # Making stochastic matrix with node links to have the weight of the product of end node's expressions 
    stoc.mat <- sapply(1:N.net, function(x) exp.vec)
    stoc.mat <- stoc.mat*t(stoc.mat)*adj.mat 

    # Normalized stochastic matrix (D^(-1/2) %*% W %*% D^(-1/2)) 
    # which can be modified to [D^(-1/2) * W^(1/2)]*[D^(-1/2) * W^(1/2)]^Transpose for simuation efficiency 
    norm.exp <- 1/sqrt(rowSums(stoc.mat))
    temp.stoc <- norm.exp*sqrt(stoc.mat)
    norm.stoc.mat <- temp.stoc * t(temp.stoc)

    # Calculate propagation for a patient 
    # Solve A.x=b 
    A.oper <- diag(N.net) - alphav * norm.stoc.mat 
    b.source <- (1-alphav)*source.mat 
    x.prop <- pcg(A.oper, b.source) 
    # it takes 4.84s for N.net = 10968 in MAC 2Gz

    return(prop.mat = x.prop)
}
