
rm(list = ls())

data.folder <- "/mnt/gluster_server/data/network"
setwd(data.folder)

library(igraph)
library(ggplot2)

#################################################################
######################## Functions ############################## 
score_crosstalk_types <- function(links.data, importance.cut){
  
  library(igraph)
  library(bc3net)
  sub.links.data <- dplyr::filter(links.data, min_max >= importance.cut)
  sub.links.data$p.or.l <- sapply(1:dim(sub.links.data)[1], function(x) ifelse(sub.links.data$from[x] == sub.links.data$to[x], "path", "link"))
  temp.link <- dplyr::filter(sub.links.data, p.or.l == "link")
  temp.pathway <- dplyr::filter(sub.links.data, p.or.l == "path")
  
  g1 <- graph_from_data_frame(sub.links.data, directed = F)
  sim.g1<-simplify(g1, remove.loop=T)
  
  gc.g1 <- getgcc(sim.g1)
  
  glink <- graph_from_data_frame(temp.link, directed = F)
  sim.glink <- simplify(glink, remove.loop=T)
  
  score.between <- max(degree(sim.g1))
  score.crosstalk <- diameter(sim.g1)/length(E(sim.g1))
  score.clique <- transitivity(sim.g1)
  score.close <- mean(closeness(sim.g1, V(gc.g1)$name, normalized = TRUE), na.rm = TRUE)
 
  ceb <- cluster_edge_betweenness(sim.g1, directed = FALSE)
  score.module <- modularity(sim.g1, ceb$membership) 
  
  return(list(score.between = score.between, score.crosstalk = score.crosstalk, # score.pathway excluded.
              score.clique = score.clique, score.close = score.close, score.module = score.module))

}

random.net.score <- function(Nnode, Nlink, Nrand){
  
  rand.cross.score <- as.data.frame(matrix(ncol=5, nrow = Nrand))
  colnames(rand.cross.score) <- c("score.between", "score.crosstalk", "score.clique", "score.degree", "score.module")
  
  for (ir in 1:Nrand){
    g.rand <- sample_gnm(Nnode, Nlink)
    V(g.rand)$name <- as.character(1:Nnode)
    gc.g.rand <- getgcc(g.rand)
    components.rand <- igraph::clusters(g.rand)
    N.single.node <- sum(components.rand$csize == 1) # Number of single node
    
    score.between <- max(degree(g.rand))
    score.crosstalk <- diameter(g.rand)/Nlink
    score.crosstalk <- score.crosstalk + score.crosstalk*runif(n=1,min=-1,max=1)*0.001 # Some cases have the same diameter values, need to give very small variances.
    score.clique <- transitivity(g.rand)
    score.close <- mean(closeness(g.rand, V(gc.g.rand)$name, normalized = TRUE), na.rm = TRUE)
    
    ceb <- cluster_edge_betweenness(g.rand, directed = FALSE)
    score.module <- modularity(g.rand, ceb$membership) 
    
    rand.cross.score[ir,] <- c(score.between, score.crosstalk, score.clique, score.close, score.module)
    
  }
  
  return(rand.cross.score = rand.cross.score)
  
}

columnMax <- function(data) sapply(data, max, na.rm = TRUE)
columnMin <- function(data) sapply(data, min, na.rm = TRUE)





####################### critical features ######################
#################################################################
#################################################################

# cancer.type.list <- c("KIDNEY", "COADREAD", "UCEC", "BRCA", "LGG", "LUSC", "OV", "LUAD", "LIHC", "STAD", "BLCA", "CESC")
cancer.type.list <- c("STAD", "LIHC", "OV", "BLCA", "UCEC")
# pathway.data <- readRDS(paste("TCGA-",cancer.type.list[1],"_best_features_nodes.rds",sep=""))

sig.path <- list()

for (i in 1:length(cancer.type.list)){
  # sig.path.gb <- read.table(paste("TCGA-",cancer.type.list[i],"_cut100_short_long_common.csv",sep=""), header = TRUE, sep = ",", quote="\"", dec=".",stringsAsFactors=FALSE)
  sig.path.gb <- read.table(paste("TCGA-",cancer.type.list[i],"_critical_features_short_long_common.csv",sep=""), header = TRUE, sep = ",", quote="\"", dec=".",stringsAsFactors=FALSE)
  sig.path.gb$minmax2 <- (sig.path.gb$minmax+0.1)/1.1
  
  ss<-sapply(1:length(sig.path.gb$variable), function(x) ifelse(length(ttt<-unlist(strsplit(sig.path.gb$variable[x],"P")))==2, 
                                                                paste0("P",ttt[2]),paste0("P",ttt[2]) ))
  tt<-sapply(1:length(sig.path.gb$variable), function(x) ifelse(length(ttt<-unlist(strsplit(sig.path.gb$variable[x],"P")))==2, 
                                                                paste0("P",ttt[2]),paste0("P",ttt[3]) ))
  
  sig.path.class <- data.frame(from = ss, to = tt, min_max = sig.path.gb$minmax2, class = sig.path.gb$classification)
  # sig.path.class$p.or.l <- sapply(1:dim(sig.path.class)[1], function(x) ifelse(sig.path.class$from[x] == sig.path.class$to[x], "path", "link"))
  
  sig.path[[i]] <- sig.path.class
  
}
names(sig.path) <- cancer.type.list






##############################################################################3
##############################################################################3
########## Degree distribution plot in log-log scale for 5 cancers
##############################################################################3
##############################################################################3

png(file = "Supple_power_law_fitting.png",
    width = 2500, height=2500, units = "px", res = 300)
par(mfrow=c(2,3))
# par(mar=c(0,0,0,0))
xmin_list <- c(2, 2, 5, 5, 5)

for (icancer in 1:length(cancer.type.list)){
  links.data <- sig.path[[as.character(cancer.type.list[icancer])]]
  
  importance.cut <- 0
  
  
  library(igraph)
  library(bc3net)
  sub.links.data <- dplyr::filter(links.data, min_max >= importance.cut)
  sub.links.data$p.or.l <- sapply(1:dim(sub.links.data)[1], function(x) ifelse(sub.links.data$from[x] == sub.links.data$to[x], "path", "link"))
  temp.link <- dplyr::filter(sub.links.data, p.or.l == "link")
  temp.pathway <- dplyr::filter(sub.links.data, p.or.l == "path")
  
  g1 <- graph_from_data_frame(sub.links.data, directed = F)
  sim.g1<-simplify(g1, remove.loop=T)
  
  degr <- degree(sim.g1)
  degr <- degr[degr!=0]
  kk <- table(degr)
  xx <-as.numeric(names(kk))

  # dd <- degree_distribution(sim.g1)
  # fitp <- fit_power_law(dd, implementation = c("plfit"))
  # plot(xx, dd[xx+1], log="xy", xlab = "Degree", ylab = "Density",
  #      main = paste0(cancer.type.list[icancer],"\n alpha = ",round(fitp$alpha, 2),", KS.p = ",round(fitp$KS.p,3)) )
  # lines(xx, xx^-fitp$alpha, col= "#b00606")
  
  # dd <- degree_distribution(sim.g1)
  # fitp <- fit_power_law(dd, implementation = c("plfit"))
  # plot(dd, log="xy",xlab = "Degree", ylab = "Density",
  #      main = paste0(cancer.type.list[icancer],"\n alpha = ",round(fitp$alpha, 2),", KS.p = ",round(fitp$KS.p,3)) )
  # lines(seq(dd), seq(dd)^-fitp$alpha, col= "#b00606")


  fitp <- fit_power_law(degr, implementation = c("plfit"), xmin = xmin_list[icancer])
  plot(xx,kk, log="xy", xlab = "Degree", ylab = "Count",
       main = paste0(cancer.type.list[icancer],"\n alpha = ",round(fitp$alpha, 2),", KS.p = ",round(fitp$KS.p,3)))
  lines(xx[xx>=xmin_list[icancer]], max(kk)*(xx[xx>=xmin_list[icancer]]/xx[xx==xmin_list[icancer]])^-fitp$alpha, col= "#b00606")

 
  
}

dev.off()







##############################################################################3
##############################################################################3
########## Plot radar plots of p-values for 5 cancers
##############################################################################3
##############################################################################3


Nrand <- 1000

vsn <- 10^(-100)

crosstalks.score <- as.data.frame(matrix(ncol=5, nrow = length(cancer.type.list)))
colnames(crosstalks.score)<-c("betweenness", "diameter", "transitivity", "closeness", "modularity")
rownames(crosstalks.score)<-cancer.type.list
r.crosstalks.score <- crosstalks.score
logp.crosstalks.score <- crosstalks.score
logp.crosstalks.score2 <- crosstalks.score
xmin_list <- c(2, 2, 5, 5, 5)

for (icancer in 1:length(cancer.type.list)){
  
  links.data <- sig.path[[as.character(cancer.type.list[icancer])]]
  
  importance.cut <- 0
  sct <- unlist(score_crosstalk_types(links.data, importance.cut))
  sct[is.na(sct)] <- 0 # NaN set to 0
  crosstalks.score[icancer,] <- sct
  
  g1 <- graph_from_data_frame(links.data, directed = F)
  g2<-simplify(g1, remove.loop=T)
  Nnode <- length(V(g2))
  Nlink <- length(E(g2))
  rns <- random.net.score(Nnode, Nlink, Nrand)
  rns[is.na(rns)] <- 0 # NaN set to 0
  r.crosstalks.score[icancer,]<- colMeans(rns)
  
  # g1 <- graph_from_data_frame(sub.links.data, directed = F)
  sim.g1<-simplify(g1, remove.loop=T)
  
  degr <- degree(sim.g1)
  degr <- degr[degr!=0]
  fit.pval <- fit_power_law(degr, implementation = c("plfit"), xmin = xmin_list[icancer])$KS.p
  
  for (iscore in 1:5){
    ### t-test
    greater.mu <- t.test(rns[,iscore], mu=crosstalks.score[icancer,iscore], alternative = "less")$p.value
    less.mu <- t.test(rns[,iscore], mu=crosstalks.score[icancer,iscore], alternative = "greater")$p.value
    logp.crosstalks.score[icancer,iscore] <- -log10(min(greater.mu, less.mu)+vsn)*ifelse(greater.mu<=less.mu, 1, -1)
    # p.crosstalks.score[i,iscore] <- t.test(rns[,iscore], mu=crosstalks.score[i,iscore])$p.value
    
    ### Monte Carlo p-value
    no.ge <- length(which(rns[,iscore] >= crosstalks.score[icancer,iscore])) # number of random cases greater than the observed network
    no.le <- length(which(rns[,iscore] <= crosstalks.score[icancer,iscore])) # number of random cases less than the observed network
    greater.mu2 <- (no.ge+1)/(Nrand+1)
    less.mu2 <- (no.le+1)/(Nrand+1)
    logp.crosstalks.score2[icancer,iscore] <- -log10(min(greater.mu2, less.mu2)+vsn)*ifelse(greater.mu2<=less.mu2, 1, -1)
    
  }
  
  logp.crosstalks.score[icancer, 1] <- fit.pval*200 # use power law p-value instead of degree p-value. need to scale up
  logp.crosstalks.score2[icancer, 1] <- fit.pval*6
}

save(logp.crosstalks.score, logp.crosstalks.score2, file = "network_index_t_test_Monte_p_value_n1000_240828.RData")
load("network_index_t_test_Monte_p_value_n1000_240828.RData")


########### PLOT ##############################################################

######################################
##### Network index by t-test p-value

all.p.score <- logp.crosstalks.score
png(file = "RadarPlot_network_metrics_ttest.png",
    width = 2048*2, height=2048, units = "px", res = 300)

library(ggradar)
plist <- list()
j<-0
for (icancer in 1:5){
  # for (icancer in 1:5){
  c.score <- all.p.score[icancer,]
  # rc.score <- all.rc.score[[icancer]]
  
  r.data <- rbind(c.score, c(10,-log10(0.05),-log10(0.05),-log10(0.05),-log10(0.05)),c(0,log10(0.05),log10(0.05),log10(0.05),log10(0.05))) # 0.05*200 or 0
  # r.data <- rbind(c.score, c(0.3,-log10(0.05),-log10(0.05),-log10(0.05),-log10(0.05)),c(0,log10(0.05),log10(0.05),log10(0.05),log10(0.05))) # 0.05*6 or 0
  
  # colnames(r.data) <- c("degree", "diameter", "transitivity", "closeness", "modularity")
  colnames(r.data) <- c("Deg.", "Dia.", "T", "C", "M")
  rownames(r.data) <- c("TCGA", "Random_H", "Random_L")
  r.data <- tibble::rownames_to_column(r.data, "group")
  # r.data$linkdensity <- r.data$linkdensity*300 # rescale for score.linkdensity
  j <- j+1
  plist[[j]]<-ggradar(
    r.data, 
    values.radar = c("-100", "0", "100"),
    grid.min = -100, grid.mid = 0, grid.max = 100,
    # values.radar = c("-3", "0", "3"),
    # grid.min = -3.01, grid.mid = 0, grid.max = 3.01,
    # Polygons
    group.line.width = 0.8, 
    group.point.size = 2,
    grid.label.size=3,
    axis.label.size=4,
    # group.colours = c("#d9d9d9","#8dd3c7","#b3de68", "#81b1d3","#fccde5","#feffb4" ), # alphabetical order
    # group.colours = c("#00AFBB", "#EE0099", "#FC4E07","#00BA38","#E7B800","#0000AC" ), # alphabetical order
    group.colours = c("grey80","grey90", "#FC4E07" ), # alphabetical order
    # Background and grid lines
    background.circle.colour = "white",
    gridline.mid.colour = "grey", plot.title = rownames(c.score), 
    plot.legend = FALSE, legend.text.size = 8, legend.position="right",
  )
  
  
}

library(gridExtra)
do.call("grid.arrange", c(plist, ncol=3))

dev.off()



##########################################
##### Network index by Monte Carlo p-value
all.p.score <- logp.crosstalks.score2

png(file = "RadarPlot_network_metrics_Monte_Carlo.png",
    width = 2048*2, height=2048, units = "px", res = 300)

    
library(ggradar)
plist <- list()
j<-0
for (icancer in 1:5){
  # for (icancer in 1:5){
  c.score <- all.p.score[icancer,]
  # rc.score <- all.rc.score[[icancer]]
  
  # r.data <- rbind(c.score, c(10,-log10(0.05),-log10(0.05),-log10(0.05),-log10(0.05)),c(0,log10(0.05),log10(0.05),log10(0.05),log10(0.05))) # 0.05*200 or 0
  r.data <- rbind(c.score, c(0.3,-log10(0.05),-log10(0.05),-log10(0.05),-log10(0.05)),c(0,log10(0.05),log10(0.05),log10(0.05),log10(0.05))) # 0.05*6 or 0
  
  # colnames(r.data) <- c("degree", "diameter", "transitivity", "closeness", "modularity")
  colnames(r.data) <- c("Deg.", "Dia.", "T", "C", "M")
  rownames(r.data) <- c("TCGA", "Random_H", "Random_L")
  r.data <- tibble::rownames_to_column(r.data, "group")
  # r.data$linkdensity <- r.data$linkdensity*300 # rescale for score.linkdensity
  j <- j+1
  plist[[j]]<-ggradar(
    r.data, 
    # values.radar = c("-100", "0", "100"),
    # grid.min = -100, grid.mid = 0, grid.max = 100,
    values.radar = c("-3", "0", "3"),
    grid.min = -3.01, grid.mid = 0, grid.max = 3.01,
    # Polygons
    group.line.width = 0.8, 
    group.point.size = 2,
    grid.label.size=3,
    axis.label.size=4,
    # group.colours = c("#d9d9d9","#8dd3c7","#b3de68", "#81b1d3","#fccde5","#feffb4" ), # alphabetical order
    # group.colours = c("#00AFBB", "#EE0099", "#FC4E07","#00BA38","#E7B800","#0000AC" ), # alphabetical order
    group.colours = c("grey80","grey90", "#FC4E07" ), # alphabetical order
    # Background and grid lines
    background.circle.colour = "white",
    gridline.mid.colour = "grey", plot.title = rownames(c.score), 
    plot.legend = FALSE, legend.text.size = 8, legend.position="right",
  )


}

library(gridExtra)
do.call("grid.arrange", c(plist, ncol=3))

dev.off()


