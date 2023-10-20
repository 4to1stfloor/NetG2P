
## UPDATED by DKS
## 1. Changing the vertex size depending on the importance values in pan-cancer network.
## 2. Fixed the following 2 lines.
##    edge.g1 <- as_data_frame(g1)
##    V(line.g1)$name <- sapply(1:dim(edge.g1)[1], function(x) paste(edge.g1[x,1], edge.g1[x,2],sep=""))
##    -->
##    V(line.g1)$name <- sapply(1:dim(sub.links.data)[1], function(x) paste(sub.links.data[x,1], sub.links.data[x,2],sep="")) 

## change the top250 for the issue on survival analysis.

rm(list = ls()) # remove all variables

normalize <- function(x) {
  return((x-min(x))/(max(x)-min(x)))
}

data.folder <- "/mnt/gluster_server/data/NetGPT/top250/"
setwd(data.folder)

library(igraph)
library(dplyr)

####################### cut 250 ######################
#################################################################
#################################################################

cancer.type.list <- c("KIDNEY", "UCEC", "BRCA", "LGG", "LUSC", "OV", "LUAD", "LIHC", "STAD", "BLCA", "CESC")

sig.path <- list()

for (i in 1:length(cancer.type.list)){
  sig.path.gb <- read.table(paste("TCGA-",cancer.type.list[i],"_top250_from_ml.csv",sep=""), header = TRUE, sep = ",", quote="\"", dec=".",stringsAsFactors=FALSE)
  
  
  ss<-sapply(1:length(sig.path.gb$variable), function(x) ifelse(length(ttt<-unlist(strsplit(sig.path.gb$variable[x],"P")))==2, 
                                                                paste0("P",ttt[2]),paste0("P",ttt[2]) ))
  tt<-sapply(1:length(sig.path.gb$variable), function(x) ifelse(length(ttt<-unlist(strsplit(sig.path.gb$variable[x],"P")))==2, 
                                                                paste0("P",ttt[2]),paste0("P",ttt[3]) ))
  sig.path.gb = sig.path.gb %>% mutate(minmax = normalize(relative_importance))
  
  # sig.path.class <- data.frame(from = ss, to = tt, min_max = sig.path.gb$minmax, class = sig.path.gb$classification)
  sig.path.class <- data.frame(from = ss, to = tt, min_max = sig.path.gb$minmax)
  # sig.path.class$p.or.l <- sapply(1:dim(sig.path.class)[1], function(x) ifelse(sig.path.class$from[x] == sig.path.class$to[x], "path", "link"))
  
  sig.path[[i]] <- sig.path.class
  
}
names(sig.path) <- cancer.type.list


# pheatmap(sub.heat)

################## Link graph #########################
################## pan-cancer network ###############
# sub.cancer.type.list <- c("KIDNEY", "UCEC", "BRCA", "LGG", "OV", "LUAD", "LIHC", "STAD", "BLCA", "CESC") # LUSC excluded
sub.cancer.type.list <- c("UCEC", "BRCA", "LGG", "OV", "LUAD", "LIHC", "STAD", "BLCA", "CESC") # KIDNEY, COADREDAD, LUSC excluded
# sub.cancer.type.list <- c("UCEC", "BRCA", "OV", "LIHC", "STAD", "BLCA") # 

default.link.weight <- 1 
within.group.weight <- 100 ## link weight between a virtual node and other nodes in a cancer group 
btw.group.weight <- 0.01 ## link weight between virtual nodes

group.ids <- list()
# icancer = 1
for (icancer in 1:length(sub.cancer.type.list)){
  # for (icancer in 2:5){
  sub.links.data <- sig.path[[which(names(sig.path) == sub.cancer.type.list[[icancer]])]]
  g1 <- graph_from_data_frame(sub.links.data, directed = F)
  line.g1 <- make_line_graph(g1)
  
  V(line.g1)$name <- sapply(1:dim(sub.links.data)[1], function(x) {
    if (sub.links.data[x, 1] == sub.links.data[x, 2]) {
      # print(sub.links.data[x, 1])
      # print(sub.links.data[x, 2])
      return(sub.links.data[x, 1])
    } else {
      return(paste(sub.links.data[x, 1], sub.links.data[x, 2], sep = ""))
    }
  })
  
  edge.line.g1 <- igraph::as_data_frame(line.g1)
  edge.line.g1$weight <- rep(default.link.weight, dim(edge.line.g1)[1]) # add weight = 1 for all links
  ## add virtual node
  virt<-data.frame(from = rep(sub.cancer.type.list[icancer],length(V(line.g1)$name)), to = V(line.g1)$name, weight = rep(within.group.weight,length(V(line.g1)$name)))
  v.edge.line.g1 <- rbind(edge.line.g1, virt)
  
  if (icancer == 1) {
    total.line.g = v.edge.line.g1
  } else {
    total.line.g <- rbind(total.line.g, v.edge.line.g1)
  }
  
  group.ids[[icancer]] <- unique(c(edge.line.g1$from, edge.line.g1$to, sub.cancer.type.list[icancer]))
  
}
names(group.ids) <- sub.cancer.type.list

total.line.graph <- dplyr::distinct(total.line.g) # remove duplicated links

temp.g <- as.data.frame(t(combn(sub.cancer.type.list, 2))) # add links between virtual nodes
temp.g$weight <- rep(btw.group.weight, dim(temp.g)[1]) # add their weight
colnames(temp.g) <- c("from", "to", "weight")
total.line.graph.2 <- rbind(total.line.graph, temp.g)

total.g <- graph_from_data_frame(total.line.graph.2, directed = FALSE)

# setting colors
library(RColorBrewer)
# group_color <- brewer.pal(length(group.ids), 'Set3') # 'Set3' 
library(ggsci)
group_color <- pal_npg("nrc")(length(group.ids))
# the fill gets an additional alpha value for transparency:
# group_color_fill <- paste0(group_color, '20')
group_color_fill <-pal_npg("nrc", alpha = 0.2)(length(group.ids))

## Remove all the virtual nodes and related links
lay <- layout_nicely(total.g)
rownames(lay) <- V(total.g)$name
total.g <- total.g - vertices(sub.cancer.type.list)
lay <- lay[which(!rownames(lay) %in% sub.cancer.type.list),]
group.ids.2 <- lapply(1:length(group.ids), function(x) setdiff(group.ids[[x]], sub.cancer.type.list[x]))
names(group.ids.2) <- names(group.ids)


divided_list = list()
cancer_type= names(group.ids.2)

tmp_boolen = c()

for (each_cancer in cancer_type) {
  for (cancer_features in group.ids.2[[each_cancer]]) {
    
    split_features <- strsplit(cancer_features, "P")
    feature1 <- paste0("P", split_features[[1]][2])
    feature2 <- paste0("P", split_features[[1]][3])
    if (feature1 == feature2) {
      cancer_features_for_minmax = feature1
      
      for (all_cancer in cancer_type) {
        
        importance_cancer <- read.table(paste("TCGA-",all_cancer,"_top250_from_ml.csv",sep=""), header = TRUE, sep = ",", quote="\"", dec=".",stringsAsFactors=FALSE)
        importance_cancer = importance_cancer %>% mutate(minmax = normalize(relative_importance))
        
        if (cancer_features %in% group.ids.2[[all_cancer]]) {
          
          if (importance_cancer[which(cancer_features_for_minmax == importance_cancer$variable),]$minmax == 0) {
            tmp_boolen = c(tmp_boolen ,1.0e-10)
          }else {
            tmp_boolen = c(tmp_boolen ,importance_cancer[which(cancer_features_for_minmax == importance_cancer$variable),]$minmax)
          }
          
          # tmp_boolen = c(tmp_boolen , 1)
        } else {
          tmp_boolen = c(tmp_boolen , 0)
        }
        
      }
      
      divided_list[[cancer_features_for_minmax]] = tmp_boolen
      tmp_boolen = c()
      
    } else {
      cancer_features_for_minmax = cancer_features
      
      for (all_cancer in cancer_type) {
        
        importance_cancer <- read.table(paste("TCGA-",all_cancer,"_top250_from_ml.csv",sep=""), header = TRUE, sep = ",", quote="\"", dec=".",stringsAsFactors=FALSE)
        importance_cancer = importance_cancer %>% mutate(minmax = normalize(relative_importance))
        
        if (cancer_features %in% group.ids.2[[all_cancer]]) {
          
          if (importance_cancer[which(cancer_features_for_minmax == importance_cancer$variable),]$minmax == 0) {
            tmp_boolen = c(tmp_boolen ,1.0e-10)
          }else {
            tmp_boolen = c(tmp_boolen ,importance_cancer[which(cancer_features_for_minmax == importance_cancer$variable),]$minmax)
          }
          
          # tmp_boolen = c(tmp_boolen , 1)
        } else {
          tmp_boolen = c(tmp_boolen , 0)
        }
        
      }
      
      divided_list[[cancer_features_for_minmax]] = tmp_boolen
      tmp_boolen = c()
      
    }
    
  }
}


total.g_filt = total.g - names(V(total.g))[which(!names(V(total.g)) %in% names(divided_list))]
lay_filt <- lay[which(rownames(lay) %in% names(divided_list)),]

order_name = V(total.g_filt)$name
tbl_count = table(unlist(group.ids.2))
tbl_count_ordered <- as.numeric(tbl_count[order(match(names(tbl_count), order_name))])

divided_list_filt <- divided_list[order(match(names(divided_list), names(V(total.g_filt))))]

plot(total.g_filt, 
     layout=lay_filt,  
     vertex.size = tbl_count_ordered, 
     vertex.label=NA,
     edge.color = rgb(0.5, 0.5, 0.5, 0.2),
     vertex.shape="pie",
     vertex.pie=divided_list_filt,
     mark.groups = group.ids.2,
     mark.col = group_color_fill,
     mark.border = group_color,
     vertex.pie.color= list(group_color))


legend('topright', legend=names(group.ids.2), 
       col = group_color,
       pch=15, bty = "n", pt.cex = 1.5, cex = 0.8,
       text.col = "black", horiz = FALSE)


#################################################################
####################### only critical features###################
#################################################################


sig.path_critical <- list()

for (i in 1:length(cancer.type.list)){
  each_sig_path_gb <- read.table(paste("TCGA-",cancer.type.list[i],"_critical_features.csv",sep=""), header = TRUE, sep = ",", quote="\"", dec=".",stringsAsFactors=FALSE)
  
  ss_cri<-sapply(1:length(each_sig_path_gb$variable), function(x) ifelse(length(ttt<-unlist(strsplit(each_sig_path_gb$variable[x],"P")))==2, 
                                                                         paste0("P",ttt[2]),paste0("P",ttt[2]) ))
  tt_cri <-sapply(1:length(each_sig_path_gb$variable), function(x) ifelse(length(ttt<-unlist(strsplit(each_sig_path_gb$variable[x],"P")))==2, 
                                                                          paste0("P",ttt[2]),paste0("P",ttt[3]) ))
  
  each_sig_path_gb = each_sig_path_gb %>% mutate(minmax = normalize(relative_importance))
  
  sig.path.class_cri <- data.frame(from = ss_cri, to = tt_cri, min_max = each_sig_path_gb$minmax)
  # sig.path.class$p.or.l <- sapply(1:dim(sig.path.class)[1], function(x) ifelse(sig.path.class$from[x] == sig.path.class$to[x], "path", "link"))
  
  sig.path_critical[[i]] <- sig.path.class_cri
  
}
names(sig.path_critical) <- cancer.type.list

# sub.cancer.type.list <- c("KIDNEY", "COADREAD", "UCEC", "BRCA", "LGG", "OV", "LUAD", "LIHC", "STAD", "BLCA", "CESC") # LUSC excluded
# sub.cancer.type.list <- c("KIDNEY", "UCEC", "BRCA", "LGG", "OV", "LUAD", "LIHC", "STAD", "BLCA", "CESC") # LUSC excluded
sub.cancer.type.list <- c("UCEC", "BRCA", "LGG", "OV", "LUAD","LIHC","STAD", "BLCA", "CESC") # KIDNEY, COADREDAD, LUSC excluded
# sub.cancer.type.list <- c("UCEC", "BRCA", "OV", "LIHC", "STAD", "BLCA") # 

default.link.weight <- 1 
within.group.weight <- 1000 ## link weight between a virtual node and other nodes in a cancer group 
btw.group.weight <- 0.00001 ## link weight between virtual nodes


group.ids_cri <- list()

for (icancer in 1:length(sub.cancer.type.list)){
  # for (icancer in 2:5){
  sub.links.data_cri <- sig.path_critical[[which(names(sig.path_critical) == sub.cancer.type.list[[icancer]])]]
  g1_cri <- graph_from_data_frame(sub.links.data_cri, directed = F)
  line.g1_cri <- make_line_graph(g1_cri)
  
  V(line.g1_cri)$name <- sapply(1:dim(sub.links.data_cri)[1], function(x) {
    if (sub.links.data_cri[x, 1] == sub.links.data_cri[x, 2]) {
      
      return(sub.links.data_cri[x, 1])
    } else {
      return(paste(sub.links.data_cri[x, 1], sub.links.data_cri[x, 2], sep = ""))
    }
  })
  
  
  
  edge.line.g1_cri <- igraph::as_data_frame(line.g1_cri)
  edge.line.g1_cri$weight <- rep(default.link.weight, dim(edge.line.g1_cri)[1]) # add weight = 1 for all links
  ## add virtual node
  virt_cri<-data.frame(from = rep(sub.cancer.type.list[icancer],length(V(line.g1_cri)$name)), to = V(line.g1_cri)$name, weight = rep(within.group.weight,length(V(line.g1_cri)$name)))
  v.edge.line.g1_cri <- rbind(edge.line.g1_cri, virt_cri)
  
  if (icancer == 1) {
    total.line.g_cri = v.edge.line.g1_cri
  } else {
    total.line.g_cri <- rbind(total.line.g_cri, v.edge.line.g1_cri)
  }
  
  group.ids_cri[[icancer]] <- unique(c(edge.line.g1_cri$from, edge.line.g1_cri$to, sub.cancer.type.list[icancer]))
  
}

names(group.ids_cri) <- sub.cancer.type.list

total.line.graph_cri <- dplyr::distinct(total.line.g_cri) # remove duplicated links
temp.g_cri <- as.data.frame(t(combn(sub.cancer.type.list, 2))) # add links between virtual nodes
temp.g_cri$weight <- rep(btw.group.weight, dim(temp.g_cri)[1]) # add their weight
colnames(temp.g_cri) <- c("from", "to", "weight")
total.line.graph_cri.2 <- rbind(total.line.graph_cri, temp.g_cri)

total.g_cri <- graph_from_data_frame(total.line.graph_cri.2, directed = FALSE)

## Remove all the virtual nodes and related links
lay_cri <- layout_nicely(total.g_cri)
rownames(lay_cri) <- V(total.g_cri)$name
total.g_cri <- total.g_cri - vertices(sub.cancer.type.list)
lay_cri <- lay_cri[-which(rownames(lay_cri) %in% sub.cancer.type.list),]
group.ids_cri.2 <- lapply(1:length(group.ids_cri), function(x) setdiff(group.ids_cri[[x]], sub.cancer.type.list[x]))
names(group.ids_cri.2) <- names(group.ids_cri)

divided_list_cri = list()
cancer_type_cri = names(group.ids_cri.2)

tmp_boolen_cri = c()
for (each_cancer in cancer_type_cri) {
  for (cancer_features in group.ids_cri.2[[each_cancer]]) {
    
    split_features <- strsplit(cancer_features, "P")
    feature1 <- paste0("P", split_features[[1]][2])
    feature2 <- paste0("P", split_features[[1]][3])
    if (feature1 == feature2) {
      cancer_features_for_minmax = feature1
      
      for (all_cancer in cancer_type) {
        
        importance_cancer <- read.table(paste("TCGA-",all_cancer,"_critical_features.csv",sep=""), 
                                        header = TRUE, sep = ",", quote="\"", dec=".",stringsAsFactors=FALSE)
        importance_cancer = importance_cancer %>% mutate(minmax = normalize(relative_importance))
        
        if (cancer_features %in% group.ids_cri.2[[all_cancer]] & cancer_features_for_minmax %in% importance_cancer$variable ) {
          
          if (importance_cancer[which( importance_cancer$variable == cancer_features_for_minmax ),]$minmax == 0) {
            tmp_boolen_cri = c(tmp_boolen_cri ,1.0e-10)
          }else {
            tmp_boolen_cri = c(tmp_boolen_cri ,importance_cancer[which(cancer_features_for_minmax == importance_cancer$variable),]$minmax)
          }
          
          # tmp_boolen = c(tmp_boolen , 1)
        } else {
          tmp_boolen_cri = c(tmp_boolen_cri , 0)
        }
        
      }
      
      divided_list_cri[[cancer_features_for_minmax]] = tmp_boolen_cri
      tmp_boolen_cri = c()
      
    } else {
      cancer_features_for_minmax = cancer_features
      for (all_cancer in cancer_type) {
        
        importance_cancer <- read.table(paste("TCGA-",all_cancer,"_critical_features.csv",sep=""), 
                                        header = TRUE, sep = ",", quote="\"", dec=".",stringsAsFactors=FALSE)
        
        importance_cancer = importance_cancer %>% mutate(minmax = normalize(relative_importance))
        
        if (cancer_features %in% group.ids_cri.2[[all_cancer]] & cancer_features_for_minmax %in% importance_cancer$variable ) {
          
          if (importance_cancer[which(importance_cancer$variable == cancer_features_for_minmax),]$minmax == 0) {
            tmp_boolen_cri = c(tmp_boolen_cri ,1.0e-10)
          }else {
            tmp_boolen_cri = c(tmp_boolen_cri ,importance_cancer[which(cancer_features_for_minmax == importance_cancer$variable),]$minmax)
          }
          
          # tmp_boolen = c(tmp_boolen , 1)
        } else {
          tmp_boolen_cri = c(tmp_boolen_cri , 0)
        }
        
      }
      
      divided_list_cri[[cancer_features_for_minmax]] = tmp_boolen_cri
      tmp_boolen_cri = c()
    }
    
  }
}

total.g_filt_cri = total.g_cri - names(V(total.g_cri))[which(!names(V(total.g_cri)) %in% names(divided_list_cri))]

without_pathway = rownames(lay_cri)[which(!rownames(lay_cri) %in% names(V(total.g_filt_cri)) )]

lay_filt_cri <- lay_cri[!(rownames(lay_cri) %in% without_pathway), ]


order_name_cri = V(total.g_filt_cri)$name
tbl_count_cri = table(unlist(group.ids_cri.2))
tbl_count_ordered_cri <- as.numeric(tbl_count_cri[order(match(names(tbl_count_cri), order_name_cri))])

divided_list_filt_cri <- divided_list_cri[order(match(names(divided_list_cri), names(V(total.g_filt_cri))))]
# names(V(total.g_filt_cri)) %in% rownames(lay_filt_cri)

plot(total.g_filt_cri,
     layout=lay_filt_cri,
     vertex.size = tbl_count_ordered_cri,
     vertex.label=NA,
     edge.color = rgb(0.5, 0.5, 0.5, 0.2),
     vertex.shape="pie",
     vertex.pie=divided_list_filt_cri,
     mark.groups = group.ids_cri.2,
     mark.col = group_color_fill,
     mark.border = group_color,
     vertex.pie.color= list(group_color)
)

legend('topright', legend=names(group.ids_cri.2), 
       col = group_color,
       pch=15, bty = "n", pt.cex = 1.5, cex = 0.8,
       text.col = "black", horiz = FALSE)


#################################################################
##### critical features on a cut 250 network#####################
#################################################################

## Remove all the virtual nodes and related links
# lay_cri <- layout_nicely(total.g_cri)
# rownames(lay_cri) <- V(total.g_cri)$name
# total.g_cri <- total.g_cri - vertices(sub.cancer.type.list)
# lay_cri <- lay_cri[-which(rownames(lay_cri) %in% sub.cancer.type.list),]
# group.ids_cri.2 <- lapply(1:length(group.ids_cri), function(x) setdiff(group.ids_cri[[x]], sub.cancer.type.list[x]))
# names(group.ids_cri.2) <- names(group.ids_cri)

communal_tmp = tbl_count[names(tbl_count) %in% names(tbl_count_cri)] +2
only_cut_tmp = tbl_count[!names(tbl_count) %in% names(tbl_count_cri)]


tbl_count_communal = c(communal_tmp,only_cut_tmp)

order_name = V(total.g_filt)$name

tbl_count_communal_ordered <- as.numeric(tbl_count_communal[order(match(names(tbl_count_communal), order_name))])

plot(total.g_filt, 
     layout=lay_filt,  
     vertex.size = tbl_count_communal_ordered, 
     vertex.label=NA,
     edge.color = rgb(0.5, 0.5, 0.5, 0.2),
     vertex.shape="pie",
     vertex.pie=divided_list_filt,
     mark.groups = group.ids.2,
     mark.col = group_color_fill,
     mark.border = group_color,
     vertex.pie.color= list(group_color))


legend('topright', legend=names(group.ids.2), 
       col = group_color,
       pch=15, bty = "n", pt.cex = 1.5, cex = 0.8,
       text.col = "black", horiz = FALSE)
# 

