library(org.Hs.eg.db)
library(data.table)
library(scater)
library(scran)
library(pheatmap)
library(devtools)
library(data.table)
library(scuttle)
library(scran)

# zscore normalization

tcga.calc.zscore = function(sce, target.genes){
  message("Calculating z-score with respect to all diploid cells. Version 2023.05.03")
  common.genes = intersect(rownames(sce), target.genes)
  if (length(common.genes) == 0) {
    stop("None of the target genes found in sce. Please check your nomenclature")
  } else if (length(common.genes) != length(target.genes)) {
    message("Some of the genes from query does not exist in this cancer type. It will result in NAs")
    message("Missing genes are: ", paste(setdiff(target.genes, common.genes), collapse = ", "))
  }
  sce.sub = subset(sce, rownames(sce) %in% common.genes,)
  #i am not checking assay names.
  count.mat = assay(sce.sub, 1)
  cnv.mat = assay(sce.sub, 3)
  z.mat = matrix(data = NA, nrow = nrow(count.mat), ncol = ncol(count.mat))
  colnames(z.mat) = colnames(count.mat)
  rownames(z.mat) = rownames(count.mat)
  for (i in 1:nrow(count.mat)) {
    idx.di = which(cnv.mat[i,] == 2)
    query.mean = mean(count.mat[i, idx.di], na.rm = T)
    query.sd = sd(count.mat[i, idx.di], na.rm = T)
    z.mat[i,] = (count.mat[i,] - query.mean)/query.sd
  }
  return(z.mat)
}
filepath = "/home/seokwon/nas/"
ref_path = paste0(filepath, "99.reference/")

Cancerlist = dir(paste0(filepath, "/00.data/filtered_TCGA/"))
setwd("~/nas/04.Results/compare_surv/pathway/")

types = c("exp","mut","")

# for all
for (type in types) {
  if (type == "") {
    
    for (num_CancerType in Cancerlist) {
      main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
      CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
      
      # call input
      data_tc_link = readRDS(file = paste0(main.path_tc,"/",CancerType,"_pathwaylink_all_log.rds"))
      data_tc_each = readRDS(file = paste0(main.path_tc,"/",CancerType,"_pathwayeach_all_log.rds"))
      
      data_tc_link_filt = data_tc_link[,which(colSums(data_tc_link[-ncol(data_tc_link)]) != 0)]
      data_tc_each_filt = data_tc_each[,which(colSums(data_tc_each[-ncol(data_tc_each)]) != 0)]
      
      data_tc_link_filt$vitalstatus = data_tc_link$vitalstatus
      data_tc_each_filt$vitalstatus = data_tc_each$vitalstatus
      
      data_tc_link_filt = data_tc_link_filt[intersect(rownames(data_tc_link_filt), rownames(data_tc_each_filt)),]
      data_tc_each_filt = data_tc_each_filt[intersect(rownames(data_tc_link_filt), rownames(data_tc_each_filt)),]
      
      if (all.equal(rownames(data_tc_link_filt), rownames(data_tc_each_filt))) {
        data_tc_link_filt$vitalstatus = NULL
        data_tc = cbind(data_tc_link_filt,data_tc_each_filt)
      }
      
      data_tc = data_tc[which(!is.na(data_tc$vitalstatus)),]
      data_tc = data_tc[which(!data_tc$vitalstatus %in% c("Not Reported")),]
      
      data_tc_tmp = data_tc[,which(!colnames(data_tc) %in% c("vitalstatus"))]
      data_tc_meta = data_tc[,c("vitalstatus")]
      
      # # Convert the matrix to a numeric matrix
      # data_tc_numeric <- matrix(as.numeric(unlist(data_tc_tmp)), nrow = nrow(data_tc_tmp))
      # 
      # # Convert the numeric matrix to a vector
      # vec <- as.vector(data_tc_numeric)
      # 
      # # Create the histogram
      # hist(vec)

      data_tc_tmp[data_tc_tmp > 15] = 15
      
      data_tc_filt = cbind(data_tc_tmp,vitalstatus = data_tc_meta)
      
      # 
      data_tc_filt_df = data_tc_filt[order(data_tc_filt$vitalstatus),]
      annotation_df <- data.frame(vitalstatus = data_tc_filt_df$vitalstatus)
      rownames(annotation_df) <- rownames(data_tc_filt_df)
      
      # Create a named color vector for the unique values of vital_status
      vital_status_colors <- c("Alive" = "yellow", "Dead" = "black")
      names(vital_status_colors) <- unique(annotation_df$vitalstatus)
      
      png(filename = paste0(CancerType,"_net_pathway_clusterT_surv.png"),
          width = 25, height = 25,  units = "cm" ,pointsize = 12,
          bg = "white", res = 1200, family = "")
      
      tmp = pheatmap::pheatmap(as.matrix(t(data_tc_filt_df[,which(!colnames(data_tc_filt_df) %in% c("vitalstatus"))])),
                               annotation_col = annotation_df,
                               annotation_colors = list(vitalstatus = vital_status_colors),
                               cluster_cols = T)
      
      print(tmp)
      
      dev.off()
      
      png(filename = paste0(CancerType,"_net_pathway_clusterF_surv.png"),
          width = 25, height = 25,  units = "cm" ,pointsize = 12,
          bg = "white", res = 1200, family = "")
      
      tmp2 = pheatmap::pheatmap(as.matrix(t(data_tc_filt_df[,which(!colnames(data_tc_filt_df) %in% c("vitalstatus"))])),
                                annotation_col = annotation_df,
                                annotation_colors = list(vitalstatus = vital_status_colors),
                                cluster_cols = F)
      print(tmp2)
      
      dev.off()
      
      png(filename = paste0(CancerType,"_net_pathway_surv_complex.png"),
          width = 25, height = 25,  units = "cm" ,pointsize = 12,
          bg = "white", res = 1200, family = "")
      tmp3 = ComplexHeatmap::pheatmap(as.matrix(t(data_tc_filt_df[,which(!colnames(data_tc_filt_df) %in% c("vitalstatus"))])),
                                      column_split = annotation_df$vitalstatus,
                                      annotation_col = annotation_df,
                                      annotation_colors = list(vitalstatus = vital_status_colors),
                                      cluster_cols = T)
      
      print(tmp3)
      
      dev.off()
      
      remove(tmp,tmp2,tmp3,data_tc_filt_df,data_tc_filt,data_tc)
    }
  
   } else {
     
     for (num_CancerType in Cancerlist) {
       
       main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
       CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
       
       # call input
       data_tc_link = readRDS(file = paste0(main.path_tc,"/",CancerType,"_",type,"_pathwaylink_all_log.rds"))
       data_tc_each = readRDS(file = paste0(main.path_tc,"/",CancerType,"_",type,"_pathwayeach_all_log.rds"))
       
       data_tc_link_filt = data_tc_link[,which(colSums(data_tc_link[-ncol(data_tc_link)]) != 0)]
       data_tc_each_filt = data_tc_each[,which(colSums(data_tc_each[-ncol(data_tc_each)]) != 0)]
       
       data_tc_link_filt$vitalstatus = data_tc_link$vitalstatus
       data_tc_each_filt$vitalstatus = data_tc_each$vitalstatus
       
       data_tc_link_filt = data_tc_link_filt[intersect(rownames(data_tc_link_filt), rownames(data_tc_each_filt)),]
       data_tc_each_filt = data_tc_each_filt[intersect(rownames(data_tc_link_filt), rownames(data_tc_each_filt)),]
       
       if (all.equal(rownames(data_tc_link_filt), rownames(data_tc_each_filt))) {
         data_tc_link_filt$vitalstatus = NULL
         data_tc = cbind(data_tc_link_filt,data_tc_each_filt)
       }
       data_tc = data_tc[which(!is.na(data_tc$vitalstatus)),]
       data_tc = data_tc[which(!data_tc$vitalstatus %in% c("Not Reported")),]
       
       data_tc_tmp = data_tc[,which(!colnames(data_tc) %in% c("vitalstatus"))]
       data_tc_meta = data_tc[,c("vitalstatus")]
       # # Convert the matrix to a numeric matrix
       # data_tc_numeric <- matrix(as.numeric(unlist(data_tc_tmp)), nrow = nrow(data_tc_tmp))
       # 
       # # Convert the numeric matrix to a vector
       # vec <- as.vector(data_tc_numeric)
       # 
       # # Create the histogram
       # hist(vec)
       # 
       if (type == "exp") {
         data_tc_tmp[data_tc_tmp > 4] = 4
       } else {
         data_tc_tmp[data_tc_tmp > 1] = 1
       }
       
       data_tc_filt = cbind(data_tc_tmp,vitalstatus = data_tc_meta)

       # 
       data_tc_filt_df = data_tc_filt[order(data_tc_filt$vitalstatus),]
       annotation_df <- data.frame(vitalstatus = data_tc_filt_df$vitalstatus)
       rownames(annotation_df) <- rownames(data_tc_filt_df)
      
       # Create a named color vector for the unique values of vital_status
       vital_status_colors <- c("Alive" = "yellow", "Dead" = "black")
       names(vital_status_colors) <- unique(annotation_df$vitalstatus)
       
       png(filename = paste0(CancerType,"_",type,"_pathway_clusterT_surv.png"),
           width = 25, height = 25,  units = "cm" ,pointsize = 12,
           bg = "white", res = 1200, family = "")
       
       tmp = pheatmap::pheatmap(as.matrix(t(data_tc_filt_df[,which(!colnames(data_tc_filt_df) %in% c("vitalstatus"))])),
                                annotation_col = annotation_df,
                                annotation_colors = list(vitalstatus = vital_status_colors),
                                cluster_cols = T)
       
       print(tmp)
       
       dev.off()
       
       png(filename = paste0(CancerType,"_",type,"_pathway_clusterF_surv.png"),
           width = 25, height = 25,  units = "cm" ,pointsize = 12,
           bg = "white", res = 1200, family = "")
       
       tmp2 = pheatmap::pheatmap(as.matrix(t(data_tc_filt_df[,which(!colnames(data_tc_filt_df) %in% c("vitalstatus"))])),
                                annotation_col = annotation_df,
                                annotation_colors = list(vitalstatus = vital_status_colors),
                                cluster_cols = F)
       print(tmp2)
       
       dev.off()
       
       png(filename = paste0(CancerType,"_",type,"_pathway_surv_complex.png"),
           width = 25, height = 25,  units = "cm" ,pointsize = 12,
           bg = "white", res = 1200, family = "")
       tmp3 = ComplexHeatmap::pheatmap(as.matrix(t(data_tc_filt_df[,which(!colnames(data_tc_filt_df) %in% c("vitalstatus"))])),
                                       column_split = annotation_df$vitalstatus,
                                       annotation_col = annotation_df,
                                       annotation_colors = list(vitalstatus = vital_status_colors),
                                       cluster_cols = T)

       print(tmp3)
       
       dev.off()
       
       remove(tmp,tmp2,tmp3,data_tc_filt_df,data_tc_filt,data_tc)
     }
     
     
   }
}

