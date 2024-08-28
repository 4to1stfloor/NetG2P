library(data.table)
library(org.Hs.eg.db)
library(data.table)
library(scater)
library(scran)
library(pheatmap)
library(devtools)
library(data.table)
library(scuttle)
library(scran)
library(RColorBrewer)

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
cli = fread("~/nas/99.reference/all_clin_indexed.csv")
setwd("~/nas/04.Results/subtypes/tmn/pathway/")

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
      
      if (CancerType == "TCGA-COADREAD") {
        cli_surv = cli[cli$project %in% c("TCGA-COAD","TCGA-READ"),
                       c("submitter_id",
                         "ajcc_pathologic_t",
                         "ajcc_pathologic_n",
                         "ajcc_pathologic_m")]
        
      } else if (CancerType == "TCGA-KIDNEY") {
        cli_surv = cli[cli$project %in% c("TCGA-KIRP","TCGA-KIRC","TCGA-KICH"),
                       c("submitter_id",
                         "ajcc_pathologic_t",
                         "ajcc_pathologic_n",
                         "ajcc_pathologic_m")]
        
      } else {
        cli_surv = cli[cli$project == CancerType,
                       c("submitter_id",
                         "ajcc_pathologic_t",
                         "ajcc_pathologic_n",
                         "ajcc_pathologic_m")]
      }
      
      if (sum(!is.na(cli_surv$ajcc_pathologic_t)) == 0 || sum(!is.na(cli_surv$ajcc_pathologic_n)) == 0 || sum(!is.na(cli_surv$ajcc_pathologic_m)) == 0) {
        print(paste0( "there are not tmn data in ",CancerType))
        next
      }
      
      cli_surv_filtered = cli_surv[which(!is.na(cli_surv$ajcc_pathologic_t)),]
      cli_surv_filtered = cli_surv_filtered[which(!is.na(cli_surv_filtered$ajcc_pathologic_n)),]
      cli_surv_filtered = cli_surv_filtered[which(!is.na(cli_surv_filtered$ajcc_pathologic_m)),]
      
      cli_surv_filtered$ajcc_pathologic_t <- ifelse(grepl("^T1", cli_surv_filtered$ajcc_pathologic_t), "T1", cli_surv_filtered$ajcc_pathologic_t)
      cli_surv_filtered$ajcc_pathologic_t <- ifelse(grepl("^T2", cli_surv_filtered$ajcc_pathologic_t), "T2", cli_surv_filtered$ajcc_pathologic_t)
      cli_surv_filtered$ajcc_pathologic_t <- ifelse(grepl("^T3", cli_surv_filtered$ajcc_pathologic_t), "T3", cli_surv_filtered$ajcc_pathologic_t)
      cli_surv_filtered$ajcc_pathologic_t <- ifelse(grepl("^T4", cli_surv_filtered$ajcc_pathologic_t), "T4", cli_surv_filtered$ajcc_pathologic_t)
      
      cli_surv_filtered$ajcc_pathologic_n <- ifelse(grepl("^N1", cli_surv_filtered$ajcc_pathologic_n), "N1", cli_surv_filtered$ajcc_pathologic_n)
      cli_surv_filtered$ajcc_pathologic_n <- ifelse(grepl("^N2", cli_surv_filtered$ajcc_pathologic_n), "N2", cli_surv_filtered$ajcc_pathologic_n)
      cli_surv_filtered$ajcc_pathologic_n <- ifelse(grepl("^N3", cli_surv_filtered$ajcc_pathologic_n), "N3", cli_surv_filtered$ajcc_pathologic_n)
      
      cli_surv_filtered$ajcc_pathologic_m <- ifelse(grepl("^M1", cli_surv_filtered$ajcc_pathologic_m), "M1", cli_surv_filtered$ajcc_pathologic_m)
      
      data_tc$vitalstatus = NULL
      #
      data_tc = data_tc[intersect(rownames(data_tc) , cli_surv_filtered$submitter_id),]
      cli_surv_filtered = cli_surv_filtered[which(cli_surv_filtered$submitter_id %in% intersect(rownames(data_tc) , cli_surv_filtered$submitter_id)),]
      
      row_names_data_tc <- rownames(data_tc)
      cli_surv_filtered_ordered <- cli_surv_filtered[match(row_names_data_tc, cli_surv_filtered$submitter_id), ]
      
      if (all.equal(rownames(data_tc),  cli_surv_filtered_ordered$submitter_id )) {

        data_tc$pat_t = cli_surv_filtered_ordered$ajcc_pathologic_t 
        data_tc$pat_n = cli_surv_filtered_ordered$ajcc_pathologic_n 
        data_tc$pat_m = cli_surv_filtered_ordered$ajcc_pathologic_m 
      }
      
      data_tc_tmp = data_tc[,which(!colnames(data_tc) %in% c("pat_t","pat_n","pat_m"))]
      data_tc_meta = data_tc[,c("pat_t","pat_n","pat_m")]
      
      # # Convert the matrix to a numeric matrix
      # data_tc_numeric <- matrix(as.numeric(unlist(data_tc_tmp)), nrow = nrow(data_tc_tmp))
      # 
      # # Convert the numeric matrix to a vector
      # vec <- as.vector(data_tc_numeric)
      # 
      # # Create the histogram
      # hist(vec)

      data_tc_tmp[data_tc_tmp > 15] = 15
      
      data_tc_filt = cbind(data_tc_tmp, data_tc_meta)
      ######
      
      data_tc_filt_df = data_tc_filt[order(data_tc_filt$pat_t),]
      annotation_df <- data.frame(pat_t = data_tc_filt_df$pat_t,
                                  pat_n = data_tc_filt_df$pat_n,
                                  pat_m = data_tc_filt_df$pat_m)
      rownames(annotation_df) <- rownames(data_tc_filt_df)
      
      # Create a named color vector for the unique values of vital_status
      
      num_pat_t = brewer.pal(length(unique(data_tc_filt_df$pat_t)), "RdGy")
      col_pat_t = setNames(num_pat_t, unique(data_tc_filt_df$pat_t))
      
      num_pat_n = brewer.pal(length(unique(data_tc_filt_df$pat_n)), "Spectral")
      col_pat_n = setNames(num_pat_n, unique(data_tc_filt_df$pat_n))
      
      num_pat_m = brewer.pal(length(unique(data_tc_filt_df$pat_m)), "Set3")
      col_pat_m = setNames(num_pat_m, unique(data_tc_filt_df$pat_m))
      
      tnm_colors = list(pat_t = col_pat_t , pat_n = col_pat_n, pat_m = col_pat_m)
      
      png(filename = paste0(CancerType,"_net_pathway_clusterT_tnm.png"),
          width = 25, height = 25,  units = "cm" ,pointsize = 12,
          bg = "white", res = 1200, family = "")
      
      tmp = pheatmap::pheatmap(as.matrix(t(data_tc_filt_df[,which(!colnames(data_tc_filt_df) %in% c("pat_t","pat_n","pat_m"))])),
                               annotation_col = annotation_df,
                               annotation_colors = tnm_colors,
                               cluster_cols = T)
      print(tmp)
      
      dev.off()
      
      png(filename = paste0(CancerType,"_net_pathway_clusterF_tnm.png"),
          width = 25, height = 25,  units = "cm" ,pointsize = 12,
          bg = "white", res = 1200, family = "")
      
      tmp2 = pheatmap::pheatmap(as.matrix(t(data_tc_filt_df[,which(!colnames(data_tc_filt_df) %in% c("pat_t","pat_n","pat_m"))])),
                                annotation_col = annotation_df,
                                annotation_colors = tnm_colors,
                                cluster_cols = F)
      print(tmp2)
      
      dev.off()
      
      png(filename = paste0(CancerType,"_net_pathway_t_complex.png"),
          width = 25, height = 25,  units = "cm" ,pointsize = 12,
          bg = "white", res = 1200, family = "")
      tmp3 = ComplexHeatmap::pheatmap(as.matrix(t(data_tc_filt_df[,which(!colnames(data_tc_filt_df) %in% c("pat_t","pat_n","pat_m"))])),
                                      column_split = annotation_df$pat_t,
                                      annotation_col = annotation_df,
                                      annotation_colors = tnm_colors,
                                      cluster_cols = T)
      print(tmp3)
      
      dev.off()
      
      png(filename = paste0(CancerType,"_net_pathway_n_complex.png"),
          width = 25, height = 25,  units = "cm" ,pointsize = 12,
          bg = "white", res = 1200, family = "")
      
      tmp4 = ComplexHeatmap::pheatmap(as.matrix(t(data_tc_filt_df[,which(!colnames(data_tc_filt_df) %in% c("pat_t","pat_n","pat_m"))])),
                                      column_split = annotation_df$pat_n,
                                      annotation_col = annotation_df,
                                      annotation_colors = tnm_colors,
                                      cluster_cols = T)
      
      print(tmp4)
      
      dev.off()
      
      png(filename = paste0(CancerType,"_net_pathway_m_complex.png"),
          width = 25, height = 25,  units = "cm" ,pointsize = 12,
          bg = "white", res = 1200, family = "")
      tmp5 = ComplexHeatmap::pheatmap(as.matrix(t(data_tc_filt_df[,which(!colnames(data_tc_filt_df) %in% c("pat_t","pat_n","pat_m"))])),
                                      column_split = annotation_df$pat_m,
                                      annotation_col = annotation_df,
                                      annotation_colors = tnm_colors,
                                      cluster_cols = T)
      
      print(tmp5)
      
      dev.off()
      
      remove(tmp,tmp2,tmp3,tmp4,tmp5,data_tc_filt_df,data_tc_filt,data_tc)
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
      
      if (CancerType == "TCGA-COADREAD") {
        cli_surv = cli[cli$project %in% c("TCGA-COAD","TCGA-READ"),
                       c("submitter_id",
                         "ajcc_pathologic_t",
                         "ajcc_pathologic_n",
                         "ajcc_pathologic_m")]
        
      } else if (CancerType == "TCGA-KIDNEY") {
        cli_surv = cli[cli$project %in% c("TCGA-KIRP","TCGA-KIRC","TCGA-KICH"),
                       c("submitter_id",
                         "ajcc_pathologic_t",
                         "ajcc_pathologic_n",
                         "ajcc_pathologic_m")]
        
      } else {
        cli_surv = cli[cli$project == CancerType,
                       c("submitter_id",
                         "ajcc_pathologic_t",
                         "ajcc_pathologic_n",
                         "ajcc_pathologic_m")]
      }
      
      if (sum(!is.na(cli_surv$ajcc_pathologic_t)) == 0 || sum(!is.na(cli_surv$ajcc_pathologic_n)) == 0 || sum(!is.na(cli_surv$ajcc_pathologic_m)) == 0) {
        print(paste0( "there are not tmn data in ",CancerType))
        next
      }
      
      cli_surv_filtered = cli_surv[which(!is.na(cli_surv$ajcc_pathologic_t)),]
      cli_surv_filtered = cli_surv_filtered[which(!is.na(cli_surv_filtered$ajcc_pathologic_n)),]
      cli_surv_filtered = cli_surv_filtered[which(!is.na(cli_surv_filtered$ajcc_pathologic_m)),]
      
      cli_surv_filtered$ajcc_pathologic_t <- ifelse(grepl("^T1", cli_surv_filtered$ajcc_pathologic_t), "T1", cli_surv_filtered$ajcc_pathologic_t)
      cli_surv_filtered$ajcc_pathologic_t <- ifelse(grepl("^T2", cli_surv_filtered$ajcc_pathologic_t), "T2", cli_surv_filtered$ajcc_pathologic_t)
      cli_surv_filtered$ajcc_pathologic_t <- ifelse(grepl("^T3", cli_surv_filtered$ajcc_pathologic_t), "T3", cli_surv_filtered$ajcc_pathologic_t)
      cli_surv_filtered$ajcc_pathologic_t <- ifelse(grepl("^T4", cli_surv_filtered$ajcc_pathologic_t), "T4", cli_surv_filtered$ajcc_pathologic_t)
      
      cli_surv_filtered$ajcc_pathologic_n <- ifelse(grepl("^N1", cli_surv_filtered$ajcc_pathologic_n), "N1", cli_surv_filtered$ajcc_pathologic_n)
      cli_surv_filtered$ajcc_pathologic_n <- ifelse(grepl("^N2", cli_surv_filtered$ajcc_pathologic_n), "N2", cli_surv_filtered$ajcc_pathologic_n)
      cli_surv_filtered$ajcc_pathologic_n <- ifelse(grepl("^N3", cli_surv_filtered$ajcc_pathologic_n), "N3", cli_surv_filtered$ajcc_pathologic_n)
      
      cli_surv_filtered$ajcc_pathologic_m <- ifelse(grepl("^M1", cli_surv_filtered$ajcc_pathologic_m), "M1", cli_surv_filtered$ajcc_pathologic_m)
      
      data_tc$vitalstatus = NULL
      #
      data_tc = data_tc[intersect(rownames(data_tc) , cli_surv_filtered$submitter_id),]
      cli_surv_filtered = cli_surv_filtered[which(cli_surv_filtered$submitter_id %in% intersect(rownames(data_tc) , cli_surv_filtered$submitter_id)),]
      
      row_names_data_tc <- rownames(data_tc)
      cli_surv_filtered_ordered <- cli_surv_filtered[match(row_names_data_tc, cli_surv_filtered$submitter_id), ]
      
      if (all.equal(rownames(data_tc),  cli_surv_filtered_ordered$submitter_id )) {
        
        data_tc$pat_t = cli_surv_filtered_ordered$ajcc_pathologic_t 
        data_tc$pat_n = cli_surv_filtered_ordered$ajcc_pathologic_n 
        data_tc$pat_m = cli_surv_filtered_ordered$ajcc_pathologic_m 
      }
      
      data_tc_tmp = data_tc[,which(!colnames(data_tc) %in% c("pat_t","pat_n","pat_m"))]
      data_tc_meta = data_tc[,c("pat_t","pat_n","pat_m")]
      
      # # Convert the matrix to a numeric matrix
      # data_tc_numeric <- matrix(as.numeric(unlist(data_tc_tmp)), nrow = nrow(data_tc_tmp))
      # 
      # # Convert the numeric matrix to a vector
      # vec <- as.vector(data_tc_numeric)
      # 
      # # Create the histogram
      # hist(vec)
      
      data_tc_tmp[data_tc_tmp > 15] = 15
      
      data_tc_filt = cbind(data_tc_tmp, data_tc_meta)
      ######
      
      data_tc_filt_df = data_tc_filt[order(data_tc_filt$pat_t),]
      annotation_df <- data.frame(pat_t = data_tc_filt_df$pat_t,
                                  pat_n = data_tc_filt_df$pat_n,
                                  pat_m = data_tc_filt_df$pat_m)
      rownames(annotation_df) <- rownames(data_tc_filt_df)
      
      # Create a named color vector for the unique values of vital_status
      
      num_pat_t = brewer.pal(length(unique(data_tc_filt_df$pat_t)), "RdGy")
      col_pat_t = setNames(num_pat_t, unique(data_tc_filt_df$pat_t))
      
      num_pat_n = brewer.pal(length(unique(data_tc_filt_df$pat_n)), "Spectral")
      col_pat_n = setNames(num_pat_n, unique(data_tc_filt_df$pat_n))
      
      num_pat_m = brewer.pal(length(unique(data_tc_filt_df$pat_m)), "Set3")
      col_pat_m = setNames(num_pat_m, unique(data_tc_filt_df$pat_m))
      
      tnm_colors = list(pat_t = col_pat_t , pat_n = col_pat_n, pat_m = col_pat_m)
      
      png(filename = paste0(CancerType,"_",type,"_pathway_clusterT_tnm.png"),
          width = 25, height = 25,  units = "cm" ,pointsize = 12,
          bg = "white", res = 1200, family = "")
      
      tmp = pheatmap::pheatmap(as.matrix(t(data_tc_filt_df[,which(!colnames(data_tc_filt_df) %in% c("pat_t","pat_n","pat_m"))])),
                               annotation_col = annotation_df,
                               annotation_colors = tnm_colors,
                               cluster_cols = T)
      print(tmp)
      
      dev.off()
      
      png(filename = paste0(CancerType,"_",type,"_pathway_clusterF_tnm.png"),
          width = 25, height = 25,  units = "cm" ,pointsize = 12,
          bg = "white", res = 1200, family = "")
      
      tmp2 = pheatmap::pheatmap(as.matrix(t(data_tc_filt_df[,which(!colnames(data_tc_filt_df) %in% c("pat_t","pat_n","pat_m"))])),
                                annotation_col = annotation_df,
                                annotation_colors = tnm_colors,
                                cluster_cols = F)
      print(tmp2)
      
      dev.off()
      
      png(filename = paste0(CancerType,"_",type,"_pathway_t_complex.png"),
          width = 25, height = 25,  units = "cm" ,pointsize = 12,
          bg = "white", res = 1200, family = "")
      tmp3 = ComplexHeatmap::pheatmap(as.matrix(t(data_tc_filt_df[,which(!colnames(data_tc_filt_df) %in% c("pat_t","pat_n","pat_m"))])),
                                      column_split = annotation_df$pat_t,
                                      annotation_col = annotation_df,
                                      annotation_colors = tnm_colors,
                                      cluster_cols = T)
      print(tmp3)
      
      dev.off()
      
      png(filename = paste0(CancerType,"_",type,"_pathway_n_complex.png"),
          width = 25, height = 25,  units = "cm" ,pointsize = 12,
          bg = "white", res = 1200, family = "")
      
      tmp4 = ComplexHeatmap::pheatmap(as.matrix(t(data_tc_filt_df[,which(!colnames(data_tc_filt_df) %in% c("pat_t","pat_n","pat_m"))])),
                                      column_split = annotation_df$pat_n,
                                      annotation_col = annotation_df,
                                      annotation_colors = tnm_colors,
                                      cluster_cols = T)
      
      print(tmp4)
      
      dev.off()
      
      png(filename = paste0(CancerType,"_",type,"_pathway_m_complex.png"),
          width = 25, height = 25,  units = "cm" ,pointsize = 12,
          bg = "white", res = 1200, family = "")
      tmp5 = ComplexHeatmap::pheatmap(as.matrix(t(data_tc_filt_df[,which(!colnames(data_tc_filt_df) %in% c("pat_t","pat_n","pat_m"))])),
                                      column_split = annotation_df$pat_m,
                                      annotation_col = annotation_df,
                                      annotation_colors = tnm_colors,
                                      cluster_cols = T)
      
      print(tmp5)
      
      dev.off()
      
      remove(tmp,tmp2,tmp3,tmp4,tmp5,data_tc_filt_df,data_tc_filt,data_tc)
      }
  }
}


