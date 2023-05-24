library(readxl)
library(devtools)
library(CMSclassifier)
library(org.Hs.eg.db)
library(data.table)
library(scater)
library(scran)
library(pheatmap)
require(TCGAbiolinks)
require(SummarizedExperiment)
require(SingleCellExperiment)
library(dplyr)
library(RColorBrewer)
setwd("~/nas/04.Results/subtypes/")

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
sce_path = "/mnt/gluster_server/data/raw/TCGA_data/00.data/"

cli = fread("~/nas/99.reference/all_clin_indexed.csv")

num_CancerType = "30.TCGA-BRCA"

# for all

main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))

# call input

sce = readRDS(paste0(sce_path,num_CancerType,"/", CancerType,".sce_raw_snv_cnv.rds"))
mut = readRDS(paste0(sce_path,num_CancerType,"/", CancerType,"_maf.rds"))
net = readRDS(paste0(main.path_tc, "/net_prop_total_966.rds"))
short_long = readRDS(paste0(filepath,"04.Results/short_long/",CancerType,"_best_features_short_long.rds"))

# call pam50 genes
library(genefu)
data("pam50")
pam50genes = rownames(pam50$centroids)

# genes name are changed 

pam50genes = pam50genes[!pam50genes %in% c("CDCA1","KNTC2","ORC6L")]
pam50genes = c(pam50genes , "BRCA1", "BRCA2", "NUF2", "NDC80", "ORC6")

sce_exp = tcga.calc.zscore(sce = sce, pam50genes)
sce_exp = as.data.frame(sce_exp)

# # Convert the matrix to a numeric matrix
# sce_exp_numeric <- matrix(as.numeric(unlist(sce_exp)), nrow = nrow(sce_exp))
# 
# # Convert the numeric matrix to a vector
# vec <- as.vector(sce_exp_numeric)
# 
# # Create the histogram
# hist(vec)

sce_exp[sce_exp > 5] <- 5
sce_exp[sce_exp < -5] <- -5

sce_exp_filt = sce_exp[,grep("*-01A", colnames(sce_exp))]

for (dup_pat in colnames(sce_exp_filt)[grep("*-01A.1", colnames(sce_exp_filt))]) {
  dup_pat_filt = gsub("\\.1","", dup_pat)
  first_tmp = as.data.frame(sce_exp_filt[,dup_pat_filt])
  second_tmp = as.data.frame(sce_exp_filt[,dup_pat])
  
  tmp_mean = as.data.frame(apply(X = cbind(first_tmp,second_tmp), 1, mean))
  rownames(tmp_mean) = rownames(sce_exp_filt)
  colnames(tmp_mean) = dup_pat_filt
  sce_exp_filt[,dup_pat] = NULL
  sce_exp_filt[,dup_pat_filt] = NULL
  
  sce_exp_filt = cbind(sce_exp_filt,tmp_mean)
  
}

sce_exp_filt_df = as.data.frame(t(sce_exp_filt))

sce_data = as.data.frame(colData(sce))
sce_data_filt = sce_data[rownames(sce_exp_filt_df),]

if (all.equal(rownames(sce_exp_filt_df) , rownames(sce_data_filt)) ) {
  sce_exp_filt_df$pam50 = sce_data_filt$paper_BRCA_Subtype_PAM50
}
sce_exp_filt_df_wona = sce_exp_filt_df[!is.na(sce_exp_filt_df$pam50),]
rownames(sce_exp_filt_df_wona) = substr(rownames(sce_exp_filt_df_wona), 1,12)
sce_exp_filt_df_wona = sce_exp_filt_df_wona[intersect(rownames(sce_exp_filt_df_wona), rownames(short_long)),]
short_long_filt = short_long[intersect(rownames(sce_exp_filt_df_wona), rownames(short_long)),]

if (all.equal(rownames(sce_exp_filt_df_wona), rownames(short_long_filt))) {
  sce_exp_filt_df_wona$cluster = short_long_filt$cluster
}

sce_exp_filt_df_wona = sce_exp_filt_df_wona[order(sce_exp_filt_df_wona$pam50),]

annotation_df <- data.frame(pam50 = sce_exp_filt_df_wona$pam50,
                            cluster = sce_exp_filt_df_wona$cluster)
rownames(annotation_df) <- rownames(sce_exp_filt_df_wona)

# Create a named color vector for the unique values of vital_status
num_pam50 = brewer.pal(length(unique(sce_exp_filt_df_wona$pam50)), "Spectral")
col_pam50 = setNames(num_pam50, unique(sce_exp_filt_df_wona$pam50))

cluster_colors <- c("short" = "yellow", "long" = "black")
subtypes_colors = list(pam50 = col_pam50 , cluster = cluster_colors)

png(filename = paste0(CancerType,"_exp_w_short_long_clusterT_pam50.png"),
    width = 25, height = 25,  units = "cm" ,pointsize = 12,
    bg = "white", res = 1200, family = "")

tmp = pheatmap::pheatmap(as.matrix(t(sce_exp_filt_df_wona[,which(!colnames(sce_exp_filt_df_wona) %in% c("pam50","cluster"))])),
                         annotation_col = annotation_df,
                         annotation_colors = subtypes_colors,
                         cluster_cols = T)
print(tmp)

dev.off()

png(filename = paste0(CancerType,"_exp_w_short_long_clusterF_pam50.png"),
    width = 25, height = 25,  units = "cm" ,pointsize = 12,
    bg = "white", res = 1200, family = "")

tmp2 = pheatmap::pheatmap(as.matrix(t(sce_exp_filt_df_wona[,which(!colnames(sce_exp_filt_df_wona) %in% c("pam50","cluster"))])),
                          annotation_col = annotation_df,
                          annotation_colors = subtypes_colors,
                          cluster_cols = F)
print(tmp2)

dev.off()

png(filename = paste0(CancerType,"_exp_w_short_long_cluster_complex_pam50.png"),
    width = 25, height = 25,  units = "cm" ,pointsize = 12,
    bg = "white", res = 1200, family = "")
tmp3 = ComplexHeatmap::pheatmap(as.matrix(t(sce_exp_filt_df_wona[,which(!colnames(sce_exp_filt_df_wona) %in% c("pam50","cluster"))])),
                                column_split = annotation_df$pam50,
                                annotation_col = annotation_df,
                                annotation_colors = subtypes_colors,
                                cluster_cols = T)

print(tmp3)
dev.off()

## mut
mut_filtered <- subsetMaf(maf = mut, query = "Variant_Classification != 'Nonsense_Mutation'")
mut_filtered@clinical.data$pam50 = NA
mut_filtered@clinical.data$pam50= as.character(mut_filtered@clinical.data$pam50)
mut_filtered@clinical.data$cluster = NA
mut_filtered@clinical.data$cluster= as.character(mut_filtered@clinical.data$cluster)

mut_filtered@clinical.data$submitter_id = substr(mut_filtered@clinical.data$Tumor_Sample_Barcode,1,12)

for (maf_patients in substr(mut_filtered@clinical.data$Tumor_Sample_Barcode,1,12)) {
  if (maf_patients  %in% rownames(sce_exp_filt_df_wona)) {
    mut_filtered@clinical.data[which(mut_filtered@clinical.data$submitter_id == maf_patients),]$pam50 = sce_exp_filt_df_wona[which(rownames(sce_exp_filt_df_wona) == maf_patients),]$pam50
    mut_filtered@clinical.data[which(mut_filtered@clinical.data$submitter_id == maf_patients),]$cluster = sce_exp_filt_df_wona[which(rownames(sce_exp_filt_df_wona) == maf_patients),]$cluster
  } 
  
}

png(filename = paste0(CancerType,"_mut_w_short_long_onco_pam50.png"),
    width = 25, height = 25,  units = "cm" ,pointsize = 12,
    bg = "white", res = 1200, family = "")

tmp_onco = oncoplot(maf = mut_filtered,
                    genes = pam50genes,
                    clinicalFeatures = c("pam50","cluster"),
                    sortByAnnotation = TRUE
                    )
print(tmp_onco)
dev.off()

#
mut_count = mutCountMatrix(mut_filtered)

mut_count_filtered = mut_count[which(rownames(mut_count) %in% pam50genes),]

mut_count_filtered_pat = mut_count_filtered[,grep("*-01A", colnames(mut_count_filtered))]
mut_count_filtered_tdf = as.data.frame(t(mut_count_filtered_pat))

dup_df = mut_count_filtered_tdf[which(duplicated(substr(rownames(mut_count_filtered_tdf),1,12))),]

for (dup_pat in rownames(dup_df)) {
  tmp = grep(paste0(substr(dup_pat,1,12),"-*"), rownames(mut_count_filtered_tdf), value = TRUE)
  
  mean_df <- data.frame(lapply(seq_along(mut_count_filtered_tdf[tmp[1],]), function(i) {
    round((mut_count_filtered_tdf[tmp[1],][[i]] + mut_count_filtered_tdf[tmp[2],][[i]]) / 2)
  }))
  colnames(mean_df) <- colnames(mut_count_filtered_tdf[tmp[1],])
  rownames(mean_df) = substr(dup_pat,1,12)
  # remove two row and replace new mean row
  
  mut_count_filtered_tdf <- subset(mut_count_filtered_tdf, rownames(mut_count_filtered_tdf) != tmp[1])
  mut_count_filtered_tdf <- subset(mut_count_filtered_tdf, rownames(mut_count_filtered_tdf) != tmp[2])
  
  mut_count_filtered_tdf = rbind(mut_count_filtered_tdf, mean_df)
  
}

rownames(mut_count_filtered_tdf) = substr(rownames(mut_count_filtered_tdf),1,16) 
common_pat = intersect( rownames(sce_data_filt) , rownames(mut_count_filtered_tdf))
mut_count_filtered_tdf = mut_count_filtered_tdf[common_pat,]
sce_data_filt_w = sce_data_filt[common_pat,]

if (all.equal(rownames(sce_data_filt_w) , rownames(mut_count_filtered_tdf))) {
  mut_count_filtered_tdf$pam50 = sce_data_filt_w$paper_BRCA_Subtype_PAM50
  
}
mut_count_filtered_tdf = mut_count_filtered_tdf[which(!is.na(mut_count_filtered_tdf$pam50)),]

mut_count_filtered_tdf = mut_count_filtered_tdf[order(mut_count_filtered_tdf$pam50),]

# # Convert the matrix to a numeric matrix
# mut_pam50_numeric <- matrix(as.numeric(unlist(mut_pam50_wo_pam50)), nrow = nrow(mut_pam50_wo_pam50))
# 
# # Convert the numeric matrix to a vector
# vec <- as.vector(mut_pam50_numeric)
# 
# # Create the histogram
# hist(vec)

mut_tmp_pam50 = mut_count_filtered_tdf[,c("pam50")]
mut_pam50_wo_pam50 = mut_count_filtered_tdf[,which(!colnames(mut_count_filtered_tdf) %in% "pam50")]
mut_pam50_wo_pam50[mut_pam50_wo_pam50 >= 1] <- 1
mut_count_filtered_tdf = cbind(mut_pam50_wo_pam50 , pam50 = mut_tmp_pam50)
rownames(mut_count_filtered_tdf) = substr(rownames(mut_count_filtered_tdf), 1,12)
mut_count_filtered_tdf$cluster = NA

for (maf_patients in rownames(mut_count_filtered_tdf)) {
  if (maf_patients  %in% rownames(sce_exp_filt_df_wona)) {
    mut_count_filtered_tdf[which(rownames(mut_count_filtered_tdf) == maf_patients),]$cluster = sce_exp_filt_df_wona[which(rownames(sce_exp_filt_df_wona) == maf_patients),]$cluster
  } 
  
}
mut_count_filtered_tdf = mut_count_filtered_tdf[which(!is.na(mut_count_filtered_tdf$cluster)),]

annotation_df <- data.frame(pam50 = mut_count_filtered_tdf$pam50,
                            cluster = mut_count_filtered_tdf$cluster)
rownames(annotation_df) <- rownames(mut_count_filtered_tdf)

# Create a named color vector for the unique values of vital_status
num_pam50 = brewer.pal(length(unique(mut_count_filtered_tdf$pam50)), "Spectral")
col_pam50 = setNames(num_pam50, unique(mut_count_filtered_tdf$pam50))

cluster_colors <- c("short" = "yellow", "long" = "black")
subtypes_colors = list(pam50 = col_pam50 , cluster = cluster_colors)

png(filename = paste0(CancerType,"_mut_w_short_long_pam50.png"),
    width = 25, height = 25,  units = "cm" ,pointsize = 12,
    bg = "white", res = 1200, family = "")

tmp = pheatmap::pheatmap(as.matrix(t(mut_count_filtered_tdf[,which(!colnames(mut_count_filtered_tdf) %in% c("pam50","cluster"))])),
                         annotation_col = annotation_df,
                         annotation_colors = subtypes_colors,
                         cluster_cols = T)
print(tmp)

dev.off()

png(filename = paste0(CancerType,"_mut_w_short_long_clusterF_pam50.png"),
    width = 25, height = 25,  units = "cm" ,pointsize = 12,
    bg = "white", res = 1200, family = "")

tmp2 = pheatmap::pheatmap(as.matrix(t(mut_count_filtered_tdf[,which(!colnames(mut_count_filtered_tdf) %in% c("pam50","cluster"))])),
                          annotation_col = annotation_df,
                          annotation_colors = subtypes_colors,
                          cluster_cols = F)
print(tmp2)

dev.off()

png(filename = paste0(CancerType,"_mut_w_short_long_cluster_complex_pam50.png"),
    width = 25, height = 25,  units = "cm" ,pointsize = 12,
    bg = "white", res = 1200, family = "")
tmp3 = ComplexHeatmap::pheatmap(as.matrix(t(mut_count_filtered_tdf[,which(!colnames(mut_count_filtered_tdf) %in% c("pam50","cluster"))])),
                                column_split = annotation_df$pam50,
                                annotation_col = annotation_df,
                                annotation_colors = subtypes_colors,
                                cluster_cols = T)

print(tmp3)
dev.off()

# net

net_pam50 = net[which(rownames(net) %in% pam50genes),]

net_pam50_filtered_pat = net_pam50[,grep("*-01A", colnames(net_pam50))]
net_pam50_filtered_tdf = as.data.frame(t(net_pam50_filtered_pat))

dup_df = net_pam50_filtered_tdf[which(duplicated(substr(rownames(net_pam50_filtered_tdf),1,12))),]

if (nrow(dup_df) != 0) {
  
  for (dup_pat in rownames(dup_df)) {
    tmp = grep(paste0(substr(dup_pat,1,12),"-*"), rownames(net_pam50_filtered_tdf), value = TRUE)
    
    mean_df <- data.frame(lapply(seq_along(net_pam50_filtered_tdf[tmp[1],]), function(i) {
      round((net_pam50_filtered_tdf[tmp[1],][[i]] + net_pam50_filtered_tdf[tmp[2],][[i]]) / 2)
    }))
    colnames(mean_df) <- colnames(net_pam50_filtered_tdf[tmp[1],])
    rownames(mean_df) = substr(dup_pat,1,12)
    # remove two row and replace new mean row
    
    net_pam50_filtered_tdf <- subset(net_pam50_filtered_tdf, rownames(net_pam50_filtered_tdf) != tmp[1])
    net_pam50_filtered_tdf <- subset(net_pam50_filtered_tdf, rownames(net_pam50_filtered_tdf) != tmp[2])
    
    net_pam50_filtered_tdf = rbind(net_pam50_filtered_tdf, mean_df)
    
  }
  
}
rownames(net_pam50_filtered_tdf) = substr(rownames(net_pam50_filtered_tdf),1,16) 

common_pat = intersect( rownames(sce_data_filt) , rownames(net_pam50_filtered_tdf))
net_pam50_filtered_tdf = net_pam50_filtered_tdf[common_pat,]
sce_data_filt_w = sce_data_filt[common_pat,]

if (all.equal(rownames(sce_data_filt_w) , rownames(net_pam50_filtered_tdf))) {
  net_pam50_filtered_tdf$pam50 = sce_data_filt_w$paper_BRCA_Subtype_PAM50
  
}
net_pam50_filtered_tdf = net_pam50_filtered_tdf[which(!is.na(net_pam50_filtered_tdf$pam50)),]
net_pam50_filtered_tdf = net_pam50_filtered_tdf[order(net_pam50_filtered_tdf$pam50),]


# Convert the matrix to a numeric matrix
net_pam50_numeric <- matrix(as.numeric(unlist(net_pam50_wo_pam50)), nrow = nrow(net_pam50_wo_pam50))

# Convert the numeric matrix to a vector
vec <- as.vector(net_pam50_numeric)

# Create the histogram
hist(vec)

tmp_pam50 = net_pam50_filtered_tdf[,c("pam50")]
net_pam50_wo_pam50 = net_pam50_filtered_tdf[,which(!colnames(net_pam50_filtered_tdf) %in% "pam50")]
net_pam50_wo_pam50[net_pam50_wo_pam50 > 0.005] <- 0.005
net_pam50_filtered_tdf = cbind(net_pam50_wo_pam50 , pam50 = tmp_pam50)
rownames(net_pam50_filtered_tdf) = substr(rownames(net_pam50_filtered_tdf),1,12)
net_pam50_filtered_tdf$cluster = NA
for (net_patients in rownames(net_pam50_filtered_tdf)) {
  if (net_patients  %in% rownames(sce_exp_filt_df_wona)) {
    net_pam50_filtered_tdf[which(rownames(net_pam50_filtered_tdf) == net_patients),]$cluster = sce_exp_filt_df_wona[which(rownames(sce_exp_filt_df_wona) == net_patients),]$cluster
  } 
  
}

annotation_df <- data.frame(pam50 = net_pam50_filtered_tdf$pam50,
                            cluster = net_pam50_filtered_tdf$cluster)
rownames(annotation_df) <- rownames(net_pam50_filtered_tdf)

# Create a named color vector for the unique values of vital_status
num_pam50 = brewer.pal(length(unique(net_pam50_filtered_tdf$pam50)), "Spectral")
col_pam50 = setNames(num_pam50, unique(net_pam50_filtered_tdf$pam50))

cluster_colors <- c("short" = "yellow", "long" = "black")
subtypes_colors = list(pam50 = col_pam50 , cluster = cluster_colors)

png(filename = paste0(CancerType,"_net_w_short_long_clusterT_pam50.png"),
    width = 25, height = 25,  units = "cm" ,pointsize = 12,
    bg = "white", res = 1200, family = "")

tmp = pheatmap::pheatmap(as.matrix(t(net_pam50_filtered_tdf[,which(!colnames(net_pam50_filtered_tdf) %in% c("pam50","cluster"))])),
                         annotation_col = annotation_df,
                         annotation_colors = subtypes_colors,
                         cluster_cols = T)
print(tmp)

dev.off()

png(filename = paste0(CancerType,"_net_w_short_long_clusterF_pam50.png"),
    width = 25, height = 25,  units = "cm" ,pointsize = 12,
    bg = "white", res = 1200, family = "")

tmp2 = pheatmap::pheatmap(as.matrix(t(net_pam50_filtered_tdf[,which(!colnames(net_pam50_filtered_tdf) %in% c("pam50","cluster"))])),
                          annotation_col = annotation_df,
                          annotation_colors = subtypes_colors,
                          cluster_cols = F)
print(tmp2)

dev.off()

png(filename = paste0(CancerType,"_net_w_short_long_cluster_complex_pam50.png"),
    width = 25, height = 25,  units = "cm" ,pointsize = 12,
    bg = "white", res = 1200, family = "")
tmp3 = ComplexHeatmap::pheatmap(as.matrix(t(net_pam50_filtered_tdf[,which(!colnames(net_pam50_filtered_tdf) %in% c("pam50","cluster"))])),
                                column_split = annotation_df$pam50,
                                annotation_col = annotation_df,
                                annotation_colors = subtypes_colors,
                                cluster_cols = T)

print(tmp3)
dev.off()
