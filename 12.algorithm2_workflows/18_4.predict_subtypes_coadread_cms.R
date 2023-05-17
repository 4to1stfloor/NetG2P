library(readxl)
library(devtools)
library(CMSclassifier)
library(org.Hs.eg.db)
library(data.table)
library(scater)
library(scran)
library(pheatmap)
library(dplyr)
library(stringr)
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

cms_genes = read_xlsx("~/nas/00.data/filtered_TCGA/34.TCGA-COADREAD/CMS_genes_merge.xlsx")
cms = readRDS("~/nas/00.data/filtered_TCGA/34.TCGA-COADREAD/TCGA-COADREAD_CMS.rds")

num_CancerType = "34.TCGA-COADREAD"

# for all

main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))

# call input

sce = readRDS(paste0(sce_path,num_CancerType,"/", CancerType,".sce_raw_snv_cnv.rds"))
mut = readRDS(paste0(sce_path,num_CancerType,"/", CancerType,"_maf.rds"))
net = readRDS(paste0(main.path_tc, "/net_prop_total_587.rds"))

cms_genes_sym = cms_genes$SYMBOL

cms_genes_sym = cms_genes_sym[!cms_genes_sym %in% c("WARS","MLL2","MLL3", "MLL","HLAA", "HLAB", "AKD1", "MLLT4", "WHSC1" , "ERBB2IP", "PARK2", "WHSC1L1")]
cms_genes_sym = c(cms_genes_sym , "WARS1", "KMT2D", "KMT2C", "KMT2A", "HLA-A", "HLA-B", "AK9", "AFDN", "NSD2", "ERBIN", "PRKN", "NSD3")

sce_exp = tcga.calc.zscore(sce = sce, cms_genes_sym)
sce_exp = as.data.frame(sce_exp)
sce_exp[sce_exp > 5] <- 5
sce_exp[sce_exp < -5] <- -5

# Convert the matrix to a numeric matrix
sce_exp_numeric <- matrix(as.numeric(unlist(sce_exp)), nrow = nrow(sce_exp))

# Convert the numeric matrix to a vector
vec <- as.vector(sce_exp_numeric)

# Create the histogram
hist(vec)

sce_exp = sce_exp[complete.cases(sce_exp),]

sce_exp_df = as.data.frame(t(sce_exp))
if (all.equal(rownames(sce_exp_df) , rownames(cms)) ) {
  sce_exp_df$CMS = cms$RF.details.RF.nearestCMS
}

sce_exp_df = sce_exp_df[order(sce_exp_df$CMS),]

annotation_df <- data.frame(CMS = sce_exp_df$CMS)
rownames(annotation_df) <- rownames(sce_exp_df)

# Create a named color vector for the unique values of vital_status
CMS_colors <- c("CMS1" = "yellow", "CMS2" = "black", "CMS3" = "green", "CMS4" = "red")
names(CMS_colors) <- unique(annotation_df$CMS)


png(filename = paste0(CancerType,"_exp_clusterT_CMS.png"),
    width = 35, height = 35,  units = "cm" ,pointsize = 12,
    bg = "white", res = 1200, family = "")

tmp = pheatmap::pheatmap(as.matrix(t(sce_exp_df[,-ncol(sce_exp_df)])),
                         annotation_col = annotation_df,
                         annotation_colors = list(CMS = CMS_colors),
                         cluster_cols = T)
print(tmp)
dev.off()

png(filename = paste0(CancerType,"_exp_column_clusterF_CMS.png"),
    width = 35, height = 35,  units = "cm" ,pointsize = 12,
    bg = "white", res = 1200, family = "")

tmp2 = pheatmap::pheatmap(as.matrix(t(sce_exp_df[,-ncol(sce_exp_df)])),
                          annotation_col = annotation_df,
                          annotation_colors = list(CMS = CMS_colors),
                          cluster_cols = F)
print(tmp2)

dev.off()

png(filename = paste0(CancerType,"_exp_cluster_complex_CMS.png"),
    width = 35, height = 35,  units = "cm" ,pointsize = 12,
    bg = "white", res = 1200, family = "")
tmp3 = ComplexHeatmap::pheatmap(as.matrix(t(sce_exp_df[,-ncol(sce_exp_df)])),
                                column_split = annotation_df$CMS,
                                annotation_col = annotation_df,
                                annotation_colors = list(CMS = CMS_colors),
                                cluster_cols = T)

print(tmp3)
dev.off()

## mut
mut_filtered <- subsetMaf(maf = mut, query = "Variant_Classification != 'Nonsense_Mutation'")

mut_count = mutCountMatrix(mut_filtered)

mut_count_filtered = mut_count[cms_genes_sym,]

mut_count_filtered_tdf = as.data.frame(t(mut_count_filtered))

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

rownames(mut_count_filtered_tdf) = substr(rownames(mut_count_filtered_tdf),1,12) 

cms$submitter_id = substr(rownames(cms), 1, 16)
cms$forfilt =  rownames(cms)
cms_01A = cms %>% filter(str_detect(forfilt, "-01A-")) 
cms_01A$forfilt = NULL

dup_cms_df = cms_01A[which(duplicated(substr(rownames(cms_01A),1,12))),]
cms_filtered = cms_01A

for (dup_pat in rownames(dup_cms_df)) {
  tmp_mut = grep(paste0(substr(dup_pat,1,12),"-*"), rownames(cms_01A), value = TRUE)
  if (cms_filtered[tmp_mut,][1,]$RF.details.RF.nearestCMS == cms_filtered[tmp_mut,][2,]$RF.details.RF.nearestCMS) {
    cms_filtered <- subset(cms_filtered, rownames(cms_filtered) != tmp_mut[2])
  } else {
    print(tmp_mut)
  }
}


rownames(cms_filtered) = substr(rownames(cms_filtered),1,12)

common_pat = intersect(rownames(mut_count_filtered_tdf) , rownames(cms_filtered))
mut_count_filtered_tdf_common = mut_count_filtered_tdf[common_pat,]
cms_filtered_common = cms_filtered[common_pat,]
mut_count_filtered_tdf_common[mut_count_filtered_tdf_common >5] = 5 
if (all.equal(rownames(mut_count_filtered_tdf_common), rownames(cms_filtered_common))) {
  mut_count_filtered_tdf_common$CMS = cms_filtered_common$RF.details.RF.nearestCMS
}

mut_count_filtered_tdf_common = mut_count_filtered_tdf_common[order(mut_count_filtered_tdf_common$CMS),]

annotation_df <- data.frame(CMS = mut_count_filtered_tdf_common$CMS)
rownames(annotation_df) <- rownames(mut_count_filtered_tdf_common)

# Create a named color vector for the unique values of vital_status
CMS_colors <- c("CMS1" = "yellow", "CMS2" = "black", "CMS3" = "green", "CMS4" = "red")
names(CMS_colors) <- unique(annotation_df$CMS)


png(filename = paste0(CancerType,"_mut_clusterT_CMS_rowgene_col_pat.png"),
    width = 25, height = 25,  units = "cm" ,pointsize = 12,
    bg = "white", res = 1200, family = "")

tmp = pheatmap::pheatmap(as.matrix(t(mut_count_filtered_tdf_common[,-ncol(mut_count_filtered_tdf_common)])),
                         annotation_col = annotation_df,
                         annotation_colors = list(CMS = CMS_colors),
                         cluster_cols = T,show_rownames = F,show_colnames = F)
print(tmp)

dev.off()

png(filename = paste0(CancerType,"_mut_clusterF_CMS_rowgene_col_pat.png"),
    width = 25, height = 25,  units = "cm" ,pointsize = 12,
    bg = "white", res = 1200, family = "")

tmp2 = pheatmap::pheatmap(as.matrix(t(mut_count_filtered_tdf_common[,-ncol(mut_count_filtered_tdf_common)])),
                          annotation_col = annotation_df,
                          annotation_colors = list(CMS = CMS_colors),
                          cluster_cols = F,show_rownames = F,show_colnames = F)
print(tmp2)

dev.off()

png(filename = paste0(CancerType,"_mut_cluster_complex_CMS.png"),
    width = 25, height = 25,  units = "cm" ,pointsize = 12,
    bg = "white", res = 1200, family = "")
tmp3 = ComplexHeatmap::pheatmap(as.matrix(t(mut_count_filtered_tdf_common[,-ncol(mut_count_filtered_tdf_common)])),
                                column_split = annotation_df$CMS,
                                annotation_col = annotation_df,
                                annotation_colors = list(CMS = CMS_colors),
                                cluster_cols = T)

print(tmp3)
dev.off()

mut_filtered@clinical.data$CMS = NA
mut_filtered@clinical.data$CMS = as.character(mut_filtered@clinical.data$CMS)

mut_filtered@clinical.data$submitter_id = substr(mut_filtered@clinical.data$Tumor_Sample_Barcode,1,12)

for (maf_patients in mut_filtered@clinical.data$submitter_id ) {
  if (maf_patients  %in% rownames(cms_filtered_common)) {
    mut_filtered@clinical.data[which(mut_filtered@clinical.data$submitter_id == maf_patients),]$CMS = cms_filtered_common[which(rownames(cms_filtered_common) == maf_patients),]$RF.details.RF.nearestCMS
  } 
  
} 
mut_filtered@clinical.data = mut_filtered@clinical.data[which(!is.na(mut_filtered@clinical.data$CMS)),]

png(filename = paste0(CancerType,"_mut_oncoplot_CMS.png"),
    width = 30, height = 30,  units = "cm" ,pointsize = 12,
    bg = "white", res = 1200, family = "")


tmp4 = oncoplot(maf = mut_filtered,
                genes = cms_genes_sym,
                clinicalFeatures = "CMS",
                sortByAnnotation = TRUE)

print(tmp4)
dev.off()

## net
net_cms = net[which(rownames(net) %in% cms_genes_sym),]

net_cms_tdf = as.data.frame(t(net_cms))

dup_net_df = net_cms_tdf[which(duplicated(substr(rownames(net_cms_tdf),1,12))),]

for (dup_pat in rownames(dup_net_df)) {
  tmp = grep(paste0(substr(dup_pat,1,12),"-*"), rownames(net_cms_tdf), value = TRUE)
  
  mean_df <- data.frame(lapply(seq_along(net_cms_tdf[tmp[1],]), function(i) {
    round((net_cms_tdf[tmp[1],][[i]] + net_cms_tdf[tmp[2],][[i]]) / 2)
  }))
  colnames(mean_df) <- colnames(net_cms_tdf[tmp[1],])
  rownames(mean_df) = substr(dup_pat,1,12)
  # remove two row and replace new mean row
  
  net_cms_tdf <- subset(net_cms_tdf, rownames(net_cms_tdf) != tmp[1])
  net_cms_tdf <- subset(net_cms_tdf, rownames(net_cms_tdf) != tmp[2])
  
  net_cms_tdf = rbind(net_cms_tdf, mean_df)
  
}

rownames(net_cms_tdf) = substr(rownames(net_cms_tdf),1,12) 

common_net_pat = intersect(rownames(net_cms_tdf) , rownames(cms_filtered))
net_cms_tdf_common = net_cms_tdf[common_net_pat,]
cms_net_filtered_common = cms_filtered[common_net_pat,]

if (all.equal(rownames(net_cms_tdf_common), rownames(cms_net_filtered_common))) {
  net_cms_tdf_common$CMS = cms_net_filtered_common$RF.details.RF.nearestCMS
}

net_cms_tdf_common = net_cms_tdf_common[order(net_cms_tdf_common$CMS),]

annotation_df <- data.frame(CMS = net_cms_tdf_common$CMS)
rownames(annotation_df) <- rownames(net_cms_tdf_common)

# Create a named color vector for the unique values of vital_status
CMS_colors <- c("CMS1" = "yellow", "CMS2" = "black", "CMS3" = "green", "CMS4" = "red")
names(CMS_colors) <- unique(annotation_df$CMS)


png(filename = paste0(CancerType,"_net_clusterT_CMS.png"),
    width = 35, height = 35,  units = "cm" ,pointsize = 12,
    bg = "white", res = 1200, family = "")

tmp = pheatmap::pheatmap(as.matrix(t(net_cms_tdf_common[,-ncol(net_cms_tdf_common)])),
                         annotation_col = annotation_df,
                         annotation_colors = list(CMS = CMS_colors),
                         cluster_cols = T)
print(tmp)

dev.off()

png(filename = paste0(CancerType,"_net_clusterF_CMS.png"),
    width = 35, height = 35,  units = "cm" ,pointsize = 12,
    bg = "white", res = 1200, family = "")

tmp2 = pheatmap::pheatmap(as.matrix(t(net_cms_tdf_common[,-ncol(net_cms_tdf_common)])),
                          annotation_col = annotation_df,
                          annotation_colors = list(CMS = CMS_colors),
                          cluster_cols = F)
print(tmp2)

dev.off()

png(filename = paste0(CancerType,"_net_cluster_complex_CMS.png"),
    width = 35, height = 35,  units = "cm" ,pointsize = 12,
    bg = "white", res = 1200, family = "")
tmp3 = ComplexHeatmap::pheatmap(as.matrix(t(net_cms_tdf_common[,-ncol(net_cms_tdf_common)])),
                                column_split = annotation_df$CMS,
                                annotation_col = annotation_df,
                                annotation_colors = list(CMS = CMS_colors),
                                cluster_cols = T)

print(tmp3)
dev.off()

#
png(filename = paste0(CancerType,"_net_clusterT_CMS_rowgenes_colpat.png"),
    width = 25, height = 25,  units = "cm" ,pointsize = 12,
    bg = "white", res = 1200, family = "")

tmp = pheatmap::pheatmap(as.matrix(t(net_cms_tdf_common[,-ncol(net_cms_tdf_common)])),
                         annotation_col = annotation_df,
                         annotation_colors = list(CMS = CMS_colors),
                         cluster_cols = T,show_rownames = F,show_colnames = F)
print(tmp)

dev.off()

png(filename = paste0(CancerType,"_net_clusterF_CMS_rowgenes_colpat.png"),
    width = 35, height = 35,  units = "cm" ,pointsize = 12,
    bg = "white", res = 1200, family = "")

tmp2 = pheatmap::pheatmap(as.matrix(t(net_cms_tdf_common[,-ncol(net_cms_tdf_common)])),
                          annotation_col = annotation_df,
                          annotation_colors = list(CMS = CMS_colors),
                          cluster_cols = F,show_rownames = F,show_colnames = F)
print(tmp2)

dev.off()

png(filename = paste0(CancerType,"_net_cluster_complex_CMS_rowgenes_colpat.png"),
    width = 35, height = 35,  units = "cm" ,pointsize = 12,
    bg = "white", res = 1200, family = "")
tmp3 = ComplexHeatmap::pheatmap(as.matrix(t(net_cms_tdf_common[,-ncol(net_cms_tdf_common)])),
                                column_split = annotation_df$CMS,
                                annotation_col = annotation_df,
                                annotation_colors = list(CMS = CMS_colors),
                                cluster_cols = T,show_rownames = F,show_colnames = F)

print(tmp3)
dev.off()
