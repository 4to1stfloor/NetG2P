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

# call pam50 genes
library(genefu)
data("pam50")
pam50genes = rownames(pam50$centroids)

pam50genes = pam50genes[!pam50genes %in% c("CDCA1","KNTC2","ORC6L")]
pam50genes = c(pam50genes , "BRCA1", "BRCA2", "NUF2", "NDC80", "ORC6")

sce_exp = tcga.calc.zscore(sce = sce, pam50genes)
sce_exp = as.data.frame(sce_exp)

# Convert the matrix to a numeric matrix
sce_exp_numeric <- matrix(as.numeric(unlist(sce_exp)), nrow = nrow(sce_exp))

# Convert the numeric matrix to a vector
vec <- as.vector(sce_exp_numeric)

# Create the histogram
hist(vec)

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

sce_exp_filt_df_wona = sce_exp_filt_df_wona[order(sce_exp_filt_df_wona$pam50),]

annotation_df <- data.frame(pam50 = sce_exp_filt_df_wona$pam50)
rownames(annotation_df) <- rownames(sce_exp_filt_df_wona)

# Create a named color vector for the unique values of vital_status
pam50_colors <- c("Basal" = "yellow", "Her2" = "black", "LumA" = "green", "LumB" = "red","Normal" = "blue" )
names(pam50_colors) <- unique(annotation_df$pam50)

png(filename = paste0(CancerType,"_exp_clusterT_pam50.png"),
    width = 25, height = 25,  units = "cm" ,pointsize = 12,
    bg = "white", res = 1200, family = "")

tmp = pheatmap::pheatmap(as.matrix(t(sce_exp_filt_df_wona[,-ncol(sce_exp_filt_df_wona)])),
                         annotation_col = annotation_df,
                         annotation_colors = list(pam50 = pam50_colors),
                         cluster_cols = T)
print(tmp)

dev.off()

png(filename = paste0(CancerType,"_exp_clusterF_pam50.png"),
    width = 25, height = 25,  units = "cm" ,pointsize = 12,
    bg = "white", res = 1200, family = "")

tmp2 = pheatmap::pheatmap(as.matrix(t(sce_exp_filt_df_wona[,-ncol(sce_exp_filt_df_wona)])),
                         annotation_col = annotation_df,
                         annotation_colors = list(pam50 = pam50_colors),
                         cluster_cols = F)
print(tmp2)

dev.off()

png(filename = paste0(CancerType,"_exp_cluster_complex_pam50.png"),
    width = 25, height = 25,  units = "cm" ,pointsize = 12,
    bg = "white", res = 1200, family = "")
tmp3 = ComplexHeatmap::pheatmap(as.matrix(t(sce_exp_filt_df_wona[,-ncol(sce_exp_filt_df_wona)])),
                                column_split = annotation_df$pam50,
                                annotation_col = annotation_df,
                                annotation_colors = list(pam50 = pam50_colors),
                                cluster_cols = T)

print(tmp3)
dev.off()

## mut
mut_filtered <- subsetMaf(maf = mut, query = "Variant_Classification != 'Nonsense_Mutation'")

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
annotation_df <- data.frame(pam50 = mut_count_filtered_tdf$pam50)
rownames(annotation_df) <- rownames(mut_count_filtered_tdf)

# Create a named color vector for the unique values of vital_status
pam50_colors <- c("Basal" = "yellow", "Her2" = "black", "LumA" = "green", "LumB" = "red","Normal" = "blue" )
names(pam50_colors) <- unique(annotation_df$pam50)


png(filename = paste0(CancerType,"_mut_clusterT_pam50.png"),
    width = 25, height = 25,  units = "cm" ,pointsize = 12,
    bg = "white", res = 1200, family = "")

tmp = pheatmap::pheatmap(as.matrix(t(mut_count_filtered_tdf[,-ncol(mut_count_filtered_tdf)])),
                         annotation_col = annotation_df,
                         annotation_colors = list(pam50 = pam50_colors),
                         cluster_cols = T)
print(tmp)

dev.off()

png(filename = paste0(CancerType,"_mut_clusterF_pam50.png"),
    width = 25, height = 25,  units = "cm" ,pointsize = 12,
    bg = "white", res = 1200, family = "")

tmp2 = pheatmap::pheatmap(as.matrix(t(mut_count_filtered_tdf[,-ncol(mut_count_filtered_tdf)])),
                         annotation_col = annotation_df,
                         annotation_colors = list(pam50 = pam50_colors),
                         cluster_cols = F)
print(tmp2)

dev.off()

png(filename = paste0(CancerType,"_mut_cluster_complex_pam50.png"),
    width = 25, height = 25,  units = "cm" ,pointsize = 12,
    bg = "white", res = 1200, family = "")
tmp3 = ComplexHeatmap::pheatmap(as.matrix(t(mut_count_filtered_tdf[,-ncol(mut_count_filtered_tdf)])),
                                column_split = annotation_df$pam50,
                                annotation_col = annotation_df,
                                annotation_colors = list(pam50 = pam50_colors),
                                cluster_cols = T)

print(tmp3)
dev.off()
