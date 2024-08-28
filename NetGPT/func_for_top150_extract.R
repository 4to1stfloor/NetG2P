library(pheatmap)
library(survminer) 
library(ggplot2) 
library(pheatmap)
library(survival) 
library(data.table) 
library(RColorBrewer) 
library(NbClust)
library(cluster) 
library(factoextra) 
library(dplyr)
library(stringr)
library(h2o)
h2o.shutdown(prompt = F)
filepath = "/home/seokwon/"
ref_path = paste0(filepath, "99.reference/")
Cancername = "TCGA-BRCA" # you can change this
main.path_tc = paste0(filepath, "00.data/","30.", Cancername,"/")

localH2O = h2o.init(ip = "localhost", port = sample(x=00000:65536,size=1,replace=F), startH2O = TRUE,min_mem_size = "50G",nthreads = 96,enable_assertions = FALSE)

best_model = h2o.loadModel(path = paste0(main.path_tc , "dl_grid_model_13_0.930851063829787/dl_grid_model_13"))

cut_top_num = c(5,6,7,8,9,10,25,50,100,150)
for (top_cut in cut_top_num) {
  assign(paste0("top_",top_cut,"_vari"),as.data.frame(h2o.varimp(best_model))$variable[1:top_cut]) 
}

phyper = readRDS(paste0(main.path_tc, "/phyper_enrichment_54_KEGG_net_GC_20220117.rds"))
filt_phyper = readRDS(file = paste0(main.path_tc,"/",Cancername,"_pathwaylink_all_0or1.rds"))

filt_phyper = filt_phyper[,-ncol(filt_phyper)]
filt_phyper_t = t(filt_phyper)

for (phyper_num in cut_top_num) {
  assign(paste0("filt_phyper_top_",phyper_num),filt_phyper_t[get(paste0("top_",phyper_num,"_vari")),]) 
  
}

for (df_num in cut_top_num) {
  assign(paste0("df_top",df_num,"_vari"),as.data.frame(h2o.varimp(best_model)[1:df_num,])) 
}


anno_sig <- function(top_vari) {
  
  top_vari$sig = NA
  top_vari[which(top_vari$relative_importance > quantile(top_vari$relative_importance, 0.9) ),]$sig = "10%"
  top_vari[which(quantile(top_vari$relative_importance, 0.9) > top_vari$relative_importance &
                   top_vari$relative_importance > quantile(top_vari$relative_importance, 0.8) ),]$sig = "20%"
  top_vari[which(quantile(top_vari$relative_importance, 0.8) > top_vari$relative_importance &
                   top_vari$relative_importance > quantile(top_vari$relative_importance, 0.7) ),]$sig = "30%"
  top_vari[which(quantile(top_vari$relative_importance, 0.7) > top_vari$relative_importance &
                   top_vari$relative_importance > quantile(top_vari$relative_importance, 0.6) ),]$sig = "40%"
  top_vari[which(quantile(top_vari$relative_importance, 0.6) > top_vari$relative_importance &
                   top_vari$relative_importance > quantile(top_vari$relative_importance, 0.5) ),]$sig = "50%"
  top_vari[which(is.na(top_vari$sig)),]$sig = 0
  return(top_vari)
  
}


for (df_num_filt in cut_top_num) {
  if (df_num_filt %in% c(5,6,7,8,9)) {
    print(df_num_filt)
  } else {
    assign(paste0("df_top",df_num_filt,"_vari_filt"), anno_sig(get(paste0("df_top",df_num_filt,"_vari"))))
  }
  
}
df_top9_vari$sig = NA
df_top9_vari$sig = c("10%","20%","30%","40%","50%","other","other","other","other")
df_top9_vari_filt = df_top9_vari

df_top8_vari$sig = NA
df_top8_vari$sig = c("10%","20%","30%","40%","50%","other","other","other")
df_top8_vari_filt = df_top8_vari

df_top7_vari$sig = NA
df_top7_vari$sig = c("10%","20%","30%","40%","50%","other","other")
df_top7_vari_filt = df_top7_vari

df_top6_vari$sig = NA
df_top6_vari$sig = c("10%","20%","30%","40%","50%","other")
df_top6_vari_filt = df_top6_vari

df_top5_vari$sig = NA
df_top5_vari$sig = c("10%","20%","30%","40%","50%")
df_top5_vari_filt = df_top5_vari

########################
## clinical data upload#
########################

cli = fread(paste0(ref_path,"all_clin_indexed.csv"))

cli_surv = cli[cli$project == Cancername,
               c("submitter_id",
                 "vital_status",
                 "days_to_death",
                 "days_to_last_follow_up")]
                    
# cli_surv_COAD = cli[cli$project == "TCGA-COAD",
#                c("submitter_id",
#                  "vital_status",
#                  "days_to_death",
#                  "days_to_last_follow_up")]
# cli_surv_READ = cli[cli$project == "TCGA-READ",
#                     c("submitter_id",
#                       "vital_status",
#                       "days_to_death",
#                       "days_to_last_follow_up")]
# cli_surv = rbind(cli_surv_COAD,cli_surv_READ)
cli_surv$deceased = cli_surv$vital_status == "Dead"
cli_surv$overall_survival = ifelse(cli_surv$deceased,
                                   cli_surv$days_to_death,
                                   cli_surv$days_to_last_follow_up)

cli_surv = cli_surv[!is.na(cli_surv$vital_status),]

cli_surv$status = NA
cli_surv$status[which(cli_surv$vital_status == "Dead")] = 1
cli_surv$status[which(cli_surv$vital_status == "Alive")] = 0

cli_surv$status = as.numeric(cli_surv$status)


# origin

setwd(paste0(main.path_tc,"topfeature_fig/"))
top_num = 150
for (top_num in cut_top_num) {
  
  filt_phyper_tmp = get(paste0("filt_phyper_top_",as.character(top_num)))
  df_top_vari_filt = get(paste0("df_top",as.character(top_num),"_vari_filt"))
  df_top_vari = get(paste0("df_top",as.character(top_num),"_vari"))
  
  out = pheatmap(filt_phyper_tmp, cluster_cols = T,
                cluster_rows = T, labels_cols = "",
                show_rownames = T, show_colnames = F)
  
  remove(dat.comb)
  for (k in 2:10) {
    
    subcluster.dat = as.data.frame(cutree(out$tree_col, k), row.names = out[["tree_col"]][["labels"]])
    colnames(subcluster.dat) = "cluster"
    subcluster.dat$submitter_id = substr(rownames(subcluster.dat), 1, 12)
    subcluster.dat = subcluster.dat[order(subcluster.dat$cluster, decreasing = FALSE),]
    
    pat.int = intersect(subcluster.dat$submitter_id, cli_surv$submitter_id)
    pat.idx = match(pat.int, subcluster.dat$submitter_id)
    pat.idx2 = match(pat.int, cli_surv$submitter_id)
    
    subcluster.dat.filt = subcluster.dat[pat.idx, ]
    cli_surv_filt = cli_surv[pat.idx2,]
    
    remove(pat.int, pat.idx, pat.idx2)
    
    cli_surv_filt$cluster = subcluster.dat.filt$cluster
    
    fit = survfit(Surv(overall_survival, status) ~ cluster, data = cli_surv_filt)
    dat = as.data.frame(cbind(k_values = k, pval = surv_pvalue(fit)[2]))
    #print(surv_pvalue(fit)[2])
    remove (fit)  
    # ggsurvplot(fit, data = cli_surv_filt, risk.table = TRUE,
    #            palette = "jco", pval = TRUE, surv.median.line = "hv", xlab = "month")
    if ( k == 2 ){
      dat.comb = dat
    } else {
      dat.comb = rbind(dat.comb, dat)
    }
  }
  dat.comb
  #########################
  # drawing the survival plot 
  ###########################
  
  which.cluster = dat.comb$k_values[which(dat.comb$pval == min(dat.comb$pval))]
  cat(which.cluster, "\n")
  
  subcluster.dat = as.data.frame(cutree(out$tree_col, which.cluster), row.names = out[["tree_col"]][["labels"]])
  
  colnames(subcluster.dat) = "cluster"
  subcluster.dat$submitter_id = substr(rownames(subcluster.dat), 1, 12)
  subcluster.dat = subcluster.dat[order(subcluster.dat$cluster, decreasing = FALSE),]
  
  pat.int = intersect(subcluster.dat$submitter_id, cli_surv$submitter_id)
  pat.idx = match(pat.int, subcluster.dat$submitter_id)
  pat.idx2 = match(pat.int, cli_surv$submitter_id)
  
  subcluster.dat.filt = subcluster.dat[pat.idx, ]
  cli_surv_filt = cli_surv[pat.idx2,]
  
  remove(pat.int, pat.idx, pat.idx2)
  
  cli_surv_filt$cluster = subcluster.dat.filt$cluster
  
  fit = survfit(Surv(overall_survival, status) ~ cluster, data = cli_surv_filt)
  
  setwd(paste0(main.path_tc,"topfeature_fig/"))
  
  png(filename = paste0("survival_",top_num ,"orgin2.png"),
      width = 1000, height = 1000, units = "px", pointsize = 12,
      bg = "white", res = NA, family = "")
  
  ggsurvplot(fit, data = cli_surv_filt, risk.table = TRUE,
             palette = "jco", pval = TRUE, surv.median.line = "hv", xlab = "month")
  
  dev.off()
  ### alive dead add
  
  flit_phyper_add_anno = as.data.frame(t(filt_phyper_tmp))
  flit_phyper_add_anno$vital_status = NA
  
  for (pat_id in rownames(flit_phyper_add_anno)) {
    flit_phyper_add_anno[pat_id,]$vital_status = cli_surv_filt[which(cli_surv_filt$submitter_id == pat_id),]$vital_status
  }
  
  row_ann = as.data.frame(df_top_vari_filt$sig)
  rownames(row_ann) = df_top_vari$variable
  colnames(row_ann) = "sig"
  
  col_ann = as.data.frame(flit_phyper_add_anno$vital_status)
  colnames(col_ann) = "vital_status"
  rownames(col_ann) = rownames(flit_phyper_add_anno)
  
  subcluster.dat.filt_col = subcluster.dat.filt[order(subcluster.dat.filt$submitter_id),]
  if (unique(rownames(col_ann) == subcluster.dat.filt[order(subcluster.dat.filt$submitter_id),]$submitter_id)) {
    col_ann$cluster = subcluster.dat.filt_col$cluster
    col_ann$cluster = as.character(col_ann$cluster) 
  } else {
    print("It doesnot match")
  }
  
  library(RColorBrewer)
  
  if (length(unique(col_ann$cluster)) == 2) {
    col_colors = list(cluster = c("#FF0000","#00FFFF"), vital_status = c("black", "yellow"), sig = brewer.pal(length(unique(row_ann$sig)),"Dark2"))
  } else {
    # col_colors = list(cluster = brewer.pal(length(unique(col_ann$cluster)), "Set3") )
    col_colors = list(cluster = rainbow(length(unique(col_ann$cluster))), vital_status = c("black", "yellow"), sig = brewer.pal(length(unique(row_ann$sig)),"Dark2"))
  }
  
  names(col_colors$cluster) = sort(unique(col_ann$cluster))
  names(col_colors$vital_status) = unique(col_ann$vital_status)
  names(col_colors$sig) = unique(row_ann$sig)
  
  library(dendsort)
  
  sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
  # mat_cluster_cols <- sort_hclust(hclust(dist(t(filt_phyper_tmp))))
  mat_cluster_rows <- sort_hclust(hclust(dist(filt_phyper_tmp)))
  
  filt_phyper_origin = as.data.frame(t(filt_phyper_tmp))
  
  cli_surv_filt_order = data.frame()
  
  for (patients_id in rownames(filt_phyper_origin) ) {
    cli_surv_filt_order = rbind(cli_surv_filt_order,cli_surv_filt[which(cli_surv_filt$submitter_id == patients_id),])
  }
  
  if (sum(rownames(filt_phyper_origin) == cli_surv_filt_order$submitter_id) == nrow(filt_phyper_origin)) {
    filt_phyper_origin$vitalstatus = cli_surv_filt_order$vital_status
    filt_phyper_origin$cluster = cli_surv_filt_order$cluster
    
  }
 
  saveRDS(filt_phyper_origin, file = paste0(main.path_tc,"topfeature_fig/" ,Cancername , top_num , "_matrix_origin_with_cluster.rds")) 
  
  saveRDS(filt_phyper_origin, file = paste0(main.path_tc,"topfeature_fig/" ,Cancername , top_num , "_matrix_adjust_with_cluster.rds")) 
  
  filt_phyper_cluster_order = filt_phyper_origin[order(filt_phyper_origin$cluster),]
  
  filt_phyper_cluster_order_t = t(filt_phyper_cluster_order[,-c(length(filt_phyper_cluster_order),length(filt_phyper_cluster_order)-1)])
  
  png(filename = paste0("dendsort_with_cluster_top",top_num ,"orgin2.png"),
      width = 2000, height = 2000, units = "px", pointsize = 12,
      bg = "white", res = NA, family = "")
  
  png(filename = paste0("dendsort_with_cluster_top",top_num ,"adjust.png"),
      width = 2000, height = 2000, units = "px", pointsize = 12,
      bg = "white", res = NA, family = "")
  
  pheatmap(filt_phyper_cluster_order_t, 
           labels_cols = "",
           show_rownames = T, 
           show_colnames = F, 
           annotation_col = col_ann ,
           annotation_colors =  col_colors,
           annotation_row = row_ann,
           clustering_method = "average",
           cluster_rows = mat_cluster_rows,
           cluster_cols = F)
  
  # cluster_cols = mat_cluster_cols)
  # cluster_arrange_top150_origin
  
  dev.off()
}

# adjust
subcluster.dat.filt$cluster
subcluster.dat.filt[which(subcluster.dat.filt$cluster == 1 ),]$cluster = 10
subcluster.dat.filt[which(subcluster.dat.filt$cluster == 2 ),]$cluster = 20

subcluster.dat.filt[which(subcluster.dat.filt$cluster == 3 ),]$cluster = 2
subcluster.dat.filt[which(subcluster.dat.filt$cluster == 4 ),]$cluster = 1
subcluster.dat.filt[which(subcluster.dat.filt$cluster == 5 ),]$cluster = 1
subcluster.dat.filt[which(subcluster.dat.filt$cluster == 6 ),]$cluster = 2
subcluster.dat.filt[which(subcluster.dat.filt$cluster == 7 ),]$cluster = 1
subcluster.dat.filt[which(subcluster.dat.filt$cluster == 8 ),]$cluster = 2
subcluster.dat.filt[which(subcluster.dat.filt$cluster == 9 ),]$cluster = 1

# subcluster.dat.filt[which(subcluster.dat.filt$cluster == 10 ),]$cluster = 5

subcluster.dat.filt[which(subcluster.dat.filt$cluster == 10 ),]$cluster = 1
subcluster.dat.filt[which(subcluster.dat.filt$cluster == 20 ),]$cluster = 1

cli_surv_filt_adjust = cli_surv_filt
cli_surv_filt_adjust$cluster = subcluster.dat.filt$cluster
fit_adjust = survfit(Surv(overall_survival, status) ~ cluster, data = cli_surv_filt_adjust)
ggsurvplot(fit_adjust, data = cli_surv_filt_adjust, risk.table = TRUE,conf.int = FALSE,
           palette = "jco", pval = TRUE ,pval.method = TRUE ,surv.median.line = "hv", xlab = "month")

# fin 


# test 
set.seed(123)
mat2 = matrix(rnorm(50*50), nrow = 50)
split = rep(1:5, each = 10)

ha = HeatmapAnnotation(
  empty = anno_empty(border = FALSE, height = unit(8, "mm")),
  foo = anno_block(gp = gpar(fill = 2:6), labels = LETTERS[1:5])
)
Heatmap(filt_phyper_cluster_order_t, name = "filt_phyper_cluster_order_t", column_split = split, top_annotation = ha, 
        column_title = NULL)

library(GetoptLong)  # for the function qq()
group_block_anno = function(group, empty_anno, gp = gpar(), 
                            label = NULL, label_gp = gpar()) {
  
  seekViewport(qq("annotation_@{empty_anno}_@{min(group)}"))
  loc1 = deviceLoc(x = unit(0, "npc"), y = unit(0, "npc"))
  seekViewport(qq("annotation_@{empty_anno}_@{max(group)}"))
  loc2 = deviceLoc(x = unit(1, "npc"), y = unit(1, "npc"))
  
  seekViewport("global")
  grid.rect(loc1$x, loc1$y, width = loc2$x - loc1$x, height = loc2$y - loc1$y, 
            just = c("left", "bottom"), gp = gp)
  if(!is.null(label)) {
    grid.text(label, x = (loc1$x + loc2$x)*0.5, y = (loc1$y + loc2$y)*0.5, gp = label_gp)
  }
}

group_block_anno(1:3, "empty", gp = gpar(fill = "red"), label = "group 1")
group_block_anno(4:5, "empty", gp = gpar(fill = "blue"), label = "group 2")




