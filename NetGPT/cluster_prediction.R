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

out = pheatmap(phyper.dat.2, cluster_cols = T,
               cluster_rows = F, labels_cols = "",
               show_rownames = T)


remove(k, kk,total_results)
total_results =""

df.results = data.frame()


for (cutoff in seq(0.1,0.9,0.1)) {
  for (k in 3:10) {
    #print(k)
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
    
    # ggsurvplot(fit, data = cli_surv_filt, risk.table = TRUE,
    #            palette = "jco", pval = TRUE, surv.median.line = "hv", xlab = "month")
    # ########################################################################################
    ########## survival rate classification ###############################################
    ########################################################################################
    cluster_bad = ""
    cluster_good = ""
    df.cluster = data.frame()
    start_cluster_num = 0 
    end_cluster_num = fit[["strata"]][["cluster=1"]]
    # print("hear")
    kk=k
    cluster_surv = vector()
    for (cluster_k in 1 : kk) {
      # print(cluster_k)
      cluster_surv = c(cluster_surv,fit[["surv"]][end_cluster_num])
      start_cluster_num = start_cluster_num + 1
      if (sum(fit[["surv"]][start_cluster_num:end_cluster_num] < cutoff) != 0) { # here can be changed the value of survival rate to 
        #print( "bad")
        df.cluster = rbind(df.cluster, "bad")
        #cluster_bad = paste0(cluster_bad, "cluster", as.character(cluster_k) )
      } else {
        #print("good")
        df.cluster = rbind(df.cluster, "good")
        #cluster_good = paste0(cluster_good, "cluster", as.character(cluster_k))
        
      } 
      
      if (cluster_k == kk) {
        break
      } else {
        start_cluster_num = end_cluster_num 
        end_cluster_num = end_cluster_num + fit[["strata"]][[paste0("cluster=", as.character(cluster_k+1))]]
      }
      
      
      
    }
    
    dimnames(df.cluster) <- list(paste("cluster", seq(1:kk)), "classification")
    
    if (length(unique(df.cluster$classification)) != 2 && unique(df.cluster$classification) == "bad") {
      df.cluster[which(cluster_surv == max(cluster_surv)),] = "good"
    } else if (length(unique(df.cluster$classification)) != 2 && unique(df.cluster$classification) == "good") {
      df.cluster$classification = "bad"
      df.cluster[which(cluster_surv == max(cluster_surv)),] = "good"
    }
    
    #df.cluster
    
    # print("hear")
    ###########################################################################################
    ### get submitter_id of good or bad cluster ###############################################
    ###########################################################################################
    
    good.pat = subset(df.cluster, df.cluster$classification == 'good')
    bad.pat = subset(df.cluster, df.cluster$classification == 'bad')
    
    tmp2.subcluster = data.frame()
    
    df.bad = data.frame()
    for (num_cluster in gsub('\\D',"",rownames(bad.pat))) {
      a1.subcluster = subset(subcluster.dat.filt ,subcluster.dat.filt$cluster == num_cluster)
      tmp2.subcluster = a1.subcluster
      df.bad= rbind(df.bad,tmp2.subcluster)
      remove(tmp2.subcluster, a1.subcluster)
      
    }
    #df.bad
    
    df.good = data.frame()
    for (num_cluster in gsub('\\D',"",rownames(good.pat))) {
      a1.subcluster = subset(subcluster.dat.filt ,subcluster.dat.filt$cluster == num_cluster)
      tmp2.subcluster = a1.subcluster
      df.good= rbind(df.good,tmp2.subcluster)
      remove(tmp2.subcluster, a1.subcluster)
    }
    #df.good
    
    ########################################################################################
    #########set reference graph for good and bad cluster ##################################
    ########################################################################################
    
    df.phyper.bm = as.data.frame(phyper.dat.bm)
    
    ### good ###
    
    df.good.filt = data.frame(matrix(nrow = length(rownames(df.phyper.bm))))
    tmp.cl = data.frame()
    # df.phyper.bm = as.data.frame(phyper.dat.bm)
    for (gd.id in df.good$submitter_id) {
      
      tmp.cl= df.phyper.bm[which(colnames(df.phyper.bm) == gd.id)]
      tmp2.cl = tmp.cl
      df.good.filt = cbind(df.good.filt,tmp2.cl)
      
      remove(tmp.cl,gd.id)
    }
    
    df.good.filt = df.good.filt[,-1]
    # df.good.filt
    ### bad ###
    
    df.bad.filt = data.frame(matrix(nrow = length(rownames(df.phyper.bm))))
    tmp.cl = data.frame()
    
    for (bd.id in df.bad$submitter_id) {
      
      tmp.cl= df.phyper.bm[which(colnames(df.phyper.bm) == bd.id)]
      tmp2.cl = tmp.cl
      df.bad.filt = cbind(df.bad.filt,tmp2.cl)
      
      remove(tmp.cl,bd.id)
    }
    
    df.bad.filt = df.bad.filt[,-1]
    #df.bad.filt
    
    #subcluster.dat.filt
    ########################################################################################
    ########### ratio calculation ##################################################
    ########################################################################################
    
    # df.bad.filt
    
    bad.ratio= as.data.frame(rowSums(df.bad.filt) / ncol(df.bad.filt))
    good.ratio =as.data.frame(rowSums(df.good.filt)/ ncol(df.good.filt))
    
    total.ratio = cbind(bad.ratio,good.ratio)
    colnames(total.ratio) = c("bad","good")
    total.ratio = as.data.frame(total.ratio)
    
    ## it would be required because of maintaining dataframe form.
    
    total.ratio$difference = total.ratio$bad - total.ratio$good
    total.ratio$ab.diff = abs(total.ratio$difference)
    
    total.ratio$judge = total.ratio$ab.diff < sd(total.ratio$ab.diff)
    
    # total.ratio$spe = ifelse(total.ratio$judge,
    #                          NA,
    #                          total.ratio$difference)
    # total.ratio.spe = na.omit(total.ratio)
    
    ### criteria changed  negative value to  individual specificity of each column 
    # data collect below 1Q 
    # 
    
    total.ratio.good = total.ratio[which(total.ratio$good < fivenum(total.ratio$good)[2]),]
    
    total.ratio.bad = total.ratio[which(total.ratio$bad < fivenum(total.ratio$bad)[2]),]
    total.ratio.good = total.ratio.good[,-1]
    total.ratio.bad = total.ratio.bad [,-2]
    # if (sum(rownames(total.ratio.bad) %in% rownames(total.ratio.good)) != 0) {
    #   
    #   if (length(total.ratio.good) < length(total.ratio.bad)) {
    #     if ( length(total.ratio.good) == 1) {
    #       total.ratio.bad = total.ratio.bad[-which(rownames(total.ratio.bad) %in% rownames(total.ratio.good)),]
    #     } else {
    #       
    #     }
    #   }
    # }
    
    correct = 0
    wrong = 0 
    for (cluster_num in 1:kk) {
      #print(kk)
      
      for (clu_sub_id in cli_surv_filt[which(cli_surv_filt$cluster == cluster_num),]$submitter_id) {
        pat.phyper.dat.bm = as.data.frame(phyper.dat.bm[,clu_sub_id])
        #print(pat.phyper.dat.bm)
        colnames(pat.phyper.dat.bm) = clu_sub_id
        
        ratio.good = 0
        ratio.bad = 0
        
        for (gd.p in rownames(total.ratio.good)) {
          if (pat.phyper.dat.bm[gd.p,] == 0) {
            ratio.good = ratio.good + 1 
            
          }
          
        }
        
        for (bd.p in rownames(total.ratio.bad)) {
          if (pat.phyper.dat.bm[bd.p,] == 0) {
            ratio.bad = ratio.bad + 1 
            
          }
          
        }
        ai = NA
        per.good = ratio.good / length(rownames(total.ratio.good))  
        per.bad = ratio.bad / length(rownames(total.ratio.bad)) 
        if (per.good > per.bad) {
          ai = "Alive"
        } else if (per.good < per.bad) {
          ai = "Dead"
        }
        
        if (is.na(ai)) {
          next
        } else if (cli_surv_filt[which(cli_surv_filt$submitter_id == clu_sub_id),]$vital_status == ai) {
          correct = correct + 1
        } else {
          wrong = wrong + 1
        }
        
        ai = NA
        per.good =0
        per.bad = 0
        ratio.bad = 0
        ratio.good =0
        
      }
      
      
    }
    
    #cli_surv_filt$submitter_id
    
    #ratio.bad
    correct.ratio= correct/length(cli_surv_filt$submitter_id)
    #correct.ratio
    wrong.ratio = wrong/length(cli_surv_filt$submitter_id)
    
    df.tmp= data.frame(k,correct.ratio, wrong.ratio)
    df.results = rbind(df.results, df.tmp)
    #df.classification = data.frame(correct.ratio,wrong.ratio, row.names = CancerType)
    
    wrong.ratio = 0
    correct.ratio = 0
  }
  df.results
  
  write.csv(df.results, paste0("clusters_prediction_results_",cutoff,".csv"), row.names = F)
  assign(paste0("results_",cutoff),read.csv(paste0("clusters_prediction_results_",cutoff,".csv"), row.names=1))
  df.results = data.frame()
  df.tmp = data.frame()
}


max_cut = vector()
tmp_cut = 0
for (cut_num2 in seq(0.1,0.9,0.1)) {
  
  tmp_cut = max(get(paste0("results_",cut_num2)))
  
  max_cut = c(max_cut, tmp_cut)
}

names(max_cut) = c(seq(0.1,0.9,0.1))

best_results = get(paste0("results_",names(max_cut[which(max_cut == max(max_cut))])))

in.max.ratio = which(get(paste0("results_",names(max_cut[which(max_cut == max(max_cut))]))) == max(max_cut))
cl_max = length(get(paste0("results_",names(max_cut[which(max_cut == max(max_cut))])))$correct.ratio)
best_cutoff = names(max_cut[which(max_cut == max(max_cut))])

if (length(in.max.ratio) != 1 ) {
  in.max.ratio = in.max.ratio[1]
  if ( in.max.ratio > cl_max) {
    best_results_out = best_results[in.max.ratio - cl_max,]
    best_results_out$best_cutoff = best_cutoff
  } else {
    best_results_out = best_results[in.max.ratio,]
  }
} else {
  if ( in.max.ratio > cl_max) {
    best_results_out = best_results[in.max.ratio - cl_max,]
    best_results_out$best_cutoff = best_cutoff
  } else {
    best_results_out = best_results[in.max.ratio,]
  }
  
}

best_results_out

