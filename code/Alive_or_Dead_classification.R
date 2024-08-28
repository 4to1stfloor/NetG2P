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

CancerType = "TCGA-COAD"

filepath ="C:/Users/user/Desktop/data/"

surv.pvals = readRDS(paste0(filepath, CancerType,"/surv_pvals_54_KEGG_net_GC_20220117_filt5.rds"))
phyper = readRDS(paste0(filepath, CancerType, "/phyper_enrichment_54_KEGG_net_GC_20220117.rds"))

mut.mat = readRDS(paste0(filepath , CancerType,"/", CancerType,"_mut_count_filt_data.rds"))
mut.count = as.data.frame(colSums(mut.mat))
colnames(mut.count) = "counts"

########
## clinical data upload
####

cli = fread("C:/Users/user/Desktop/data/99.reference/all_clin_indexed.csv")
cli_surv = cli[cli$project == CancerType,
               c("submitter_id",
                 "vital_status",
                 "days_to_death",
                 "days_to_last_follow_up")]

cli_surv$deceased = cli_surv$vital_status == "Dead"
cli_surv$overall_survival = ifelse(cli_surv$deceased,
                                   cli_surv$days_to_death,
                                   cli_surv$days_to_last_follow_up)

cli_surv = cli_surv[!is.na(cli_surv$vital_status),]

cli_surv$status[which(cli_surv$vital_status == "Dead")] = 1
cli_surv$status[which(cli_surv$vital_status == "Alive")] = 0

cli_surv$status = as.numeric(cli_surv$status)
cli_surv
########################
#### 01A filtering #####
########################

mut.count$patient_id = substr(rownames(mut.count), 1, 16)
mut.count$forfilt =  rownames(mut.count)
mut.count.filt = mut.count %>% filter(str_detect(forfilt, "-01A-")) 
mut.filt = mut.count.filt[which(mut.count.filt$counts < 1000),]
mut.filt$forfilt = NULL

# dim(mut.filt)
# colnames(phyper)
mut.filt
phyper
int.id = intersect(colnames(phyper), mut.filt$patient_id)
idx.id = match(int.id, colnames(phyper))
phyper.filt = phyper[,idx.id]

surv = surv.pvals
surv$hallmark = rownames(surv)
surv.sig = surv$hallmark[which(surv[,1] < 0.05)]
surv.sig

##########################################################
#### below code make pvals for 0 or 1 to make binary matrix
###########################################################

for (n in 1:length(surv.sig)) {
  phyper.each = phyper.filt[which(rownames(phyper.filt) == surv.sig[n]), ]
  phyper.dat.tmp = as.data.frame(cbind(id = substr(names(phyper.each), 1, 12), pvals = as.numeric(phyper.each)))
  phyper.dat.tmp$group = ""
  phyper.dat.tmp$group[which(phyper.dat.tmp$pvals <= 0.05)] = 0
  phyper.dat.tmp$group[which(phyper.dat.tmp$pvals > 0.05)] = 1
  if (n == 1) {
    phyper.dat.bm= as.numeric(phyper.dat.tmp$group)
  } else {
    phyper.dat.bm = rbind(phyper.dat.bm, as.numeric(phyper.dat.tmp$group))
  }
}

rownames(phyper.dat.bm) = surv.sig
colnames(phyper.dat.bm) = substr(colnames(phyper.filt), 1,12)

dim(phyper.dat.bm)

for (n in 1:length(surv.sig)) {
  phyper.each = phyper.filt[which(rownames(phyper.filt) == surv.sig[n]), ]
  phyper.dat = as.data.frame(cbind(id = substr(names(phyper.each), 1, 12), pvals = as.numeric(phyper.each)))
  
  if (n == 1) {
    phyper.dat.2= as.numeric(phyper.dat$pvals)
  } else {
    phyper.dat.2 = rbind(phyper.dat.2, as.numeric(phyper.dat$pvals))
  }
}

phyper.dat.2
#phyper.dat.2 = t(as.matrix(phyper.dat.2))
rownames(phyper.dat.2) = surv.sig
colnames(phyper.dat.2) = substr(colnames(phyper.filt),1,12)

##########################################################################
#################### test ################################################
##########################################################################

cli_surv_Alive = cli_surv[which(cli_surv$vital_status == "Alive"),]
cli_surv_Dead = cli_surv[which(cli_surv$vital_status == "Dead"),]
df.phyper.bm = as.data.frame(phyper.dat.bm)

df.phyper.al= data.frame(matrix(nrow = length(rownames(df.phyper.bm))))
for (alive_sub in cli_surv_Alive$submitter_id) {
  if (alive_sub %in% colnames(phyper.dat.bm)) {
    df.phyper.al = cbind(df.phyper.al,df.phyper.bm[alive_sub])
    }
  
}  
df.phyper.al = df.phyper.al[,-1]

df.phyper.dd = data.frame(matrix(nrow = length(rownames(df.phyper.bm))))
for (dead_sub in cli_surv_Dead$submitter_id) {
  if (dead_sub %in% colnames(phyper.dat.bm)) {
    df.phyper.dd = cbind(df.phyper.dd,df.phyper.bm[dead_sub])
  }
  
}  
df.phyper.dd= df.phyper.dd[,-1]

########################################################################################
########### ratio calculation ##########################################################
########################################################################################

# df.bad.filt

bad.ratio= as.data.frame(rowSums(df.phyper.dd) / ncol(df.phyper.dd))
good.ratio =as.data.frame(rowSums(df.phyper.al)/ ncol(df.phyper.al))

total.ratio = cbind(bad.ratio,good.ratio)
colnames(total.ratio) = c("bad","good")
total.ratio = as.data.frame(total.ratio)
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

for (row.name in rownames(total.ratio.good)) {
  if (row.name %in% rownames(total.ratio.bad)) {
    if (total.ratio.bad[row.name,]$bad > total.ratio.good[row.name,]$good) {
      total.ratio.bad = total.ratio.bad[-which(rownames(total.ratio.bad) == row.name),]
    } else if (total.ratio.good[row.name,]$good > total.ratio.bad[row.name,]$bad) {
      total.ratio.good = total.ratio.good[-which(rownames(total.ratio.good) == row.name),]
    }
  }
}
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
df.results = data.frame()

for (clu_sub_id in cli_surv$submitter_id) {
  if (clu_sub_id %in% colnames(phyper.dat.bm)) {
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
    per.good = ratio.good / length(colnames(phyper.dat.bm))  
    per.bad = ratio.bad / length(colnames(phyper.dat.bm)) 
    if (per.good > per.bad) {
      ai = "Alive"
    } else if (per.good < per.bad) {
      ai = "Dead"
    } else {
      print(per.bad)
      print(per.good)
    }
    
    if (is.na(ai)) {
      next
    } else if (cli_surv[which(cli_surv$submitter_id == clu_sub_id),]$vital_status == ai) {
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
correct.ratio= correct/length(colnames(phyper.dat.bm))
#correct.ratio
wrong.ratio = wrong/length(colnames(phyper.dat.bm))
correct.ratio + wrong.ratio
df.tmp= data.frame(correct.ratio, wrong.ratio)

df.results = rbind(df.results, df.tmp)
#df.classification = data.frame(correct.ratio,wrong.ratio, row.names = CancerType)

wrong.ratio = 0
correct.ratio = 0

df.results
write.csv(df.results, paste0("clusters_prediction_results_alordd.csv"), row.names = F)

