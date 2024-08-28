meta.tbl <- read.csv(file = '/mnt/gluster_server/data/raw/DepMap/23Q2/Model.csv', 
                     sep=',', header=T, fill=T)

total_meta = data.frame()
for (oncotree in unique(meta.tbl$OncotreeLineage)) {
  tmp_df = meta.tbl[which(meta.tbl$OncotreeLineage == oncotree & meta.tbl$PrimaryOrMetastasis == 'Primary' & meta.tbl$GrowthPattern != 'Organoid'),]
  sample_collection = unique(tmp_df$SampleCollectionSite)
  oncotreecode = unique(tmp_df$OncotreeCode)
  
  tmp_filt_meta = data.frame(oncotree = rep(oncotree , max(length(sample_collection) , length(oncotreecode))))
  if (length(sample_collection) > length(oncotreecode)) {
    length(oncotreecode) = length(sample_collection)
  } else {
    length(sample_collection) = length(oncotreecode)
  }
                                              
  tmp_filt_meta = cbind(tmp_filt_meta , sample_collection = sample_collection , oncotreecode = oncotreecode)
  
  total_meta = rbind(total_meta, tmp_filt_meta)
}


write_csv(total_meta, "~/nas/99.reference/Depmap_meta_filt.csv")
