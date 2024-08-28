library(TCGAbiolinks)

filepath = "/home/seokwon/nas/"

Cancerlist = dir(paste0(filepath, "/00.data/filtered_TCGA/"))
Cancerlist = Cancerlist[c(-11,-12)]

fig_path = paste0(filepath,"/99.reference/TCGA_clinical_drug/")
if(!dir.exists(fig_path)){
  dir.create(fig_path)
  print(paste0("Created folder: ", fig_path))
} else {
  print(paste0("Folder already exists: ", fig_path))
}

setwd(fig_path)

for (num_CancerType in Cancerlist) {
  
  main.path_tc = paste0(filepath, "00.data/filtered_TCGA/", num_CancerType)
  CancerType = gsub('[.]','',gsub('\\d','', num_CancerType))
  Cancername = gsub('TCGA-' , '', CancerType)
  
  query <- GDCquery(
    project = CancerType, 
    data.category = "Clinical",
    data.type = "Clinical Supplement", 
    data.format = "BCR Biotab"
  )
  
  GDCdownload(query)
  clinical.BCRtab.all <- GDCprepare(query)
 
  drug_info = as.data.frame(  clinical.BCRtab.all[[grep("drug", names(clinical.BCRtab.all))]])
  setwd(fig_path)

  write.csv(drug_info , paste0(Cancername , "_clinical_druginfo.csv"))
  
}




