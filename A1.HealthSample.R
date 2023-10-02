 #!/usr/bin/env Rscript
rm(list = ls())
gc()

library(dplyr)
library(Seurat)
library(DropletUtils)

setwd("G://ZXL_221024")

data_dir <- "Data//GSE157278"

list.files(data_dir)
data <- Read10X(data.dir = data_dir)
seurat_object <- CreateSeuratObject(counts = data,min.cells = 3, min.features = 200)

seurat_object@meta.data$group = row.names(seurat_object@meta.data)
seurat_object@meta.data$orig.ident = substring(seurat_object@meta.data$group,nchar(seurat_object@meta.data$group[1]))
seurat_object@meta.data <- select(seurat_object@meta.data,-group)

HC_1 <- subset(seurat_object, orig.ident == "1")
HC_3 <- subset(seurat_object, orig.ident == "3")
HC_4 <- subset(seurat_object, orig.ident == "4")

HC_1@meta.data$orig.ident = "0.HC.1"
HC_3@meta.data$orig.ident = "0.HC.3"
HC_4@meta.data$orig.ident = "0.HC.4"

if(file.exists("G://ZXL_221212//Data") == TRUE){
  print("This directory already exists.")
}else{
  dir.create("G://ZXL_221212//Data")
}

if(file.exists("G://ZXL_221212//Data//RowData") == TRUE){
  print("This directory already exists.")
}else{
  dir.create("G://ZXL_221212//Data//RowData")
}

setwd("G://ZXL_221212//Data//RowData")
saveRDS(HC_1,"0.HC.1.rds")
saveRDS(HC_3,"0.HC.3.rds")
saveRDS(HC_4,"0.HC.4.rds")

if(file.exists("G://ZXL_221212//Data//RowMatrix") == TRUE){
  print("This directory already exists.")
}else{
  dir.create("G://ZXL_221212//Data//RowMatrix")
}
setwd("G://ZXL_221212//Data//RowMatrix")

if(FALSE){
  sample_name <- c("0.HC.1", "0.HC.3", "0.HC.4")
  sample_info <- c("HC_1", "HC_3", "HC_4")
  for(num in 1:length(sample_info)){
    name <-  sample_name[num]
    sample <- get(sample_info[num])
    write10xCounts(paste0("G://ZXL_221212//Data//RowMatrix//",name,"//"),sample[["RNA"]]@counts,version = "3")
  }
}