#!/usr/bin/env Rscript
rm(list = ls())
gc()

library(dplyr)
library(Seurat)
library(ggplot2)
library(SoupX)

setwd("H://CloudData//01_ZXL_scRNAseq")

Clean_Soup <- function(path){
  out <- load10X(path) %>%
    autoEstCont() %>%
    adjustCounts()
  return (out)
}

File_path <- c('SA1', 'SA2', 'SA3', 'SA4', 'SA5', 'SA6')
Sample_info <- c('1.SLE.1', '1.SLE.2', '1.SLE.3', 
                 '2.SEP.1', '2.SEP.2', '2.SEP.3')

for(num in 1:length(File_path)){
  sample <- Clean_Soup(File_path[num]) %>%
    CreateSeuratObject(project = Sample_info[num], min.cells = 3, min.features = 200)
  saveRDS(sample, paste0("G://ZXL_221212//Data//RowData//",Sample_info[num],".rds"))
}
