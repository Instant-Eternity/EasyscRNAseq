#!/usr/bin/env Rscript
rm(list = ls())
gc()

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(harmony)

setwd("G://ZXL_221212")
if(file.exists("A7.SubCluster") == TRUE){
  print("This directory already exists.")
}else{
  dir.create("A7.SubCluster")
}

setwd("A7.SubCluster")

if(FALSE){
  immune.combined <- readRDS("A6.Harmony//R6.SeuObjAfterCluster.rds") 
  
  sce.all.list <- SplitObject(immune.combined , split.by = "seurat_clusters")
  sce.all.list 
  names(sce.all.list)
  for (i in names(sce.all.list)) {
    epi_mat = sce.all.list[[i]]@assays$RNA@counts
    epi_phe = sce.all.list[[i]]@meta.data
    sce=CreateSeuratObject(counts = epi_mat, 
                           meta.data = epi_phe )
    sce
    table(sce@meta.data$orig.ident) 
    save(sce,file = paste0(i,'.Rdata'))
  }
}