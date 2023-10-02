#!/usr/bin/env Rscript
rm(list = ls())
gc()

library(dplyr)
library(Seurat)
library(patchwork)
library(stringr)
library(DropletUtils)

setwd("G://ZXL_221212//Data//RowData")

if(file.exists("G://ZXL_221212//Data//SplitMatrix") == TRUE){
  print("This directory already exists.")
}else{
  dir.create("G://ZXL_221212//Data//SplitMatrix")
}

sample_info <- list.files(getwd())

for(num in 1:length(sample_info)){
  sample_name <-  str_sub(sample_info[num], start = 3, end = nchar(sample_info[num]) - 4)
  sample <- readRDS(sample_info[num])
  write10xCounts(paste0("G://ZXL_221212//Data//SplitMatrix//",sample_name,"//"),sample[["RNA"]]@counts,version = "3")
}

#python Scrublet.py