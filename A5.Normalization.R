#!/usr/bin/env Rscript
rm(list = ls())
gc()

library(dplyr)
library(Seurat)
library(ggplot2)

setwd("G://ZXL_221212")

scRNA <- readRDS(file = "A4.QualityControl//R2.SeuObjAfterQC.rds")

if(file.exists("A5.Normalization") == TRUE){
  print("This directory already exists.")
}else{
  dir.create("A5.Normalization")
}
setwd("A5.Normalization")

scRNA <- NormalizeData(scRNA)
scRNA  <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(scRNA), 10)
top10 
plot1 <- VariableFeaturePlot(scRNA)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
ggsave("D1.VariableFeature.pdf", plot = plot2, height = 8, width = 12)

all.genes <- rownames(scRNA)
scRNA <- ScaleData(scRNA, features = all.genes, vars.to.regress = c("percent.Rb", "percent.mt"))
saveRDS(scRNA, file = "R3.SeuObjAfterScale.rds")	
scRNA <- RunPCA(scRNA, features = VariableFeatures(object = scRNA))

plot3 <- VizDimLoadings(scRNA, dims = 1:9, reduction = "pca") & 
  theme(axis.text=element_text(size=5), axis.title=element_text(size=8,face="bold"))
ggsave("D2.VizDimLoadings.pdf", plot = plot3, height = 12, width = 12)

DimHeatmap(scRNA, dims = 1:9, nfeatures = 20, cells = 500, balanced = T, combine = TRUE)

plot4 <- DimPlot(scRNA, reduction = "pca")
ggsave("D3.PCADimPlot.pdf", plot = plot4, height = 12, width = 12)

plot5 <- ElbowPlot(scRNA, ndims = 50, reduction = "pca")
ggsave("E1.ElbowPlot.pdf", plot = plot5, height = 10, width = 12)

scRNA <- RunTSNE(scRNA, reduction = "pca", dims = 1:30, verbose = F)
scRNA <- RunUMAP(scRNA, reduction = "pca", dims = 1:30, verbose = F)

plot6 <- DimPlot(scRNA,reduction = "umap")
ggsave("D4.UMAPDimBeforeHar.pdf", plot = plot6, height = 12, width = 12)

saveRDS(scRNA, file = "R4.SeuObjBeforeHarmony.rds")