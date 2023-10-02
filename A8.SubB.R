#!/usr/bin/env Rscript
rm(list = ls())
gc()

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(harmony)

setwd("G://ZXL_221212")
load('A7.SubCluster//2.Rdata')

if(file.exists("A8.SubB") == TRUE){
  print("This directory already exists.")
}else{
  dir.create("A8.SubB")
}

setwd("A8.SubB")
set.seed(1234)

NorMethod <- "SCT"
#NorMethod <- "RNA"
pc.num = 1:15

if(NorMethod == "SCT") {
  sce <- sce %>%
    SCTransform(vars.to.regress = c("percent.Rb", "percent.mt")) %>%
    RunPCA(assay="SCT", verbose = FALSE) %>%
    RunHarmony(group.by.vars="orig.ident", assay.use="SCT", max.iter.harmony = 20) %>%
    RunUMAP(reduction = "harmony", dims = pc.num) %>%
    FindNeighbors(reduction = "harmony", dims = pc.num) %>%
    FindClusters(resolution = 0.4)
} else if (NorMethod == "RNA") {
  sce <- sce %>%
    NormalizeData(verbose = F) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000, verbose = F) %>%
    ScaleData(verbose = F) %>%
    RunPCA(npcs = 30, verbose = F) %>%
    RunUMAP(reduction = "pca", dims = 1:30, verbose = F) %>%
    FindClusters(resolution = 0.4)
}





if(FLASE){
  table(sce@meta.data$seurat_clusters, sce@meta.data$orig.ident)  
  cellNum <- as.data.frame(table(sce@meta.data$SCT_snn_res.0.4, sce@meta.data$orig.ident))
  library(tidyr)
  library(reshape2)
  cellNum <- spread(cellNum, Var2, Freq)
  row.names(cellNum) <- cellNum[,1]
  cellNum <- cellNum[,-1]
  write.csv(cellNum, "Int.cellNum.csv")
  
  library(qiime2R)
  cellPercent <- as.data.frame(make_percent(cellNum))
  write.csv(cellPercent, "Int.cellPercent.csv")
  
  D <- DimPlot(sce, reduction = "umap", label = T)
  D1 <- DimPlot(sce, reduction = "umap", label = T, split.by = "group")
  
  ggsave("Int.BCellDimplot.pdf", D, height = 9, width = 10)
  ggsave("Int.BCellDimpotGroup.pdf", D1, height = 9, width = 30)
  
  
  sce <- PrepSCTFindMarkers(object = sce, assay = "SCT", verbose = TRUE)
  sce.markers <- FindAllMarkers(sce, 
                                only.pos = TRUE, 
                                min.pct = 0.25, 
                                logfc.threshold = 0.25)
  
  sce.markers %>%
    group_by(cluster) %>%
    slice_max(n = 5, order_by = avg_log2FC) %>%
    print(n = 35)
  
  write.csv(sce.markers, "C1.BMarker.csv")
}

sce <- subset(sce, seurat_clusters != 6 & seurat_clusters != 9)
sce <- RunPCA(sce, assay="SCT", verbose = FALSE)
sce <- RunHarmony(sce, group.by.vars="orig.ident", assay.use="SCT", max.iter.harmony = 20) 
sce <- RunUMAP(sce, reduction = "harmony", dims = 1:10)
sce <- FindNeighbors(sce, reduction = "harmony", dims = 1:10)
sce <- FindClusters(sce, resolution = 0.4)

DimPlot(sce, reduction = "umap", label = T)

D <- DimPlot(sce, reduction = "umap", label = T)
D1 <- DimPlot(sce, reduction = "umap", label = T, split.by = "group")

ggsave("Int2.BCellDimplot.pdf", D, height = 9, width = 10)
ggsave("Int2.BCellDimpotGroup.pdf", D1, height = 9, width = 30)

sce <- PrepSCTFindMarkers(object = sce, assay = "SCT", verbose = TRUE)
sce.markers <- FindAllMarkers(sce, 
                              only.pos = TRUE, 
                              min.pct = 0.25, 
                              logfc.threshold = 0.25)

sce.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC) %>%
  print(n = 35)

write.csv(sce.markers, "C2.BMarker.csv")

markers <- FindMarkers(sce, ident.1 = 0, ident.2 = 1, min.pct = 0.25)
write.csv(markers, "C2.01.csv")

markers <- FindMarkers(sce, ident.1 = 3, ident.2 = 4, min.pct = 0.25)
write.csv(markers, "C2.34.csv")

markers <- FindMarkers(sce, ident.1 = 6, ident.2 = c(0, 2, 3, 5), min.pct = 0.25)
write.csv(markers, "C2.6.0235.csv")

library(RColorBrewer)

saveRDS(sce, file = "R7.SubB.rds")

FeaturePlot(immune.combined, "CD79A", label = T, order = T) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))

FeaturePlot(sce, "IGKC", label = T, order = T) + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))

FeaturePlot(sce, "NEAT1", label = T, order = T) + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))

FeaturePlot(sce, "FCER2", label = F, order = T) + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))

FeaturePlot(sce, "IGHD", label = T, order = T) + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))

FeaturePlot(sce, "CD27", label = T, order = T) + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))

VlnPlot(sce, c("TCL1A", "IGHD", "IL4R", "AIM2",
               "TNFRSF13B", "CD27", "MZB1",
               "IGHG3", "JCHAIN"), pt.size = 0)
