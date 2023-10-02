rm(list = ls())
gc()

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(harmony)

setwd("G://ZXL_221212")
load('A7.SubCluster//1.Rdata')

if(file.exists("A9.SubM") == TRUE){
  print("This directory already exists.")
}else{
  dir.create("A9.SubM")
}

setwd("A9.SubM")
set.seed(1234)

sce <- SCTransform(sce)
sce <- RunPCA(sce, assay="SCT", verbose = FALSE)
ElbowPlot(sce)
sce <- RunHarmony(sce, group.by.vars="orig.ident", assay.use="SCT", max.iter.harmony = 20) 
sce <- RunUMAP(sce, reduction = "harmony", dims = 1:15)
sce <- FindNeighbors(sce, reduction = "harmony", dims = 1:15)
sce <- FindClusters(sce, resolution = 0.4)

DimPlot(sce, reduction = "umap", label = T)
DimPlot(sce, reduction = "umap", label = T, split.by = "group")

sce <- subset(sce, seurat_clusters != 11 & seurat_clusters != 12)

sce <- RunPCA(sce, assay="SCT", verbose = FALSE)
ElbowPlot(sce)
sce <- RunHarmony(sce, group.by.vars="orig.ident", assay.use="SCT", max.iter.harmony = 20) 
sce <- RunUMAP(sce, reduction = "harmony", dims = 1:15)
sce <- FindNeighbors(sce, reduction = "harmony", dims = 1:15)
sce <- FindClusters(sce, resolution = 0.4)

DimPlot(sce, reduction = "umap", label = T)

library(RColorBrewer)
  
FeaturePlot(sce, "ITGAX", label = T, order = T) + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
 
VlnPlot(sce, c("FCGR3A", "CD14", "CD1C",
			   "IL3RA", "HLA-DRA"), pt.size = 0)

FeaturePlot(sce, "EEF1G", label = T, order = T) + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral"))) 
 
sce <- PrepSCTFindMarkers(object = sce, assay = "SCT", verbose = TRUE)
sce.markers <- FindAllMarkers(sce, 
                              only.pos = TRUE, 
                              min.pct = 0.25, 
                              logfc.threshold = 0.25)

sce.markers %>%
  group_by(cluster) %>%
  slice_max(n = 3, order_by = avg_log2FC) %>%
  print(n = 33)

write.csv(sce.markers, "C1.MMarker.csv")  