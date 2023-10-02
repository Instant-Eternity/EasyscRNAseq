rm(list = ls())
gc()

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(harmony)

setwd("G://ZXL_221212")
load('A7.SubCluster//0.Rdata')

if(file.exists("A10.SubT") == TRUE){
  print("This directory already exists.")
}else{
  dir.create("A10.SubT")
}

setwd("A10.SubT")
set.seed(2000)

sce <- NormalizeData(sce, normalization.method =  "LogNormalize",  
                     scale.factor = 1e4)
GetAssay(sce,assay = "RNA")
sce <- FindVariableFeatures(sce, 
                            selection.method = "vst", nfeatures = 2000)  
sce <- ScaleData(sce, vars.to.regress = c("percent.Rb", "percent.mt")) 
sce <- RunPCA(object = sce, pc.genes = VariableFeatures(sce)) 

sce <- RunHarmony(sce, "orig.ident")
names(sce@reductions)
Seed = 9000 #2000
sce <- sce %>% 
  RunUMAP(reduction = "harmony", dims = 1:20, verbose = F, seed.use = Seed) %>% 
  FindNeighbors(reduction = "harmony", k.param = 10, dims = 1:20) %>% 
  FindClusters(resolution = 0.2)
 
DimPlot(sce, reduction = "umap", seed = Seed, label = T) 



sce <- SCTransform(sce, vars.to.regress = c("percent.Rb", "percent.mt"))
sce <- RunPCA(sce, assay="SCT", verbose = FALSE)
ElbowPlot(sce)
sce <- RunHarmony(sce, group.by.vars="orig.ident", assay.use="SCT", max.iter.harmony = 20) 

ElbowPlot(sce, ndims = 50, reduction = "harmony")

Seed = 9000 #2000
sce <- sce %>% 
  RunUMAP(reduction = "harmony", dims = 1:20, verbose = F, seed.use = Seed) %>% 
  FindNeighbors(reduction = "harmony", k.param = 10, dims = 1:20) %>% 
  FindClusters(resolution = 0.2)
 
DimPlot(sce, reduction = "umap", seed = Seed, label = T) 


library(RColorBrewer)
  
FeaturePlot(sce, "FOXP3", label = T, order = T) + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
 
VlnPlot(sce, c("CD3D", "CD4", "CD8A",
			   "LEF1", "CCR7", "TCF7",
			   "GZMK", "CD69", "AQP3",
			   "NKG7", "GNLY", "FCGR3A"), ncol = 3, pt.size = 0)

FeaturePlot(sce, "FCGR3A", label = T, order = T) + 
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

write.csv(sce.markers, "C1.TMarker.csv")  