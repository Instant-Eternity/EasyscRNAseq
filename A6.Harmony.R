#!/usr/bin/env Rscript
rm(list = ls())
gc()

library(dplyr)
library(Seurat)
library(ggplot2)
library(harmony)

setwd("G://ZXL_221212")

scRNA <- readRDS(file = "A5.Normalization//R4.SeuObjBeforeHarmony.rds")

if(file.exists("A6.Harmony") == TRUE){
  print("This directory already exists.")
}else{
  dir.create("A6.Harmony")
}
setwd("A6.Harmony")

scRNA <- RunHarmony(scRNA, group.by.vars="orig.ident", assay.use="RNA", max.iter.harmony = 20, lambda = 1.5) 
harmony_embeddings <- Embeddings(scRNA, 'harmony')
harmony_embeddings[1:5, 1:5]

library(cowplot)
p1 <- DimPlot(object = scRNA, reduction = "harmony", pt.size = .1, group.by = "orig.ident") + NoLegend()
p2 <- VlnPlot(object = scRNA, features = "harmony_1", group.by = "orig.ident", pt.size = .1) + NoLegend()
ggsave("G1.AfterHarmony.pdf", plot = plot_grid(p1,p2), height = 6, width = 12)

scRNA <- scRNA %>% 
  RunUMAP(reduction = "harmony", dims = 1:30, verbose = F) %>% 
  FindNeighbors(reduction = "harmony", k.param = 10, dims = 1:30) %>% 
  FindClusters() %>% 
  identity()

library(patchwork)
scRNA <- SetIdent(scRNA,value = "orig.ident")
plot1 <- DimPlot(scRNA,reduction = "umap") + 
  plot_annotation(title = "PBMC cells, after integration (Harmony)")
ggsave("D1.UMAPDimAfterHar.pdf", plot = plot1, height = 12, width = 12)

plot2 <- DimPlot(scRNA, reduction = "umap", group.by = "orig.ident", pt.size = .1, split.by = 'orig.ident', ncol = 3) + NoLegend()
ggsave("D2.UMAPDimAfterHarSample.pdf", plot = plot2, height = 12, width = 12)

scRNA <- SetIdent(scRNA,value = "seurat_clusters")
plot3 <- DimPlot(scRNA,label = T) + NoLegend()
ggsave("D3.UMAPDimIdents.pdf", plot = plot3, height = 12, width = 12)

#saveRDS(scRNA, file = "R6.SeuObjAfterHarmony.rds")
saveRDS(scRNA, file = "R5.SeuObjAfterHarmony.rds")

immune.combined <- readRDS(file = "R5.SeuObjAfterHarmony.rds")

library(stringr)
immune.combined@meta.data$group <- str_sub(immune.combined@meta.data$orig.ident, start = 1, end = -3)

ElbowPlot(immune.combined, ndims = 50, reduction = "harmony")

Seed = 6000 #2000
immune.combined <- immune.combined %>% 
  RunUMAP(reduction = "harmony", dims = 1:20, verbose = F, seed.use = Seed) %>% 
  FindNeighbors(reduction = "harmony", k.param = 10, dims = 1:20) %>% 
  FindClusters(resolution = 0.02)
  
DimPlot(immune.combined, reduction = "umap", seed = Seed, label = T) 

immune.combined.markers <- FindAllMarkers(immune.combined, 
                                          only.pos = TRUE, 
                                          min.pct = 0.25, 
                                          logfc.threshold = 0.25) 

immune.combined.markers %>%
    group_by(cluster) %>%
    slice_max(n = 5, order_by = avg_log2FC) %>%
	print(n = 25)

immune.combined.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(immune.combined, features = top10$gene) + NoLegend()
ggsave("H1.DoHeatmap.pdf", height = 12, width = 18)

#cluster3.markers <- FindMarkers(immune.combined, ident.1 = 3, ident.2 = c(0, 1), min.pct = 0.25)

write.csv(immune.combined.markers, "Har.allMarker.csv")
#write.csv(cluster3.markers, "Cluster3.marker.csv")

plot4 <- DimPlot(immune.combined, reduction = "umap", seed = Seed, label = T) 
ggsave("D4.UMAPDimRe0.5.pdf", plot = plot4, height = 12, width = 12)
plot5 <- DimPlot(immune.combined, reduction = "umap", seed = Seed, label = T, split.by = "group") 
ggsave("D5.UMAPDimRe0.5group.pdf", plot = plot5, height = 9, width = 27)

cellNum <- as.data.frame(table(immune.combined@meta.data$group,immune.combined@meta.data$seurat_clusters))
library(tidyr)
library(reshape2)
cellNum <- spread(cellNum, Var1, Freq)
row.names(cellNum) <- cellNum[,1]
cellNum <- cellNum[,-1]
write.csv(cellNum, "Har.cellNum.csv")
library(qiime2R)
cellPercent <- as.data.frame(make_percent(cellNum))
write.csv(cellPercent, "Har.cellPercent.csv")
sampleNum <- as.data.frame(table(immune.combined@meta.data$orig.ident,immune.combined@meta.data$seurat_clusters))
sampleNum <- spread(sampleNum, Var1, Freq)
row.names(sampleNum) <- sampleNum[,1]
sampleNum <- sampleNum[,-1]
write.csv(sampleNum, "Har.sampleNum.csv")
samplePercent <- as.data.frame(make_percent(sampleNum))
write.csv(samplePercent, "Har.samplePercent.csv")

library(RColorBrewer)
plots <- list()
plot.featrures <- c("CD3D", "NKG7", "FCGR3A", "CD14", "CD79A", "MS4A1", "PPBP", "PF4", "PASK")
for(i in seq_along(plot.featrures)){
  plots[[i]] <- FeaturePlot(immune.combined, plot.featrures[i], label = F, order = T) + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
}
Feature <- wrap_plots(plots = plots, ncol = 3)
ggsave("F1.FeaturePlotCluster.pdf", plot = Feature, height = 15, width = 16)
			 
plots <- list()
plot.featrures <- c("CD3D", "NKG7", "FCGR3A", "CD14", "CD79A", "MS4A1", "PPBP", "PF4", "PASK")
col.data <- c("#DE5E60", "#5F9BD0", "#86BF38", "#7B2A9D", "#6D818E")
for(i in seq_along(plot.featrures)){
  plots[[i]] <- VlnPlot(immune.combined, plot.featrures[i], pt.size = 0) +
  scale_fill_manual(values = col.data ) # 
}
Vln <- wrap_plots(plots = plots, ncol = 1)
ggsave("V1.VlnPlotCluster.pdf", plot = Vln, height = 15, width = 3)

new.cluster.ids <- c("T/NK cell", "Myeloid cell", "B cell", "Platelet", "PASK+ cell")
names(new.cluster.ids) <- levels(immune.combined)
immune.combined <- RenameIdents(immune.combined, new.cluster.ids)

saveRDS(immune.combined, file = "R6.SeuObjAfterCluster.rds")