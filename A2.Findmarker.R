#!/usr/bin/env Rscript
rm(list = ls())
gc()
library(dplyr)
library(scater)
library(patchwork)
library(ggplot2)
library(RColorBrewer)
library(corrplot)
library(cowplot)
library(Seurat)
library(SoupX)
library(Matrix)
library(ggsci)
library(DoubletFinder)

set.seed(6565)


library(patchwork)

immune.combined <- readRDS("/bigdata/wangzhang_guest/chenpeng_project/06_result/20_LJX_scRNAseq/01_rds/E1.harmony.rds")

DefaultAssay(immune.combined) <- "RNA"

immune.combined.markers <- FindAllMarkers(immune.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top2.marker <- immune.combined.markers %>% group_by(cluster) %>% slice_max(n = 2, order_by = avg_log2FC)
