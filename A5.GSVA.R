#!/usr/bin/env Rscript
rm(list = ls())
gc()

library(Matrix)
library(dplyr)
library(Seurat)
library(scater)
library(patchwork)
library(ggplot2)
library(ggsci)
library(RColorBrewer)
library(monocle)
library(corrplot)
library(DoubletFinder)
library(ggnewscale)
set.seed(6565)

immune.combined <- readRDS("/bigdata/wangzhang_guest/chenpeng_project/06_result/20_LJX_scRNAseq/01_rds/E1.harmony.rds")

DefaultAssay(immune.combined) <- "RNA"

expr <- AverageExpression(immune.combined, group.by = "sample", assays = "RNA", slot = "data")[[1]]
expr <- expr[rowSums(expr)>0,]
expr <- as.matrix(expr)
head(expr)

expr0 <- AverageExpression(subset(immune.combined,subset = seurat_clusters == 0), group.by = "sample", assays = "RNA", slot = "data")[[1]]
expr0 <- expr0[rowSums(expr0)>0,]
expr0 <- as.data.frame(expr0)
expr0$cluster = '0'
expr0$gene = row.names(expr0)
expr1 <- AverageExpression(subset(immune.combined,subset = seurat_clusters == 1), group.by = "sample", assays = "RNA", slot = "data")[[1]]
expr1 <- expr1[rowSums(expr1)>0,]
expr1 <- as.data.frame(expr1)
expr1$cluster = '1'
expr1$gene = row.names(expr1)
expr2 <- AverageExpression(subset(immune.combined,subset = seurat_clusters == 2), group.by = "sample", assays = "RNA", slot = "data")[[1]]
expr2 <- expr2[rowSums(expr2)>0,]
expr2 <- as.data.frame(expr2)
expr2$cluster = '2'
expr2$gene = row.names(expr2)
expr3 <- AverageExpression(subset(immune.combined,subset = seurat_clusters == 3), group.by = "sample", assays = "RNA", slot = "data")[[1]]
expr3 <- expr3[rowSums(expr3)>0,]
expr3 <- as.data.frame(expr3)
expr3$cluster = '3'
expr3$gene = row.names(expr3)
expr4 <- AverageExpression(subset(immune.combined,subset = seurat_clusters == 4), group.by = "sample", assays = "RNA", slot = "data")[[1]]
expr4 <- expr4[rowSums(expr4)>0,]
expr4 <- as.data.frame(expr4)
expr4$cluster = '4'
expr4$gene = row.names(expr4)
all_expr <- bind_rows(expr0,expr1,expr2,expr3,expr4)
all_expr <- all_expr[,c(3,1,2,4,5)]
all_expr <- as.matrix(all_expr)

setwd("/bigdata/wangzhang_guest/chenpeng_project/06_result/20_LJX_scRNAseq/03_csv")
write.csv(all_expr,'P2.all_expr.csv')

setwd("/bigdata/wangzhang_guest/chenpeng_project/06_result/20_LJX_scRNAseq/03_csv")
write.csv(expr,'P2.expr.csv')

Macrophage <- readRDS("/bigdata/wangzhang_guest/chenpeng_project/06_result/20_LJX_scRNAseq/01_rds/E4.Macrophage.rds")
DefaultAssay(Macrophage) <- "RNA"

expr <- AverageExpression(Macrophage, group.by = "sample", assays = "RNA", slot = "data")[[1]]
expr <- expr[rowSums(expr)>0,]
expr <- as.matrix(expr)
head(expr)

setwd("/bigdata/wangzhang_guest/chenpeng_project/06_result/20_LJX_scRNAseq/03_csv")
write.csv(expr,'P4.exprMacrophage.csv')
