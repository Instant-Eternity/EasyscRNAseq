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

quickmerge <- function(df1, df2) {
    df1.names <- names(df1)
    df2.names <- names(df2)
    df2.add <- setdiff(df1.names, df2.names)
    df1.add <- setdiff(df2.names, df1.names)
    if(length(df2.add) > 0) {
        for(i in 1:length(df2.add)) {
            df2[df2.add[i]] <- NA
        }
    }
    if(length(df1.add) > 0) {
        for(i in 1:length(df1.add)) {
            df1[df1.add[i]] <- NA
        }
    }
    return(rbind(df1, df2))
}

site <- immune.combined@reductions$umap@cell.embeddings
info <- immune.combined@meta.data[,c(1,8,9)]

row.data <- as.data.frame(t(quickmerge(t(info),t(site))))

setwd("/bigdata/wangzhang_guest/chenpeng_project/06_result/20_LJX_scRNAseq/03_csv")
write.csv(row.data,"P1.A.row_data.csv")

gene_expr <- FetchData(object = immune.combined, vars = c("Adgre1","S100a4","Ms4a6c"))
write.csv(gene_expr,"P1.C.exprM.csv")

gene_expr <- FetchData(object = immune.combined, vars = c("groups","sample","seurat_clusters","Cd3d","Cd3e","Cd3g","Ncr1","Klrb1c","Cd19","Cd79a","Cd79b","Itgam","Ppbp","Itga2b"))
write.csv(gene_expr,"P1.C.expr.csv")
