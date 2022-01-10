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

Clean_Soup <- function(path){
    sc <- load10X(path)
    sc <- autoEstCont(sc)
    out <- adjustCounts(sc)
    return (out)
}

Per_Count <- function(Seu_obj){
    Hb.genes_total <- c("Hba-a1","Hba-a2","Hbb-bt","Hbb-bs")
    Hb <- match(Hb.genes_total,rownames(Seu_obj@assays$RNA))
    Hb.genes <- rownames(Seu_obj@assays$RNA)[Hb]
    Hb.genes <- Hb.genes[!is.na(Hb.genes)]
    Seu_obj[["percent.Hb"]] <- PercentageFeatureSet(Seu_obj,features = Hb.genes)
    Rb.genes_total <- rownames(Seu_obj)[grep("^Rp[sl]",rownames(Seu_obj))]
    Rb <- match(Rb.genes_total,rownames(Seu_obj@assays$RNA))
    Rb.genes <- rownames(Seu_obj@assays$RNA)[Rb]
    Rb.genes <- Rb.genes[!is.na(Rb.genes)]
    Seu_obj[["percent.Rb"]] <- PercentageFeatureSet(Seu_obj,features = Rb.genes)
    Seu_obj[["percent.mt"]] <- PercentageFeatureSet(Seu_obj, pattern = "^mt-")
    return (Seu_obj)
}

Per_Clean <- function(Seu_obj,Feature_min,Featurn_max,Count_min,Count_max,Hb,Rb,mt){
    Seu_obj <- subset(Seu_obj,subset = nFeature_RNA > Feature_min &
                                       nFeature_RNA < Featurn_max &
                                       nCount_RNA >= Count_min &
                                       nCount_RNA <= Count_max &
                                       percent.Hb < Hb &
                                       percent.Rb < Rb &
                                       percent.mt < mt)
    return (Seu_obj)
}

Per_Process <- function(Seu_obj){
    Seu_obj <- NormalizeData(Seu_obj, normalization.method = "LogNormalize", scale.factor = 10000)
    Seu_obj <- FindVariableFeatures(Seu_obj, selection.method = "vst", nfeatures = 2000)
    Seu_obj <- ScaleData(Seu_obj, features = rownames(Seu_obj))
    Seu_obj <- RunPCA(Seu_obj, features = VariableFeatures(Seu_obj),npcs = 50)
    Seu_obj <- FindNeighbors(Seu_obj, dims = 1:20)
    Seu_obj <- FindClusters(Seu_obj, resolution = 0.5)
    Seu_obj <- RunUMAP(Seu_obj, dims = 1:20)
    Seu_obj <- RunTSNE(Seu_obj, dims = 1:20)
    return (Seu_obj)
}

Clean_Doublet <- function(Seu_obj){
    sweep.res.list_pbmc <- paramSweep_v3(Seu_obj, PCs = 1:30, sct = TRUE)
    sweep.stats_pbmc <- summarizeSweep(sweep.res.list_pbmc, GT = FALSE)
    bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
    opt_pK <- as.numeric(as.vector(bcmvn_pbmc$pK[which.max(bcmvn_pbmc$BCmetric)]))
    annotations <- Seu_obj@meta.data$seurat_clusters
    homotypic.prop <- modelHomotypic(annotations)
    nExp_poi <- round(0.06*nrow(Seu_obj@meta.data))
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
    Seu_obj <- doubletFinder_v3(Seu_obj,
                                PCs = 1:30,
                                pN = 0.25,
                                pK = opt_pK,
                                nExp = nExp_poi,
                                reuse.pANN = FALSE,
                                sct = TRUE)
    names(Seu_obj@meta.data)[12] <- 'DF.classifications'
    Seu_obj <- subset(Seu_obj, DF.classifications == 'Singlet')
}

Sha <- Clean_Soup('/bigdata/wangzhang_guest/chenpeng_project/01_data/21_LJX_scRNAseq/01_sha')
Clp <- Clean_Soup('/bigdata/wangzhang_guest/chenpeng_project/01_data/21_LJX_scRNAseq/02_clp')
HDC <- Clean_Soup('/bigdata/wangzhang_guest/chenpeng_project/01_data/21_LJX_scRNAseq/03_HDC')

Sha <- CreateSeuratObject(counts = Sha, project = "Sha", min.cells = 3, min.features = 200)
Clp <- CreateSeuratObject(counts = Clp, project = "Clp", min.cells = 3, min.features = 200)
HDC <- CreateSeuratObject(counts = HDC, project = "HDC", min.cells = 3, min.features = 200)

Sha <- Per_Count(Sha)
Clp <- Per_Count(Clp)
HDC <- Per_Count(HDC)

V1 <- VlnPlot(Sha, features = c("nFeature_RNA", "nCount_RNA", "percent.Hb", "percent.Rb", "percent.mt"))
V2 <- VlnPlot(Clp, features = c("nFeature_RNA", "nCount_RNA", "percent.Hb", "percent.Rb", "percent.mt"))
V3 <- VlnPlot(HDC, features = c("nFeature_RNA", "nCount_RNA", "percent.Hb", "percent.Rb", "percent.mt"))

pdf("/bigdata/wangzhang_guest/chenpeng_project/06_result/20_LJX_scRNAseq/02_plot/E1.VlnPlot.pdf")
V1
V2
V3
dev.off()

Sha <- Per_Clean(Sha,200,4000,0,20000,0.01,50,10)
Clp <- Per_Clean(Clp,200,4000,0,20000,0.01,50,10)
HDC <- Per_Clean(HDC,200,4000,0,20000,0.01,50,10)

Sha <- Per_Process(Sha)
Clp <- Per_Process(Clp)
HDC <- Per_Process(HDC)

Sha$groups <- '1.Sha'
Clp$groups <- '2.Clp'
HDC$groups <- '3.HDC'

Sha$sample <- 'Sha.1'
Clp$sample <- 'Clp.1'
HDC$sample <- 'HDC.1'

Sha <- Clean_Doublet(Sha)
Clp <- Clean_Doublet(Clp)
HDC <- Clean_Doublet(HDC)

info.list <- c(Sha,Clp,HDC)
features <- SelectIntegrationFeatures(object.list = info.list)
immune.anchors <- FindIntegrationAnchors(object.list = info.list, anchor.features = features)
immune.combined <- IntegrateData(anchorset = immune.anchors)

saveRDS(immune.combined, "/bigdata/wangzhang_guest/chenpeng_project/06_result/20_LJX_scRNAseq/01_rds/A1.immune.combined")

library(harmony)

DefaultAssay(immune.combined) <- "integrated"

immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
sce.all.filt=immune.combined
sce.all.filt@meta.data <- sce.all.filt@meta.data[,1:10]

DefaultAssay(immune.combined) <- "RNA"
#object <- RunHarmony(object = object,group.by.vars = "var",assay.use = "SCT")
object <- RunHarmony(object = sce.all.filt,group.by.vars = "groups",assay.use = "integrated")

sce.all.int <- object
sce.all.int <- RunHarmony(sce.all.filt,c( "orig.ident" ))

names(sce.all.int@reductions)
harmony_embeddings <- Embeddings(sce.all.int, 'harmony')
harmony_embeddings[1:5, 1:5]

sce.all.int=RunTSNE(sce.all.int,reduction = "harmony", dims = 1:30)
sce.all.int=RunUMAP(sce.all.int,reduction = "harmony",dims = 1:30)

sce=sce.all.int
sce <- FindNeighbors(sce, reduction = "harmony",dims = 1:15)
sce <- FindClusters(sce, resolution = 0.015)

table(sce$orig.ident,sce$seurat_clusters)

p1=DimPlot(sce,reduction = "umap",label=T
           ,group.by = 'groups')

p2=DimPlot(sce,reduction = "umap",label=T
           ,group.by = 'orig.ident')

p3=DimPlot(sce,reduction = "umap",label=T)

library(patchwork)

pdf('/bigdata/wangzhang_guest/chenpeng_project/06_result/20_LJX_scRNAseq/02_plot/E2.DimPlot.pdf',width = 9,height = 9)
p1
p2
p3
dev.off()

immune.combined <- sce

DefaultAssay(immune.combined) <- "RNA"

immune.combined.markers <- FindAllMarkers(immune.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top2.marker <- immune.combined.markers %>% group_by(cluster) %>% slice_max(n = 2, order_by = avg_log2FC)

write.csv(immune.combined.markers,"/bigdata/wangzhang_guest/chenpeng_project/06_result/20_LJX_scRNAseq/03_csv")
saveRDS(immune.combined, "/bigdata/wangzhang_guest/chenpeng_project/06_result/20_LJX_scRNAseq/01_rds/E1.harmony.rds")
