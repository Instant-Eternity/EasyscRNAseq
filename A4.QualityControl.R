#!/usr/bin/env Rscript
rm(list = ls())
gc()

library(dplyr)
library(Seurat)
library(SoupX)
library(patchwork)

setwd("G://ZXL_221212//Data//SplitMatrix")

type = 'hsa'
if (type == "hsa"){
  Hb.gene.total <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
  Rb.gene <- "^RP[SL]"
  Mt.gene <- "^MT-"
} else if (type == "mmu") {
  Hb.gene.total <- c("Hba-a1","Hba-a2","Hbb-bt","Hbb-bs")
  Rb.gene <- "^Rp[sl]"
  Mt.gene <- "^mt-"
}

Pre_Count <- function(Seu_obj, Hb.gene.total, Rb.gene, Mt.gene){
  Hb <- match(Hb.gene.total, rownames(Seu_obj@assays$RNA))
  Hb.genes <- rownames(Seu_obj@assays$RNA)[Hb]
  Hb.genes <- Hb.genes[!is.na(Hb.genes)]
  Seu_obj[["percent.Hb"]] <- PercentageFeatureSet(Seu_obj, features = Hb.genes)
  Seu_obj[["percent.Rb"]] <- PercentageFeatureSet(Seu_obj, pattern = Rb.gene)
  Seu_obj[["percent.mt"]] <- PercentageFeatureSet(Seu_obj, pattern = Mt.gene)
  return (Seu_obj)
}

sample_info <- list.files(getwd())
scRNAlist <- list()
for(num in 1:length(sample_info)){
  doublets = read.csv(paste0(sample_info[num], "//doublet_prediction.txt"), header = T, row.names = 1)
  colnames(doublets) <- c("Doublet_score","Is_doublet")
  #row.names(doublets) <- paste0(num ,"-", row.names(doublets))
  sample <- Read10X(data.dir = paste0(sample_info[num])) %>%
    CreateSeuratObject(project = sample_info[num], min.cells = 3, min.features = 200) %>%
	AddMetaData(doublets) %>%
    RenameCells(add.cell.id = num) %>%
    Pre_Count(Hb.gene.total, Rb.gene, Mt.gene)
  scRNAlist[[num]] <- sample
}

setwd("G://ZXL_221212")
if(file.exists("A4.QualityControl") == TRUE){
  print("This directory already exists.")
}else{
  dir.create("A4.QualityControl")
}
setwd("A4.QualityControl")

names(scRNAlist) <- sample_info
system.time(saveRDS(scRNAlist, file = "R0.scRNAlist.rds"))

scRNA <- merge(scRNAlist[[1]], scRNAlist[2:length(scRNAlist)])

RowNum <- table(scRNA$orig.ident)
write.csv(as.data.frame(RowNum), file = "C1.RowCellNum.csv")

theme.set2 <- theme(axis.title.x = element_blank())
plot.featrures <- c("nFeature_RNA", "nCount_RNA", "percent.Hb", "percent.Rb", "percent.mt")
group <- "orig.ident"

plots <- list()
for(i in seq_along(plot.featrures)){
  plots[[i]] <- VlnPlot(scRNA, group.by=group, pt.size = 0,
                       features = plot.featrures[i]) + theme.set2 + NoLegend()}
violin <- wrap_plots(plots = plots, nrow=2)
ggsave("V1.VlnplotBeforeQC.pdf", plot = violin, height = 8, width = 9)

violin <- VlnPlot(scRNA, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.Rb"),ncol = 4,pt.size = 0.1) & 
  theme(plot.title = element_text(size=10))
ggsave("V2.VlnplotBeforeQC.pdf", plot = violin, height = 8, width = 12)

plot.featrures = data.frame(
  feature1 = c("nCount_RNA","nCount_RNA","nCount_RNA","percent.Rb","nFeature_RNA"),
  feature2 = c("percent.mt","nFeature_RNA","percent.Rb","percent.mt","Doublet_score")
)

plots <- list()
for(i in 1:dim(plot.featrures)[1]){
#print(i)
  plots[[i]] <- FeatureScatter(scRNA, feature1 = plot.featrures$feature1[i], 
                               feature2 = plot.featrures$feature2[i])
}

Feature <- wrap_plots(plots = plots)
ggsave("F1.FeatureScatterQC.pdf", plot = Feature, height = 8, width = 12)

#-----Doublet------
scRNA[['QC']] <- ifelse(scRNA@meta.data$Is_doublet == 'True','Doublet','Pass')
#-----HB-----
scRNA[['QC']] <- ifelse(scRNA@meta.data$percent.Hb > 0.01 &
                          scRNA@meta.data$QC == 'Pass',
                        'High_Hb',
                        scRNA@meta.data$QC)
#-----MT-----
scRNA[['QC']] <- ifelse(scRNA@meta.data$percent.mt > 15 & 
                          scRNA@meta.data$QC == 'Pass',
                        'High_MT',
                        scRNA@meta.data$QC)
scRNA[['QC']] <- ifelse(scRNA@meta.data$percent.mt > 15 & 
                          scRNA@meta.data$QC != 'Pass' &
                          scRNA@meta.data$QC != 'High_MT',
                        paste('High_MT',scRNA@meta.data$QC,sep = ' and '),
                        scRNA@meta.data$QC)
#-----RB-----
scRNA[['QC']] <- ifelse(scRNA@meta.data$percent.Rb > 50 & 
                          scRNA@meta.data$QC == 'Pass',
                        'High_Rb',
                        scRNA@meta.data$QC)
scRNA[['QC']] <- ifelse(scRNA@meta.data$percent.Rb > 50 & 
                          scRNA@meta.data$QC != 'Pass' &
                          scRNA@meta.data$QC != 'High_Rb',
                        paste('High_Rb',scRNA@meta.data$QC,sep = ' and '),
                        scRNA@meta.data$QC)
#-----Gene-----
scRNA[['QC']] <- ifelse(scRNA@meta.data$nFeature_RNA < 500 & 
                          scRNA@meta.data$QC == 'Pass',
                        'Low_nFeature',
                        scRNA@meta.data$QC)
scRNA[['QC']] <- ifelse(scRNA@meta.data$nFeature_RNA < 500 & 
                          scRNA@meta.data$QC != 'Pass' &
                          scRNA@meta.data$QC != 'Low_nFeature',
                        paste('Low_nFeature',scRNA@meta.data$QC,sep = ' and '),
                        scRNA@meta.data$QC)
scRNA[['QC']] <- ifelse(scRNA@meta.data$nFeature_RNA > 3000 & 
                          scRNA@meta.data$QC == 'Pass',
                        'High_nFeature',
                        scRNA@meta.data$QC)
scRNA[['QC']] <- ifelse(scRNA@meta.data$nFeature_RNA > 3000 & 
                          scRNA@meta.data$QC != 'Pass' &
                          scRNA@meta.data$QC != 'High_nFeature',
                        paste('High_nFeature',scRNA@meta.data$QC,sep = ' and '),
                        scRNA@meta.data$QC)
#-----UMI-----
scRNA[['QC']] <- ifelse(scRNA@meta.data$nCount_RNA > 15000 & 
                          scRNA@meta.data$QC == 'Pass',
                        'High_nCount',
                        scRNA@meta.data$QC)
scRNA[['QC']] <- ifelse(scRNA@meta.data$nCount_RNA > 15000 & 
                          scRNA@meta.data$QC != 'Pass' &
                          scRNA@meta.data$QC != 'High_nCount',
                        paste('High_nCount',scRNA@meta.data$QC,sep = ' and '),
                        scRNA@meta.data$QC)
table(scRNA[['QC']])

plot.featrures <- c("nFeature_RNA", "nCount_RNA", "percent.Hb", "percent.Rb", "percent.mt")
group <- "orig.ident"

plots <- list()
for(i in seq_along(plot.featrures)){
  plots[[i]] <- VlnPlot(subset(scRNA, subset = QC == 'Pass'), group.by=group, pt.size = 0,
                       features = plot.featrures[i]) + theme.set2 + NoLegend()}
violin <- wrap_plots(plots = plots, nrow=2)
ggsave("V1.VlnplotAfterQC.pdf", plot = violin, height = 8, width = 9)

violin <- VlnPlot(subset(scRNA, subset = QC == 'Pass'), 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.Rb"), ncol = 4, pt.size = 0.1) & 
  theme(plot.title = element_text(size=10))
ggsave("V2.VlnplotAfterQC.pdf", plot = violin, height = 8, width = 12)

QCResult <- as.data.frame(table(scRNA[['QC']]))
write.csv(QCResult, "C2.QCResult.csv")
saveRDS(scRNA, file = "R1.scRNA.rds")

scRNA <- subset(scRNA, QC == "Pass")
saveRDS(scRNA, file = "R2.SeuObjAfterQC.rds")	