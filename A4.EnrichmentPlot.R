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

new.cluster.ids <- c("Myeloid I", "T", "B", "NK", "Myeloid II", "Platelet")
names(new.cluster.ids) <- levels(immune.combined)
immune.combined <- RenameIdents(immune.combined, new.cluster.ids)
DefaultAssay(immune.combined) <- "RNA"

#-----GSEA-----
library(presto)
Sha_Clp <- subset(immune.combined,subset = groups != "3.HDC")
Clp_HDC <- subset(immune.combined,subset = groups != "1.Sha")
Sha_HDC <- subset(immune.combined,subset = groups != "2.Clp")

Sha_Clp.gene <- wilcoxauc(Sha_Clp,"groups")
Clp_HDC.gene <- wilcoxauc(Clp_HDC,"groups")
Sha_HDC.gene <- wilcoxauc(Sha_HDC,"groups")

setwd("/bigdata/wangzhang_guest/chenpeng_project/06_result/20_LJX_scRNAseq/03_csv")
write.csv(Sha_Clp.gene,"Sha_Clp.gene.csv")
write.csv(Clp_HDC.gene,"Clp_HDC.gene.csv")
write.csv(Sha_HDC.gene,"Sha_HDC.gene.csv")

Macrophage <- readRDS('/bigdata/wangzhang_guest/chenpeng_project/06_result/20_LJX_scRNAseq/01_rds/E4.Macrophage.rds')
DefaultAssay(Macrophage) <- "RNA"

library(presto)
Sha_ClpM2 <- subset(Macrophage,subset = groups != "3.HDC")
Clp_HDCM2 <- subset(Macrophage,subset = groups != "1.Sha")

Sha_ClpM2.gene <- wilcoxauc(Sha_ClpM2,"groups")
Clp_HDCM2.gene <- wilcoxauc(Clp_HDCM2,"groups")

setwd("/bigdata/wangzhang_guest/chenpeng_project/06_result/20_LJX_scRNAseq/03_csv")
write.csv(Sha_ClpM2.gene,"Sha_ClpM2.gene.csv")
write.csv(Clp_HDCM2.gene,"Clp_HDCM2.gene.csv")

#-----KEGG-----
immune.combined$celltype.stage <- paste(Idents(immune.combined), immune.combined$groups, sep = "_")
immune.combined$celltype <- Idents(immune.combined)
Idents(immune.combined) <- "celltype.stage"

M1.CS <- FindMarkers(immune.combined, ident.1 = "Myeloid I_2.Clp", ident.2 = "Myeloid I_1.Sha", verbose = FALSE)
M1.HC <- FindMarkers(immune.combined, ident.1 = "Myeloid I_3.HDC", ident.2 = "Myeloid I_2.Clp", verbose = FALSE)
M2.CS <- FindMarkers(immune.combined, ident.1 = "Myeloid II_2.Clp", ident.2 = "Myeloid II_1.Sha", verbose = FALSE)
M2.HC <- FindMarkers(immune.combined, ident.1 = "Myeloid II_3.HDC", ident.2 = "Myeloid II_2.Clp", verbose = FALSE)
P.CS <- FindMarkers(immune.combined, ident.1 = "Platelet_2.Clp", ident.2 = "Platelet_1.Sha", verbose = FALSE)
P.HC <- FindMarkers(immune.combined, ident.1 = "Platelet_3.HDC", ident.2 = "Platelet_2.Clp", verbose = FALSE)
T.CS <- FindMarkers(immune.combined, ident.1 = "T_2.Clp", ident.2 = "T_1.Sha", verbose = FALSE)
T.HC <- FindMarkers(immune.combined, ident.1 = "T_3.HDC", ident.2 = "T_2.Clp", verbose = FALSE)
NK.CS <- FindMarkers(immune.combined, ident.1 = "NK_2.Clp", ident.2 = "NK_1.Sha", verbose = FALSE)
NK.HC <- FindMarkers(immune.combined, ident.1 = "NK_3.HDC", ident.2 = "NK_2.Clp", verbose = FALSE)
B.CS <- FindMarkers(immune.combined, ident.1 = "B_2.Clp", ident.2 = "B_1.Sha", verbose = FALSE)
B.HC <- FindMarkers(immune.combined, ident.1 = "B_3.HDC", ident.2 = "B_2.Clp", verbose = FALSE)

M1.CS$Cluster = "0.MyeloidI"
M1.HC$Cluster = "0.MyeloidI"
M1.CS$gene = row.names(M1.CS)
M1.HC$gene = row.names(M1.HC)

M2.CS$Cluster = "4.MyeloidII"
M2.HC$Cluster = "4.MyeloidII"
M2.CS$gene = row.names(M2.CS)
M2.HC$gene = row.names(M2.HC)

#P.CS$Cluster = "5.Platelet"
#P.HC$Cluster = "5.Platelet"
#P.CS$gene = row.names(P.CS)
#P.HC$gene = row.names(P.HC)

T.CS$Cluster = "1.T"
T.HC$Cluster = "1.T"
T.CS$gene = row.names(T.CS)
T.HC$gene = row.names(T.HC)

NK.CS$Cluster = "3.NK"
NK.HC$Cluster = "3.NK"
NK.CS$gene = row.names(NK.CS)
NK.HC$gene = row.names(NK.HC)

B.CS$Cluster = "2.B"
B.HC$Cluster = "2.B"
B.CS$gene = row.names(B.CS)
B.HC$gene = row.names(B.HC)

all_CS <- bind_rows(M1.CS,M2.CS,T.CS,NK.CS,B.CS)
all_HC <- bind_rows(M1.HC,M2.HC,T.HC,NK.HC,B.HC)

library(Seurat)
library(gplots)
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)
Run_KEGG <- function(marker){
    #marker$gene <- row.names(marker)
    ids=bitr(marker$gene,'SYMBOL','ENTREZID','org.Mm.eg.db')
    sce.markers=merge(marker,ids,by.x='gene',by.y='SYMBOL')
    gcSample=split(sce.markers$ENTREZID, sce.markers$Cluster)
    result <- compareCluster(gcSample,
                             fun = "enrichKEGG",
                             organism = "mmu",
                             pvalueCutoff = 0.05)
    return (result)
}

CS_result <- Run_KEGG(all_CS)
HC_result <- Run_KEGG(all_HC)

setwd('/bigdata/wangzhang_guest/chenpeng_project/06_result/20_LJX_scRNAseq/03_csv')
write.csv(all_CS,"all_CS.csv")
write.csv(all_HC,"all_HC.csv")

write.csv(CS_result,"CS_result.csv")
write.csv(HC_result,"HC_result.csv")

P1 <- dotplot(CS_result)
P1 + theme(axis.text.x = element_text(
    angle = 45,
    vjust = 0.5, hjust = 0.5))

P2 <- dotplot(HC_result)
P2 + theme(axis.text.x = element_text(
    angle = 45,
    vjust = 0.5, hjust = 0.5))

setwd('/bigdata/wangzhang_guest/chenpeng_project/06_result/20_LJX_scRNAseq/02_plot')
pdf('P2.KEGG.pdf',12,8)
P1
P2
dev.off()

Run_GO <- function(marker){
    ids=bitr(marker$gene,'SYMBOL','ENTREZID','org.Mm.eg.db')
    sce.markers=merge(marker,ids,by.x='gene',by.y='SYMBOL')
    gcSample=split(sce.markers$ENTREZID, sce.markers$Cluster)
    display_number = c(18,18,18)
    result_CC <- compareCluster(gcSample,
                                fun = "enrichGO",
                                OrgDb = "org.Mm.eg.db",
                                ont = "CC",
                                pAdjustMethod = "BH",
                                pvalueCutoff = 0.01,
                                qvalueCutoff = 0.05)
    result_CC_show <- as.data.frame(result_CC)%>% group_by(Cluster) %>% do(head(., n = 3))
    result_MF <- compareCluster(gcSample,
                                fun = "enrichGO",
                                OrgDb = "org.Mm.eg.db",
                                ont = "MF",
                                pAdjustMethod = "BH",
                                pvalueCutoff = 0.01,
                                qvalueCutoff = 0.05)
    result_MF_show <- as.data.frame(result_MF)%>% group_by(Cluster) %>% do(head(., n = 3))
    result_BP <- compareCluster(gcSample,
                                fun = "enrichGO",
                                OrgDb = "org.Mm.eg.db",
                                ont = "BP",
                                pAdjustMethod = "BH",
                                pvalueCutoff = 0.01,
                                qvalueCutoff = 0.05)
    result_BP_show <- as.data.frame(result_BP)%>% group_by(Cluster) %>% do(head(., n = 3))
    num_BP <- dim(result_BP_show)[1]
    num_CC <- dim(result_CC_show)[1]
    num_MF <- dim(result_MF_show)[1]
    go_enrich_df <- data.frame(ID=c(result_BP_show$ID,
                                    result_CC_show$ID,
                                    result_MF_show$ID),
                               Description=c(result_BP_show$Description,
                                             result_CC_show$Description,
                                             result_MF_show$Description),
                               GeneNumber=c(result_BP_show$Count,
                                            result_CC_show$Count,
                                            result_MF_show$Count),
                               pvalue=c(result_BP_show$pvalue,
                                         result_CC_show$pvalue,
                                         result_MF_show$pvalue),
                               Cluster=c(paste(result_BP_show$Cluster),
                                         paste(result_CC_show$Cluster),
                                         paste(result_MF_show$Cluster)),
                               type=factor(c(rep("biological process", num_BP),
                                             rep("cellular component", num_CC),
                                             rep("molecular function", num_MF)),
                               levels=c("molecular function",
                                        "cellular component",
                                        "biological process")))
    return (go_enrich_df)
}

setwd('/bigdata/wangzhang_guest/chenpeng_project/06_result/20_LJX_scRNAseq/03_csv')
CS_GO_result <- Run_GO(all_CS)
write.csv(CAH_GO_result,"CS_GO_result.csv")
HC_GO_result <- Run_GO(all_HC)
write.csv(CAL_GO_result,"HC_GO_result.csv")

setwd('/bigdata/wangzhang_guest/chenpeng_project/06_result/20_LJX_scRNAseq/02_plot')
P9 <- ggplot(data=CS_GO_result, aes(x=Description, y=GeneNumber, fill=Cluster)) +
    geom_bar(stat="identity", width=0.8) + coord_flip() +
    #scale_fill_manual(values = cols) +
    facet_grid(type~.,scales = "free",space = "free") +
    coord_flip() +
    theme_bw() +
    xlab("GO term") +
    scale_x_discrete(labels=rev(CS_GO_result$Description)) +
    theme(axis.text=element_text(face = "bold", color="gray50")) +
    labs(title = "The Most Enriched GO Terms")

pdf('P2.GO_CS.pdf',12,12)
P9
dev.off()

P10 <- ggplot(data=HC_GO_result, aes(x=Description, y=GeneNumber, fill=Cluster)) +
    geom_bar(stat="identity", width=0.8) + coord_flip() +
    #scale_fill_manual(values = cols) +
    facet_grid(Cluster~.,scales = "free",space = "free") +
    coord_flip() +
    theme_bw() +
    xlab("GO term") +
    scale_x_discrete(labels=rev(HC_GO_result$Description)) +
    theme(axis.text=element_text(face = "bold", color="gray50")) +
    labs(title = "The Most Enriched GO Terms")

pdf('P2.GO_HC.pdf',12,12)
P10
dev.off()
