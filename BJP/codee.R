rm(list = ls())
library(Seurat)
library(dplyr)
library(ggplot2)
setwd("/data/leili/kid/GSE131882/GSE131882_RAW")
load("/data/leili/kid/GSE131882/GSE131882_RAW/seurat.RData")
#####--------------------A.umap-------------------------------
colors <- c("grey","grey","grey",
            "grey","grey","grey","grey","grey",
            "grey","grey","#dc143c","grey")
DimPlot(seurat.combined, label = F, cols = colors,repel = T)+
  NoLegend()

#######------------------B.Dotplot-----------------------------------------
DefaultAssay(seurat.combined) <- "RNA"
genes_to_check = c('ITGA8','PDGFRB')
DotPlot(seurat.combined, group.by = 'seurat_clusters',cols = c("lightgrey","red"),
        features = unique(genes_to_check)) + RotatedAxis()

#####------------------------C.Dotplot2-----------------------------------
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
MC.cells <- subset(seurat.combined, idents = "10")
Idents(MC.cells) <- factor(Idents(MC.cells), levels = c("control", "diabetes"))
markers.to.plot <- c("YAP1","TAZ", "CCN2","CCN4","PDGFRB",
                     "LATS1","LATS2", "ACTA2")
DotPlot(MC.cells, features = markers.to.plot, group.by = 'group',
             cols = c("lightgrey","red")) + RotatedAxis()
