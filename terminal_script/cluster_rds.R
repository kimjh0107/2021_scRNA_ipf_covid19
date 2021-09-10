library(Seurat)
library(here)
library(tidyverse)
library(cowplot)
library(knitr)

library(BiocManager)
library(sctransform)
library(glmGamPoi)
library(Matrix)
library(future)
library(future.apply)

refquery <- readRDS(here('05_Integrate_df_annotation.RDS'))


# Epithelial cells cluster plot 

Epithelial <- subset(refquery, cell_type_main == 'Epithelial cells')
Epithelial <- ScaleData(object = Epithelial)
Epithelial <- RunPCA(object = Epithelial)
Epithelial <- FindNeighbors(Epithelial, dims = 1:20)
Epithelial <- FindClusters(Epithelial, resolution = 0.6)
Epithelial <- RunUMAP(object = Epithelial, dims = 1:20)

saveRDS(Epithelial, file = "Epithelial_cluster.RDS")


# Immune cells

Immune <- subset(refquery, cell_type_main == 'Immune cells')
Immune <- ScaleData(object = Immune)
Immune <- RunPCA(object = Immune)
Immune <- FindNeighbors(Immune, dims = 1:20)
Immune <- FindClusters(Immune, resolution = 0.6)
Immune <- RunUMAP(object = Immune, dims = 1:20)
saveRDS(Immune, file = "Immune_cluster.RDS")


# Endothelial cells

Endothelial <- subset(refquery, cell_type_main == 'Endothelial cells')
Endothelial <- ScaleData(object = Endothelial)
Endothelial <- RunPCA(object = Endothelial)
Endothelial <- FindNeighbors(Endothelial, dims = 1:20)
Endothelial <- FindClusters(Endothelial, resolution = 0.6)
Endothelial <- RunUMAP(object = Endothelial, dims = 1:20)
saveRDS(Endothelial, file = "Endothelial_cluster.RDS")



# Mesenchymal  cells
Mesenchymal <- subset(refquery, cell_type_main == 'Mesenchymal cells')
Mesenchymal <- ScaleData(object = Mesenchymal)
Mesenchymal <- RunPCA(object = Mesenchymal)
Mesenchymal <- FindNeighbors(Mesenchymal, dims = 1:20)
Mesenchymal <- FindClusters(Mesenchymal, resolution = 0.6)
Mesenchymal <- RunUMAP(object = Mesenchymal, dims = 1:20)
saveRDS(Mesenchymal, file = "Mesenchymal_cluster.RDS")


