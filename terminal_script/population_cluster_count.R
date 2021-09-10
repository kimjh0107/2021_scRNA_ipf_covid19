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

plan('multiprocess', workers = 20)
options(future.globals.maxSize = 2e+05 * 1024^2)

# Epithelial cells
refquery <- readRDS(here('05_Integrate_df_annotation.RDS'))

refquery <- subset(refquery, cell_type_main == 'Epithelial cells')
refquery <- ScaleData(object = refquery)
refquery <- RunPCA(object = refquery)
refquery <- FindNeighbors(refquery, dims = 1:20)
refquery <- FindClusters(refquery, resolution = 0.6)
refquery <- RunUMAP(object = refquery, dims = 1:20)

unique(refquery@meta.data$integrated_snn_res.0.6)
unique(refquery@meta.data$integrated_snn_res.0.8)
unique(refquery@meta.data$seurat_clusters)

rm(refquery)


# Immune cells
refquery <- readRDS(here('05_Integrate_df_annotation.RDS'))

refquery <- subset(refquery, cell_type_main == 'Immune cells')
refquery <- ScaleData(object = refquery)
refquery <- RunPCA(object = refquery)
refquery <- FindNeighbors(refquery, dims = 1:20)
refquery <- FindClusters(refquery, resolution = 0.6)
refquery <- RunUMAP(object = refquery, dims = 1:20)

unique(refquery@meta.data$integrated_snn_res.0.6)
unique(refquery@meta.data$integrated_snn_res.0.8)
unique(refquery@meta.data$seurat_clusters)

rm(refquery)


# Endothelial cells 
refquery <- readRDS(here('05_Integrate_df_annotation.RDS'))

refquery <- subset(refquery, cell_type_main == 'Endothelial cells')
refquery <- ScaleData(object = refquery)
refquery <- RunPCA(object = refquery)
refquery <- FindNeighbors(refquery, dims = 1:20)
refquery <- FindClusters(refquery, resolution = 0.6)
refquery <- RunUMAP(object = refquery, dims = 1:20)

unique(refquery@meta.data$integrated_snn_res.0.6)
unique(refquery@meta.data$integrated_snn_res.0.8)
unique(refquery@meta.data$seurat_clusters)

rm(refquery)


# Mesenchymal cells
refquery <- readRDS(here('05_Integrate_df_annotation.RDS'))

refquery <- subset(refquery, cell_type_main == 'Mesenchymal cells')
refquery <- ScaleData(object = refquery)
refquery <- RunPCA(object = refquery)
refquery <- FindNeighbors(refquery, dims = 1:20)
refquery <- FindClusters(refquery, resolution = 0.6)
refquery <- RunUMAP(object = refquery, dims = 1:20)

unique(refquery@meta.data$integrated_snn_res.0.6)
unique(refquery@meta.data$integrated_snn_res.0.8)
unique(refquery@meta.data$seurat_clusters)

rm(refquery)