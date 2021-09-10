library(Seurat)
library(here)
library(tidyverse)
library(cowplot)
library(knitr)
library(sctransform)
library(Matrix)
library(future)
library(future.apply)

# workers
plan('multiprocess', workers = 20)
options(future.globals.maxSize = 2e+05 * 1024^2)

# Epithelial cells
refquery <- readRDS(here('05_Integrate_df_annotation.RDS'))

# subset to Epithelial cells 
refquery <- subset(refquery, cell_type_main == 'Epithelial cells')
refquery <- ScaleData(object = refquery)
refquery <- RunPCA(object = refquery)
refquery <- FindNeighbors(refquery, dims = 1:20)
refquery <- FindClusters(refquery, resolution = 0.6)
refquery <- RunUMAP(object = refquery, dims = 1:20)

DimPlot(refquery, reduction = "umap", group.by = "Diagnosis", label = TRUE, label.size = 3, repel = TRUE)
ggsave(here('figure/FindConservedMarkers_Diagnosis_proportion/Epithelial_Diagnosis.pdf'))

rm(refquery)

# Immune cells
refquery <- readRDS(here('05_Integrate_df_annotation.RDS'))

refquery <- subset(refquery, cell_type_main == 'Immune cells')
refquery <- ScaleData(object = refquery)
refquery <- RunPCA(object = refquery)
refquery <- FindNeighbors(refquery, dims = 1:20)
refquery <- FindClusters(refquery, resolution = 0.6)
refquery <- RunUMAP(object = refquery, dims = 1:20)

DimPlot(refquery, reduction = "umap", group.by = "Diagnosis", label = TRUE, label.size = 3, repel = TRUE)
ggsave(here('figure/FindConservedMarkers_Diagnosis_proportion/Immune_Diagnosis.pdf'))

rm(refquery)


# Endothelial cells
refquery <- readRDS(here('05_Integrate_df_annotation.RDS'))

refquery <- subset(refquery, cell_type_main == 'Endothelial cells')
refquery <- ScaleData(object = refquery)
refquery <- RunPCA(object = refquery)
refquery <- FindNeighbors(refquery, dims = 1:20)
refquery <- FindClusters(refquery, resolution = 0.6)
refquery <- RunUMAP(object = refquery, dims = 1:20)

DimPlot(refquery, reduction = "umap", group.by = "Diagnosis", label = TRUE, label.size = 3, repel = TRUE)
ggsave(here('figure/FindConservedMarkers_Diagnosis_proportion/Endothelial_Diagnosis.pdf'))

rm(refquery)

# Mesenchymal  cells
refquery <- readRDS(here('05_Integrate_df_annotation.RDS'))

refquery <- subset(refquery, cell_type_main == 'Mesenchymal cells')
refquery <- ScaleData(object = refquery)
refquery <- RunPCA(object = refquery)
refquery <- FindNeighbors(refquery, dims = 1:20)
refquery <- FindClusters(refquery, resolution = 0.6)
refquery <- RunUMAP(object = refquery, dims = 1:20)

DimPlot(refquery, reduction = "umap", group.by = "Diagnosis", label = TRUE, label.size = 3, repel = TRUE)
ggsave(here('figure/FindConservedMarkers_Diagnosis_proportion/Mesenchymal_Diagnosis.pdf'))

rm(refquery)