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

# load RDS 
#refquery <- readRDS(here('05_Integrate_df_annotation.RDS'))

# subset to Endothelial cells 
#refquery <- subset(refquery, cell_type_main == 'Epithelial cells')
#refquery <- ScaleData(object = refquery)
#refquery <- RunPCA(object = refquery)
#refquery <- FindNeighbors(refquery, dims = 1:20)
#refquery <- FindClusters(refquery, resolution = 0.6)
#refquery <- RunUMAP(object = refquery, dims = 1:20)

#epithelial_markers21 <- FindMarkers(refquery, ident.1 = 21, verbose = T)
#write.csv(epithelial_markers21, 'epithelial_markers21.csv', row.names = TRUE)
#rm(refquery)


# Immune cells
#refquery <- readRDS(here('05_Integrate_df_annotation.RDS'))

#refquery <- subset(refquery, cell_type_main == 'Immune cells')
#refquery <- ScaleData(object = refquery)
#refquery <- RunPCA(object = refquery)
#refquery <- FindNeighbors(refquery, dims = 1:20)
#refquery <- FindClusters(refquery, resolution = 0.6)
#refquery <- RunUMAP(object = refquery, dims = 1:20)

#Immune_markers20 <- FindMarkers(refquery, ident.1 = 20, verbose = T)
#write.csv(Immune_markers20, 'Immune_markers20.csv', row.names = TRUE)

#Immune_markers22 <- FindMarkers(refquery, ident.1 = 22, verbose = T)
#write.csv(Immune_markers22, 'Immune_markers22.csv', row.names = TRUE)

#rm(refquery)


# Endothelial cells
#refquery <- readRDS(here('05_Integrate_df_annotation.RDS'))

#refquery <- subset(refquery, cell_type_main == 'Endothelial cells')
#refquery <- ScaleData(object = refquery)
#refquery <- RunPCA(object = refquery)
#refquery <- FindNeighbors(refquery, dims = 1:20)
#refquery <- FindClusters(refquery, resolution = 0.6)
#refquery <- RunUMAP(object = refquery, dims = 1:20)

#Endothelial_markers8 <- FindMarkers(refquery, ident.1 = 8, verbose = T)
#write.csv(Endothelial_markers8, 'Endothelial_markers8.csv', row.names = TRUE)

#Endothelial_markers10 <- FindMarkers(refquery, ident.1 = 10, verbose = T)
#write.csv(Endothelial_markers10, 'Endothelial_markers10.csv', row.names = TRUE)

#Endothelial_markers12 <- FindMarkers(refquery, ident.1 = 12, verbose = T)
#write.csv(Endothelial_markers12, 'Endothelial_markers12.csv', row.names = TRUE)

#Endothelial_markers13 <- FindMarkers(refquery, ident.1 = 13, verbose = T)
#write.csv(Endothelial_markers13, 'Endothelial_markers13.csv', row.names = TRUE)

#Endothelial_markers14 <- FindMarkers(refquery, ident.1 = 14, verbose = T)
#write.csv(Endothelial_markers14, 'Endothelial_markers14.csv', row.names = TRUE)
#rm(refquery)


# Mesenchymal  cells
refquery <- readRDS(here('05_Integrate_df_annotation.RDS'))

refquery <- subset(refquery, cell_type_main == 'Mesenchymal  cells')
refquery <- ScaleData(object = refquery)
refquery <- RunPCA(object = refquery)
refquery <- FindNeighbors(refquery, dims = 1:20)
refquery <- FindClusters(refquery, resolution = 0.6)
refquery <- RunUMAP(object = refquery, dims = 1:20)

Mesenchymal_markers5 <- FindMarkers(refquery, ident.1 = 5, verbose = T)
write.csv(Mesenchymal_markers5, 'Mesenchymal_markers5_Findmarkers.csv', row.names = TRUE)

Mesenchymal_markers11 <- FindMarkers(refquery, ident.1 = 11, verbose = T)
write.csv(Mesenchymal_markers11, 'Mesenchymal_markers11_Findmarkers.csv', row.names = TRUE)

Mesenchymal_markers13 <- FindMarkers(refquery, ident.1 = 13, verbose = T)
write.csv(Mesenchymal_markers13, 'Mesenchymal_markers13_Findmarkers.csv', row.names = TRUE)
rm(refquery)
