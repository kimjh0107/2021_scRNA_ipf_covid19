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
refquery <- readRDS(here('05_Integrate_df_annotation.RDS'))

# subset to Mesenchymal cells 
refquery <- subset(refquery, cell_type_main == 'Mesenchymal cells')
refquery <- ScaleData(object = refquery)
refquery <- RunPCA(object = refquery)
refquery <- FindNeighbors(refquery, dims = 1:20)
refquery <- FindClusters(refquery, resolution = 0.6)
refquery <- RunUMAP(object = refquery, dims = 1:20)

# FindConservedMesenchymal_markers


# markers12
DefaultAssay(refquery) <- "RNA"
Mesenchymal_markers12 <- FindConservedMarkers(refquery, ident.1 = 12,  grouping.var = "Diagnosis", verbose = T)
head(Mesenchymal_markers12)
write.csv(Mesenchymal_markers12, 'Mesenchymal_markers12.csv', row.names = TRUE)

# markers13
DefaultAssay(refquery) <- "RNA"
Mesenchymal_markers13 <- FindConservedMarkers(refquery, ident.1 = 13,  grouping.var = "Diagnosis", verbose = T)
head(Mesenchymal_markers13)
write.csv(Mesenchymal_markers13, 'Mesenchymal_markers13.csv', row.names = TRUE)
