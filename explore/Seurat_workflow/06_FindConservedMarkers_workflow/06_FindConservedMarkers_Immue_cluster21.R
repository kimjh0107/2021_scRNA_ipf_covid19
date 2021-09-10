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

# subset to Immune cells 
refquery <- subset(refquery, cell_type_main == 'Immune cells')
refquery <- ScaleData(object = refquery)
refquery <- RunPCA(object = refquery)
refquery <- FindNeighbors(refquery, dims = 1:20)
refquery <- FindClusters(refquery, resolution = 0.6)
refquery <- RunUMAP(object = refquery, dims = 1:20)

# FindConservedImmune_markers
# markers21
DefaultAssay(refquery) <- "RNA"
Immune_markers21 <- FindConservedMarkers(refquery, ident.1 = 21,  grouping.var = "Diagnosis", verbose = T)
head(Immune_markers21)
write.csv(Immune_markers21, 'Immune_markers21.csv', row.names = TRUE)

# markers22
DefaultAssay(refquery) <- "RNA"
Immune_markers21 <- FindConservedMarkers(refquery, ident.1 = 22,  grouping.var = "Diagnosis", verbose = T)
head(Immune_markers22)
write.csv(Immune_markers22, 'Immune_markers22.csv', row.names = TRUE)
