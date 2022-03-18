library(Seurat)
library(here)
library(tidyverse)

# IPF Preprocessing workflow
ild <- readRDS(here('explore/01_IPF_SeuratObject.RDS'))

ild <- NormalizeData(ild)
ild <- FindVariableFeatures(ild)
ild <- ScaleData(ild)
ild <- RunPCA(ild)
ild <- FindNeighbors(ild, dims = 1:30)
ild <- FindClusters(ild)
ild <- RunUMAP(ild, dims = 1:30)

# COVID Preprocessing workflow 
covid <- readRDS(here('explore/01_COVID_h5_to_SeuratObject_merge.RDS'))

covid <- NormalizeData(covid)
covid <- FindVariableFeatures(covid)
covid <- ScaleData(covid)
covid <- RunPCA(covid)
covid <- FindNeighbors(covid, dims = 1:30)
covid <- FindClusters(covid)
covid <- RunUMAP(covid, dims = 1:30)





