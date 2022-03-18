library(Seurat)
library(here)
library(tidyverse)

ipf_control <- readRDS(here('explore/RDS/02_ild_seurat_workflow.RDS'))
covid <- readRDS(here('explore/RDS/02_covid_seurat_workflow.RDS'))

# split the dataset
ipf.list <- SplitObject(ipf_control, split.by = 'Diagnosis')
object.list <- c(ipf.list, covid)

# Prepare object list 
object.list <- lapply(X = object.list, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, selection.method = "vst", verbose = FALSE)
})

features <- SelectIntegrationFeatures(object.list = object.list)


object.list <- future_lapply(X = object.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

# Find anchors -> RPCA
anchors <- FindIntegrationAnchors(object.list = object.list, reduction = 'rpca', dims = 1:50)

# Integrate data sets
refquery <- IntegrateData(anchorset = anchors, dims = 1:50)

# Normal workflow
refquery <- ScaleData(object = refquery)
refquery <- RunPCA(object = refquery)
refquery <- FindNeighbors(refquery, dims = 1:30)
refquery <- FindClusters(refquery)
refquery <- RunUMAP(object = refquery, dims = 1:30)

saveRDS(refquery, file = "03_Integrate_df.RDS")



