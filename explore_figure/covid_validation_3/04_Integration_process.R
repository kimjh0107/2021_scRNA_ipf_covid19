library(Matrix)
library(rhdf5)
library(tidyverse)
library(glue)
library(here)
library(hdf5r)
library(Seurat)
library(knitr)
library(sctransform)
library(Matrix)
library(future)
library(future.apply)



covid <- readRDS(here('03_covid_preprocess.RDS'))
ild <- readRDS(here('03_ipf_preprocess.RDS'))


# split the dataset
ipf.list <- SplitObject(ild, split.by = 'Diagnosis')
covid.list <- SplitObject(covid, split.by = 'Diagnosis')
object.list <- c(ipf.list, covid.list)

# Prepare object list 
object.list <- lapply(X = object.list, FUN = function(x) {
  x <- NormalizeData(x, verbose = T)
  x <- FindVariableFeatures(x, selection.method = "vst", verbose = T)
})

features <- SelectIntegrationFeatures(object.list = object.list)


object.list <- future_lapply(X = object.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = T)
  x <- RunPCA(x, features = features, verbose = T)
})


# Find anchors -> RPCA
anchors <- FindIntegrationAnchors(object.list = object.list, reduction = 'rpca', dims = 1:50)

# Integrate data sets
refquery <- IntegrateData(anchorset = anchors, dims = 1:50)

# Normal workflow
refquery <- ScaleData(object = refquery)
refquery <- FindVariableFeatures(object = refquery)
refquery <- RunPCA(object = refquery)
refquery <- FindNeighbors(refquery, dims = 1:30)
refquery <- FindClusters(refquery)
refquery <- RunUMAP(object = refquery, dims = 1:30)

saveRDS(refquery, file = "04_integration.RDS")



