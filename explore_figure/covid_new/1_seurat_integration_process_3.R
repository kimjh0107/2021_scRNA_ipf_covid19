# Seurat tutorial의 정석적인 방법
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



covid_val <- readRDS(here('01_COVID_merge_h5_to_SeuratObject_merge.RDS'))
ipf_val <- readRDS(here('01_IPF_SeuratObject.RDS'))

# IPF preprocessing 
ipf_val <- PercentageFeatureSet(ipf_val, pattern = "^MT-", col.name = "percent.mt")
ipf_val <- subset(ipf_val, subset = nFeature_RNA > 1000 & percent.mt < 25)

# Covid preprocessing
covid_val <- PercentageFeatureSet(object = covid_val, pattern = "^MT-", col.name = "percent.mt")
covid_val <- subset(covid_val, subset = nFeature_RNA > 1000 & percent.mt < 25)

# split the dataset
#ipf.list <- SplitObject(ipf_val, split.by = 'Diagnosis')
#covid.list <- SplitObject(covid_val, split.by = 'Diagnosis')
object.list <- c(ipf_val, covid_val)


# Prepare object list 
object.list <- lapply(X = object.list, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, selection.method = "vst", verbose = FALSE)
})


features <- SelectIntegrationFeatures(object.list = object.list)



# Find anchors -> RPCA
anchors <- FindIntegrationAnchors(object.list = object.list, anchor.features = features, dims = 1:50)

# Integrate data sets
refquery <- IntegrateData(anchorset = anchors, dims = 1:50)

# Normal workflow
DefaultAssay(refquery) <- "integrated"
refquery <- ScaleData(object = refquery)
refquery <- RunPCA(object = refquery)
refquery <- RunUMAP(object = refquery, dims = 1:30)
refquery <- FindNeighbors(refquery, dims = 1:30)
refquery <- FindClusters(refquery)

saveRDS(refquery, file = "Integrate_df_seuratworkflow.RDS")