library(Matrix)
library(rhdf5)
library(tidyverse)
library(glue)
library(here)
library(hdf5r)
library(Seurat)

covid_val <- readRDS(here('COVID_validate_h5_to_SeuratObject_merge.RDS'))

# QC 
covid_val <- PercentageFeatureSet(covid_val, pattern = "^MT-", col.name = "percent.mt")
covid_val <- subset(covid_val, subset = nFeature_RNA > 1000 & percent.mt < 25)

# Preprocessing
covid_val <- NormalizeData(covid_val, normalization.method = "LogNormalize", scale.factor = 10000)
covid_val <- FindVariableFeatures(covid_val, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(covid_val)
covid_val <- ScaleData(covid_val, features = all.genes)
covid_val <- RunPCA(covid_val)
covid_val <- FindNeighbors(covid_val, dims = 1:30, verbose = FALSE)
covid_val <- FindClusters(covid_val, verbose = FALSE)

covid_val <- RunUMAP(covid_val, dims = 1:30, verbose = FALSE, return.model=TRUE)


saveRDS(covid_val, file = "COVID_val_Preprocessing2.RDS")







