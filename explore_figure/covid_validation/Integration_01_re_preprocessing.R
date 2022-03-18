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

plan('multiprocess', workers = 40)
options(future.globals.maxSize = 2e+05 * 1024^2)

# read RDS
covid_val <- readRDS(here('COVID_validate_h5_to_SeuratObject_merge.RDS'))
ipf_val <- readRDS(here('01_IPF_SeuratObject.RDS'))


# IPF preprocessing 
ipf_val <- PercentageFeatureSet(ipf_val, pattern = "^MT-", col.name = "percent.mt")
ipf_val <- subset(ipf_val, subset = nFeature_RNA > 1000 & percent.mt < 25)
ipf_val <- SCTransform(ipf_val, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = T)
ipf_val <- RunPCA(ipf_val)
ipf_val <- RunUMAP(ipf_val, dims = 1:30, verbose = T, return.model=TRUE)
ipf_val <- FindNeighbors(ipf_val, dims = 1:30, verbose = T)
ipf_val <- FindClusters(ipf_val, verbose = T)

# Covid preprocessing
covid_val <- PercentageFeatureSet(object = covid_val, pattern = "^MT-", col.name = "percent.mt")
covid_val <- subset(covid_val, subset = nFeature_RNA > 1000 & percent.mt < 25)
covid_val <- SCTransform(covid_val, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = T)


saveRDS(ipf_val, file = "ipf_val_preprocess.RDS")
saveRDS(covid_val, file = "covid_val_preprocess.RDS")








