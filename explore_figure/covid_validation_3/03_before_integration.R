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


covid <- readRDS(here('02_COVID_make_Diagnosis.RDS'))
ipf <- readRDS(here('01_IPF_SeuratObject.RDS'))

# IPF preprocessing 
ipf <- PercentageFeatureSet(ipf, pattern = "^MT-", col.name = "percent.mt")
ipf <- subset(ipf, subset = nFeature_RNA > 1000 & percent.mt < 25)


# Covid preprocessing
covid <- PercentageFeatureSet(object = covid, pattern = "^MT-", col.name = "percent.mt")
covid <- subset(covid, subset = nFeature_RNA > 1000 & percent.mt < 25)


saveRDS(ipf, file = "ipf_preprocess.RDS")
saveRDS(covid, file = "covid_preprocess.RDS")

