library(here)
library(Seurat)
library(tidyverse)

emt <- readRDS(here("EMT/RDS/emt_new_annotation.RDS"))
subset <- subset(emt, subset = new_submain2 %in% c("Fibroblasts", "Myofibroblasts"))
subset <- ScaleData(object = subset)
subset <- RunPCA(object = subset)
subset <- FindNeighbors(subset, dims = 1:30)
subset <- FindClusters(subset, resolution = 0.6)
subset <- RunUMAP(object = subset, dims = 1:30)
saveRDS(subset, file = here("EMT/RDS/fibro_annotation.RDS"))


emt <- readRDS(here("EMT/RDS/emt_snRNA_new_annotation.RDS"))
subset <- subset(emt, subset = new_submain3 %in% c("Fibroblasts", "Myofibroblasts"))
subset <- ScaleData(object = subset)
subset <- RunPCA(object = subset)
subset <- FindNeighbors(subset, dims = 1:30)
subset <- FindClusters(subset, resolution = 0.6)
subset <- RunUMAP(object = subset, dims = 1:30)
saveRDS(subset, file = here("EMT/RDS/snRNA_fibro_annotation.RDS"))

