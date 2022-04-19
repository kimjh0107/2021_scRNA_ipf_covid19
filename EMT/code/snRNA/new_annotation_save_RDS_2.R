library(tidyverse)
library(Seurat)
library(here)

emt <- readRDS(here("EMT/RDS/EMT.RDS"))

emt@meta.data$new_submain2 <- ifelse(emt@meta.data$seurat_clusters %in% c(10,8,17,6), "AT1", "NA")
emt@meta.data$new_submain2 <- ifelse(emt@meta.data$seurat_clusters %in% c(1,0,11,4,5,2,9,13,16), "AT2", emt@meta.data$new_submain2)
emt@meta.data$new_submain2 <- ifelse(emt@meta.data$seurat_clusters %in% c(), "Myofibroblasts", emt@meta.data$new_submain2)
emt@meta.data$new_submain2 <- ifelse(emt@meta.data$seurat_clusters %in% c(), "Fibroblasts", emt@meta.data$new_submain2)
saveRDS(emt, here("EMT/RDS/emt_new_annotation.RDS"))


emt <- readRDS(here("EMT/RDS/snRNA_EMT.RDS"))

emt@meta.data$new_submain3 <- ifelse(emt@meta.data$seurat_clusters %in% c(20,19,14,18,3), "AT1", "NA")
emt@meta.data$new_submain3 <- ifelse(emt@meta.data$seurat_clusters %in% c(16,12,0,15,11,7,2), "AT2", emt@meta.data$new_submain3)
emt@meta.data$new_submain3 <- ifelse(emt@meta.data$seurat_clusters %in% c(), "Myofibroblasts", emt@meta.data$new_submain3)
emt@meta.data$new_submain3 <- ifelse(emt@meta.data$seurat_clusters %in% c(), "Fibroblasts", emt@meta.data$new_submain3)
saveRDS(emt, here("EMT/RDS/emt_snRNA_new_annotation.RDS"))