#!/usr/bin/env Rscript

library(Seurat)
library(here)
library(tidyverse)

# Create Seruat Object
ild_data <- Read10X(here('data/GSE135893_IPF_1/GSE135893_seurat/'), gene.column = 1)
ild <- CreateSeuratObject(ild_data,
                          project = "scRNA_lung", 
                          min.features = 200)

# Add metadata information
metadata <- read_csv(here('data/IPF_metadata.csv'))
ild <- AddMetaData(ild, rownames(ild@meta.data), col.name = "cellbarcodes")
ild <- subset(x = ild, subset = cellbarcodes %in% metadata$cellbarcodes)
cells <- ild@meta.data
cells <- cells %>% left_join(metadata %>% select(Diagnosis, Sample_Name, Sample_Source, Status, celltype, population, cellbarcodes), by = "cellbarcodes")

ild <- AddMetaData(ild, cells$Diagnosis, col.name = "Diagnosis")
ild <- AddMetaData(ild, cells$Sample_Name, col.name = "Sample_Name")
ild <- AddMetaData(ild, cells$Sample_Source, col.name = "Sample_Source")
ild <- AddMetaData(ild, cells$Status, col.name = "Status")
ild <- AddMetaData(ild, cells$celltype, col.name = "celltype")
ild <- AddMetaData(ild, cells$population, col.name = "population")


saveRDS(ild, file = here("data/processed/IPF_SeuratObject.RDS"))