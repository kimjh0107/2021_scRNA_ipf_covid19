library(tidyverse)
library(Seurat)
library(here)
emt <- readRDS(here("EMT/RDS/emt_new_annotation.RDS"))
emt_sn <- readRDS(here("EMT/RDS/emt_snRNA_new_annotation.RDS"))


Idents(emt) <- "seurat_clusters"
DefaultAssay(emt) <- 'RNA'
all_markers <- FindAllMarkers(emt, logfc.threshold = 0.25, min.pct = 0.25)
all_markers$gene <- rownames(all_markers)
write_csv(all_markers, here("EMT/csv/EMT_seurat_clusters_celltype_findallmarkers.csv"))


Idents(emt_sn) <- "seurat_clusters"
DefaultAssay(emt_sn) <- 'RNA'
all_markers <- FindAllMarkers(emt_sn, logfc.threshold = 0.25, min.pct = 0.25)
all_markers$gene <- rownames(all_markers)
write_csv(all_markers, here("EMT/csv/EMT_snRNA_seurat_clusters_celltype_findallmarkers.csv"))