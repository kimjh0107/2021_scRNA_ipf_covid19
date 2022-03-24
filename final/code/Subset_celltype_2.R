library(Seurat)
library(tidyverse)
library(here)
library(ggpubr)

df <- readRDS(here("final/RDS/df.RDS"))

### 1
subset <- subset(df, subset = new_submain %in% c("Mast cells"))
subset <- ScaleData(object = subset)
subset <- RunPCA(object = subset)
subset <- FindNeighbors(subset, dims = 1:30)
subset <- FindClusters(subset, resolution = 0.6)
subset <- RunUMAP(object = subset, dims = 1:30)
saveRDS(subset, file = here("final/RDS/mast.RDS"))
DimPlot(subset, group.by = "seurat_clusters", label = T, repel = T) + NoLegend()
ggsave(here("final/figure/DimPlot_mast_celltype.pdf"), width = 7, height = 5)


### 1
subset <- subset(df, subset = new_submain %in% c("B cells"))
subset <- ScaleData(object = subset)
subset <- RunPCA(object = subset)
subset <- FindNeighbors(subset, dims = 1:30)
subset <- FindClusters(subset, resolution = 0.6)
subset <- RunUMAP(object = subset, dims = 1:30)
saveRDS(subset, file = here("final/RDS/bcell.RDS"))
DimPlot(subset, group.by = "seurat_clusters", label = T, repel = T) + NoLegend()
ggsave(here("final/figure/DimPlot_mast_celltype.pdf"), width = 7, height = 5)


### 1
subset <- subset(df, subset = new_submain %in% c("Neutrophils"))
subset <- ScaleData(object = subset)
subset <- RunPCA(object = subset)
subset <- FindNeighbors(subset, dims = 1:30)
subset <- FindClusters(subset, resolution = 0.6)
subset <- RunUMAP(object = subset, dims = 1:30)
saveRDS(subset, file = here("final/RDS/neutrophil.RDS"))
DimPlot(subset, group.by = "seurat_clusters", label = T, repel = T) + NoLegend()
ggsave(here("final/figure/DimPlot_neutrophil_celltype.pdf"), width = 7, height = 5)