library(Seurat)
library(tidyverse)
library(here)
library(ggpubr)

df <- readRDS(here("EMT/RDS/df.RDS"))
### 1
subset <- subset(df, subset = new_submain %in% c("AT1", "AT2", "Fibroblasts"))
subset <- ScaleData(object = subset)
subset <- RunPCA(object = subset)
subset <- FindNeighbors(subset, dims = 1:30)
subset <- FindClusters(subset, resolution = 0.6)
subset <- RunUMAP(object = subset, dims = 1:30)
saveRDS(subset, file = here("EMT/RDS/EMT.RDS"))
DimPlot(subset, group.by = "new_submain", label = T, repel = T) + NoLegend()
ggsave(here("EMT/figure/DimPlot_emt_celltype.pdf"), width = 7, height = 5)



### 1
subset <- subset(df, subset = new_submain %in% c("AT1", "AT2"))
subset <- ScaleData(object = subset)
subset <- RunPCA(object = subset)
subset <- FindNeighbors(subset, dims = 1:30)
subset <- FindClusters(subset, resolution = 0.6)
subset <- RunUMAP(object = subset, dims = 1:30)
saveRDS(subset, file = here("EMT/RDS/AT.RDS"))
DimPlot(subset, group.by = "new_submain", label = T, repel = T) + NoLegend()
ggsave(here("EMT/figure/DimPlot_emt_celltype.pdf"), width = 7, height = 5)



### 1
subset <- subset(df, subset = new_submain %in% c("Fibroblasts"))
subset <- ScaleData(object = subset)
subset <- RunPCA(object = subset)
subset <- FindNeighbors(subset, dims = 1:30)
subset <- FindClusters(subset, resolution = 0.6)
subset <- RunUMAP(object = subset, dims = 1:30)
saveRDS(subset, file = here("EMT/RDS/fibro.RDS"))
DimPlot(subset, group.by = "new_submain", label = T, repel = T) + NoLegend()
ggsave(here("EMT/figure/DimPlot_emt_celltype.pdf"), width = 7, height = 5)