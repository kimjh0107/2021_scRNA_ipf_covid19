library(Seurat)
library(tidyverse)
library(here)
library(ggpubr)

df <- readRDS(here("final/RDS/df.RDS"))

### 1
subset <- subset(df, subset = new_submain %in% c("Macrophages", "Monocytes"))
subset <- ScaleData(object = subset)
subset <- RunPCA(object = subset)
subset <- FindNeighbors(subset, dims = 1:30)
subset <- FindClusters(subset, resolution = 0.6)
subset <- RunUMAP(object = subset, dims = 1:30)
saveRDS(subset, file = here("final/RDS/macro_mono.RDS"))
DimPlot(subset, group.by = "new_submain", label = T, repel = T) + NoLegend()
ggsave(here("final/figure/DimPlot_macro_mono_celltype.pdf"), width = 7, height = 5)

#### 2
subset <- subset(df, subset = new_submain %in% c("Macrophages"))
subset <- ScaleData(object = subset)
subset <- RunPCA(object = subset)
subset <- FindNeighbors(subset, dims = 1:30)
subset <- FindClusters(subset, resolution = 0.6)
subset <- RunUMAP(object = subset, dims = 1:30)
saveRDS(subset, file = here("final/RDS/macro.RDS"))
DimPlot(subset, group.by = "new_submain", label = T, repel = T) + NoLegend()
ggsave(here("final/figure/DimPlot_macro_celltype.pdf"), width = 7, height = 5)

#### 3 
subset <- subset(df, subset = new_submain %in% c("Macrophages", "Monocytes", "mDCs"))
subset <- ScaleData(object = subset)
subset <- RunPCA(object = subset)
subset <- FindNeighbors(subset, dims = 1:30)
subset <- FindClusters(subset, resolution = 0.6)
subset <- RunUMAP(object = subset, dims = 1:30)
saveRDS(subset, file = here("final/RDS/macro_mono_mdcs.RDS"))
DimPlot(subset, group.by = "new_submain", label = T, repel = T) + NoLegend()
ggsave(here("final/figure/DimPlot_macro_mono_mdc_celltype.pdf"), width = 7, height = 5)

#### 3 
subset <- subset(df, subset = new_submain %in% c("CD4 Tcells", "CD8 Tcells", "NK cells"))
subset <- ScaleData(object = subset)
subset <- RunPCA(object = subset)
subset <- FindNeighbors(subset, dims = 1:30)
subset <- FindClusters(subset, resolution = 0.6)
subset <- RunUMAP(object = subset, dims = 1:30)
saveRDS(subset, file = here("final/RDS/T_NK.RDS"))
DimPlot(subset, group.by = "new_submain", label = T, repel = T) + NoLegend()
ggsave(here("final/figure/DimPlot_T_nk_celltype.pdf"), width = 7, height = 5)

#### 3 
subset <- subset(df, subset = new_submain %in% c("AT1", "AT2"))
subset <- ScaleData(object = subset)
subset <- RunPCA(object = subset)
subset <- FindNeighbors(subset, dims = 1:30)
subset <- FindClusters(subset, resolution = 0.6)
subset <- RunUMAP(object = subset, dims = 1:30)
saveRDS(subset, file = here("final/RDS/AT.RDS"))
DimPlot(subset, group.by = "new_submain", label = T, repel = T) + NoLegend()
ggsave(here("final/figure/DimPlot_at_celltype.pdf"), width = 7, height = 5)


#### 3 
subset <- subset(df, subset = new_submain %in% c("Fibroblasts"))
subset <- ScaleData(object = subset)
subset <- RunPCA(object = subset)
subset <- FindNeighbors(subset, dims = 1:30)
subset <- FindClusters(subset, resolution = 0.6)
subset <- RunUMAP(object = subset, dims = 1:30)
saveRDS(subset, file = here("final/RDS/fibro.RDS"))
DimPlot(subset, group.by = "new_submain", label = T, repel = T) + NoLegend()
ggsave(here("final/figure/DimPlot_fibro_celltype.pdf"), width = 7, height = 5)