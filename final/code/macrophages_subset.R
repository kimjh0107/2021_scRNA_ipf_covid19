library(here)
library(Seurat)
library(tidyverse)

df <- readRDS(here("final/RDS/scRNA.RDS"))

#### Macrophage & Monocytes #### 
subset <- subset(df, subset = new_submain %in% c("Monocytes", "Macrophages"))

subset <- ScaleData(object = subset)
subset <- RunPCA(object = subset)
subset <- FindNeighbors(subset, dims = 1:30)
subset <- FindClusters(subset, resolution = 0.6)
subset <- RunUMAP(object = subset, dims = 1:30)

saveRDS(subset, file = here("final/RDS/macro_mono_subset.RDS"))

DimPlot(subset, group.by = "seurat_clusters", label = T, repel = T)
ggsave(here("final/figure/DimPlot_macro_subset.pdf"), width = 7, height = 5)
