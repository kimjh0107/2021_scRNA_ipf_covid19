library(Seurat)
library(tidyverse)
library(here)
library(ggpubr)

df <- readRDS(here("EMT/RDS/df.RDS"))
### 1
subset <- subset(df, subset = new_submain %in% c("AT2", "Fibroblasts"))
subset <- ScaleData(object = subset)
subset <- RunPCA(object = subset)
subset <- FindNeighbors(subset, dims = 1:30)
subset <- FindClusters(subset, resolution = 0.6)
subset <- RunUMAP(object = subset, dims = 1:30)
saveRDS(subset, file = here("EMT/RDS/at2_fibro.RDS"))


