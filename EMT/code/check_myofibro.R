library(here)
library(Seurat)
library(tidyverse)

df <- readRDS(here("EMT/RDS/df.RDS"))
unique(df$new_submain)
### 1
subset <- subset(df, subset = new_submain %in% c("AT2", "Fibroblasts","Smooth Muscle cells"))
subset <- ScaleData(object = subset)
subset <- RunPCA(object = subset)
subset <- FindNeighbors(subset, dims = 1:30)
subset <- FindClusters(subset, resolution = 0.6)
subset <- RunUMAP(object = subset, dims = 1:30)
saveRDS(subset, file = here("EMT/RDS/at2_fibro_smc.RDS"))

"Smooth Muscle cells"