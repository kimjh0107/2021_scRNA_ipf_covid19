library(here)
library(Seurat)
library(tidyverse)
library(ggpubr)


# Sup_Fig1a 
df <- readRDS(here("EMT/RDS/at2_fibro_basal.RDS"))

DefaultAssay(df) <- "RNA"

all_markers <- FindAllMarkers(df, min.pct = 0.1, test.use = "wilcox")
all_markers$gene <- rownames(all_markers)
write_csv(all_markers, here("EMT/csv/findallmarkers_at2_basal_fibro.csv"))

