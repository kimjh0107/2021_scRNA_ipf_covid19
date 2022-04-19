library(here)
library(Seurat)
library(tidyverse)
library(ggpubr)


# Sup_Fig1a 
df <- readRDS(here("EMT/RDS/EMT.RDS"))
Idents(df) <- "new_submain"
DefaultAssay(df) <- "RNA"
unique(df$new_submain)
all_markers <- FindAllMarkers(df, logfc.threshold = 0.25, min.pct = 0.25,  test.use = "wilcox")
all_markers$gene <- rownames(all_markers)
write_csv(all_markers, here("EMT/csv/EMT_celltype_findallmarkers.csv"))
