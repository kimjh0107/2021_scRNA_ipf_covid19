library(Seurat)
library(tidyverse)
library(here)
library(SingleCellExperiment)
library(cowplot)
library(knitr)
library(sctransform)
library(Matrix)
library(future)
library(future.apply)

plan('multiprocess', workers = 40)
options(future.globals.maxSize = 2e+05 * 1024^2)


macrophage <- readRDS(here('macrophage_journal_replicaton3.RDS'))

Idents(macrophage) <- "new_submain"
DefaultAssay(macrophage) <- "integrated"

replication_markers <- FindAllMarkers(macrophage, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25, test.use = "MAST")
replication_markers$gene <- rownames(replication_markers)
write_csv(replication_markers, "replication_markers_mast.csv")








