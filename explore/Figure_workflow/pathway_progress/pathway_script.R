library(Seurat)
library(here)
library(tidyverse)
library(cowplot)
library(knitr)

library(BiocManager)
library(sctransform)
library(glmGamPoi)
library(Matrix)
library(future)
library(future.apply)

plan('multiprocess', workers = 20)
options(future.globals.maxSize = 2e+05 * 1024^2)

pathway <- readRDS(here('pathway.RDS'))

pathway_markers <- FindAllMarkers(pathway, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(pathway_markers, "pathway_markers.csv")

