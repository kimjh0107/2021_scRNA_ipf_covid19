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

macrophage <- readRDS(here('macrophage_journal_replicaton.RDS'))
control_covid <- subset(macrophage, subset = Diagnosis == "Control" | Diagnosis == "COVID")

DefaultAssay(macrophage) <- "RNA"
DefaultAssay(control_covid) <- "RNA"

control_covid_ipf_markers <- FindAllMarkers(macrophage, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25, test.use = "MAST")
control_covid_ipf_markers$gene <- rownames(control_covid_ipf_markers)
write_csv(control_covid_ipf_markers, "control_covid_ipf_markers.csv")








