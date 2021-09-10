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



df <- readRDS(here('GSE135893_ILD_annotated_fullsize.rds'))
refquery <- readRDS(here('05_Integrate_df_annotation.RDS'))

df<- AddMetaData(df, rownames(df@meta.data), col.name = "cellbarcodes")


anchors <- FindTransferAnchors(
  reference = df,
  query = refquery,
  normalization.method = "LogNormalize",
  reference.reduction = "pca",
  dims = 1:50, verbose = TRUE)
  
refquery <- MapQuery(
    anchorset = anchors,
    query = refquery,
    reference = df,
    refdata = list(
      celltype = "celltype",
    ),
    reference.reduction = "pca", 
    reduction.model = "umap", verbose = TRUE)


df$id <- 'reference'
refquery$id <- 'query'
merge <- merge(df, refquery)

merge <- ScaleData(object = merge)
merge <- RunPCA(object = merge)
merge <- FindNeighbors(merge, dims = 1:20)
merge <- FindClusters(merge, resolution = 0.6)
merge <- RunUMAP(object = merge, dims = 1:20)

saveRDS(merge, file = "Merge.RDS")
