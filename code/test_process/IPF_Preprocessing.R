```{r setup, include=FALSE}
library(Seurat)
library(here)
library(tidyverse)
library(cowplot)
library(knitr)

library(BiocManager)
library(sctransform)
library(glmGamPoi)
library(Matrix)

ild <- readRDS(here('explore/IPF_SeuratObject.RDS'))
head(ild)
```


```{r}
# QC 
ild <- PercentageFeatureSet(ild, pattern = "^MT-", col.name = "percent.mt")
ild <- subset(ild, subset = nFeature_RNA > 1000 & percent.mt < 25)

# SCTransform & Preprocessing
ild <- SCTransform(ild, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)

ild <- RunPCA(ild)
ild <- RunUMAP(ild, dims = 1:30, verbose = FALSE, return.model=TRUE)
ild <- FindNeighbors(ild, dims = 1:30, verbose = FALSE)
ild <- FindClusters(ild, verbose = FALSE)
```

```{r}
saveRDS(ild, file = "IPF_Preprocessing.RDS")
```