---
  title: "Mapping_process"
output: html_document
output: github_document
---
  
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
```

```{r}
ild<- readRDS(here('explore/RDS/IPF_Preprocessing.RDS'))
covid <- readRDS(here('explore/RDS/COVID_Preprocessing.RDS'))


anchors <- FindTransferAnchors(
  reference = ild,
  query = covid,
  normalization.method = "SCT",
  reference.reduction = "pca",
  dims = 1:50
)


covid <- MapQuery(
  anchorset = anchors,
  query = covid,
  reference = ild,
  refdata = list(
    celltype = "celltype",
    predicted_SCT = "SCT"
  ),
  reference.reduction = "pca", 
  reduction.model = "umap"
)


saveRDS(covid, file = "COVID_Mapping.RDS")
```


