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
covid <- readRDS(here('explore/RDS/COVID_Mapping.RDS'))


#merge reference and query
ild$id <- 'reference'
covid$id <- 'query'
refquery <- merge(ild, covid)
refquery[["pca"]] <- merge(ild[["pca"]], covid[["ref.pca"]])
refquery <- RunUMAP(refquery, reduction = 'pca', dims = 1:50)
DimPlot(refquery, group.by = 'id', shuffle = TRUE)
```
