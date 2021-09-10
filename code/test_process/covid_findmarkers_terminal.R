library(Seurat)
library(here)
library(tidyverse)
library(cowplot)
library(knitr)

library(BiocManager)
library(sctransform)
library(glmGamPoi)
library(Matrix)


covid <- readRDS(here('explore/RDS/COVID_Clustering.RDS'))

cluster0.markers <- FindMarkers(covid, ident.1 = 0, min.pct = 0.25)
head(clister0.markers, n =5)

cluster1.markers <- FindMarkers(covid, ident.1 = 1, min.pct = 0.25)
head(clister1.markers, n =5)

cluster2.markers <- FindMarkers(covid, ident.1 = 2, min.pct = 0.25)
head(clister2.markers, n =5)

cluster3.markers <- FindMarkers(covid, ident.1 = 3, min.pct = 0.25)
head(clister3.markers, n =5)

cluster4.markers <- FindMarkers(covid, ident.1 = 4, min.pct = 0.25)
head(clister4.markers, n =5)

cluster5.markers <- FindMarkers(covid, ident.1 = 5, min.pct = 0.25)
head(clister5.markers, n =5)

cluster6.markers <- FindMarkers(covid, ident.1 = 6, min.pct = 0.25)
head(clister6.markers, n =5)

cluster7.markers <- FindMarkers(covid, ident.1 = 7, min.pct = 0.25)
head(clister7.markers, n =5)

cluster8.markers <- FindMarkers(covid, ident.1 = 8, min.pct = 0.25)
head(clister8.markers, n =5)

cluster9.markers <- FindMarkers(covid, ident.1 = 9, min.pct = 0.25)
head(clister9.markers, n =5)

cluster10.markers <- FindMarkers(covid, ident.1 = 10, min.pct = 0.25)
head(clister10.markers, n =5)

cluster11.markers <- FindMarkers(covid, ident.11 = 11, min.pct = 0.25)
head(clister11.markers, n =5)

cluster12.markers <- FindMarkers(covid, ident.12 = 12, min.pct = 0.25)
head(clister12.markers, n =5)

cluster13.markers <- FindMarkers(covid, ident.13 = 13, min.pct = 0.25)
head(clister13.markers, n =5)

cluster14.markers <- FindMarkers(covid, ident.14 = 14, min.pct = 0.25)
head(clister14.markers, n =5)

cluster15.markers <- FindMarkers(covid, ident.15 = 15, min.pct = 0.25)
head(clister15.markers, n =5)

cluster16.markers <- FindMarkers(covid, ident.16 = 16, min.pct = 0.25)
head(clister16.markers, n =5)

cluster17.markers <- FindMarkers(covid, ident.17 = 17, min.pct = 0.25)
head(clister17.markers, n =5)

cluster18.markers <- FindMarkers(covid, ident.18 = 18, min.pct = 0.25)
head(clister18.markers, n =5)

cluster19.markers <- FindMarkers(covid, ident.19 = 19, min.pct = 0.25)
head(clister19.markers, n =5)

cluster20.markers <- FindMarkers(covid, ident.20 = 20, min.pct = 0.25)
head(clister20.markers, n =5)

cluster21.markers <- FindMarkers(covid, ident.21 = 21, min.pct = 0.25)
head(clister21.markers, n =5)

cluster22.markers <- FindMarkers(covid, ident.22 = 22, min.pct = 0.25)
head(clister22.markers, n =5)

cluster23.markers <- FindMarkers(covid, ident.23 = 23, min.pct = 0.25)
head(clister23.markers, n =5)

cluster24.markers <- FindMarkers(covid, ident.24 = 24, min.pct = 0.25)
head(clister24.markers, n =5)

cluster25.markers <- FindMarkers(covid, ident.25 = 25, min.pct = 0.25)
head(clister25.markers, n =5)

