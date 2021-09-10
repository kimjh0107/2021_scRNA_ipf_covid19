library(Seurat)
library(here)
library(tidyverse)
library(cowplot)
library(knitr)
library(sctransform)
library(Matrix)
library(future)
library(future.apply)

# workers
plan('multiprocess', workers = 20)
options(future.globals.maxSize = 2e+05 * 1024^2)

# load RDS 
refquery <- readRDS(here('05_Integrate_df_annotation.RDS'))

# subset to Immune cells 
refquery <- subset(refquery, cell_type_main == 'Immune cells')
refquery <- ScaleData(object = refquery)
refquery <- RunPCA(object = refquery)
refquery <- FindNeighbors(refquery, dims = 1:20)
refquery <- FindClusters(refquery, resolution = 0.6)
refquery <- RunUMAP(object = refquery, dims = 1:20)

# FindConservedImmune_markers
# markers0
DefaultAssay(refquery) <- "RNA"
Immune_markers0 <- FindConservedMarkers(refquery, ident.1 = 0,  grouping.var = "Diagnosis", verbose = T)
head(Immune_markers0)
write.csv(Immune_markers0, 'Immune_markers0.csv', row.names = TRUE)

# markers1
DefaultAssay(refquery) <- "RNA"
Immune_markers1 <- FindConservedMarkers(refquery, ident.1 = 1,  grouping.var = "Diagnosis", verbose = T)
head(Immune_markers1)
write.csv(Immune_markers1, 'Immune_markers1.csv', row.names = TRUE)

# markers2
DefaultAssay(refquery) <- "RNA"
Immune_markers2 <- FindConservedMarkers(refquery, ident.1 = 2,  grouping.var = "Diagnosis", verbose = T)
head(Immune_markers2)
write.csv(Immune_markers2, 'Immune_markers2.csv', row.names = TRUE)

# markers3
DefaultAssay(refquery) <- "RNA"
Immune_markers3 <- FindConservedMarkers(refquery, ident.1 = 3,  grouping.var = "Diagnosis", verbose = T)
head(Immune_markers3)
write.csv(Immune_markers3, 'Immune_markers3.csv', row.names = TRUE)

# markers4
DefaultAssay(refquery) <- "RNA"
Immune_markers4 <- FindConservedMarkers(refquery, ident.1 = 4,  grouping.var = "Diagnosis", verbose = T)
head(Immune_markers4)
write.csv(Immune_markers4, 'Immune_markers4.csv', row.names = TRUE)

# markers5
DefaultAssay(refquery) <- "RNA"
Immune_markers5 <- FindConservedMarkers(refquery, ident.1 = 5,  grouping.var = "Diagnosis", verbose = T)
head(Immune_markers5)
write.csv(Immune_markers5, 'Immune_markers5.csv', row.names = TRUE)

# markers6
DefaultAssay(refquery) <- "RNA"
Immune_markers6 <- FindConservedMarkers(refquery, ident.1 = 6,  grouping.var = "Diagnosis", verbose = T)
head(Immune_markers6)
write.csv(Immune_markers6, 'Immune_markers6.csv', row.names = TRUE)

# markers7
DefaultAssay(refquery) <- "RNA"
Immune_markers7 <- FindConservedMarkers(refquery, ident.1 = 7,  grouping.var = "Diagnosis", verbose = T)
head(Immune_markers7)
write.csv(Immune_markers7, 'Immune_markers7.csv', row.names = TRUE)

# markers8
DefaultAssay(refquery) <- "RNA"
Immune_markers8 <- FindConservedMarkers(refquery, ident.1 = 8,  grouping.var = "Diagnosis", verbose = T)
head(Immune_markers8)
write.csv(Immune_markers8, 'Immune_markers8.csv', row.names = TRUE)

# markers9
DefaultAssay(refquery) <- "RNA"
Immune_markers9 <- FindConservedMarkers(refquery, ident.1 = 9,  grouping.var = "Diagnosis", verbose = T)
head(Immune_markers9)
write.csv(Immune_markers9, 'Immune_markers9.csv', row.names = TRUE)

# markers10
DefaultAssay(refquery) <- "RNA"
Immune_markers10 <- FindConservedMarkers(refquery, ident.1 = 10,  grouping.var = "Diagnosis", verbose = T)
head(Immune_markers10)
write.csv(Immune_markers10, 'Immune_markers10.csv', row.names = TRUE)

# markers11
DefaultAssay(refquery) <- "RNA"
Immune_markers11 <- FindConservedMarkers(refquery, ident.1 = 11,  grouping.var = "Diagnosis", verbose = T)
head(Immune_markers11)
write.csv(Immune_markers11, 'Immune_markers11.csv', row.names = TRUE)

# markers12
DefaultAssay(refquery) <- "RNA"
Immune_markers12 <- FindConservedMarkers(refquery, ident.1 = 12,  grouping.var = "Diagnosis", verbose = T)
head(Immune_markers12)
write.csv(Immune_markers12, 'Immune_markers12.csv', row.names = TRUE)

# markers13
DefaultAssay(refquery) <- "RNA"
Immune_markers13 <- FindConservedMarkers(refquery, ident.1 = 13,  grouping.var = "Diagnosis", verbose = T)
head(Immune_markers13)
write.csv(Immune_markers13, 'Immune_markers13.csv', row.names = TRUE)

# markers14
DefaultAssay(refquery) <- "RNA"
Immune_markers14 <- FindConservedMarkers(refquery, ident.1 = 14,  grouping.var = "Diagnosis", verbose = T)
head(Immune_markers14)
write.csv(Immune_markers14, 'Immune_markers14.csv', row.names = TRUE)

# markers15
DefaultAssay(refquery) <- "RNA"
Immune_markers15 <- FindConservedMarkers(refquery, ident.1 = 15,  grouping.var = "Diagnosis", verbose = T)
head(Immune_markers15)
write.csv(Immune_markers15, 'Immune_markers15.csv', row.names = TRUE)

# markers16
DefaultAssay(refquery) <- "RNA"
Immune_markers16 <- FindConservedMarkers(refquery, ident.1 = 16,  grouping.var = "Diagnosis", verbose = T)
head(Immune_markers16)
write.csv(Immune_markers16, 'Immune_markers16.csv', row.names = TRUE)

# markers17
DefaultAssay(refquery) <- "RNA"
Immune_markers17 <- FindConservedMarkers(refquery, ident.1 = 17,  grouping.var = "Diagnosis", verbose = T)
head(Immune_markers17)
write.csv(Immune_markers17, 'Immune_markers17.csv', row.names = TRUE)

# markers18
DefaultAssay(refquery) <- "RNA"
Immune_markers18 <- FindConservedMarkers(refquery, ident.1 = 18,  grouping.var = "Diagnosis", verbose = T)
head(Immune_markers18)
write.csv(Immune_markers18, 'Immune_markers18.csv', row.names = TRUE)

# markers19
DefaultAssay(refquery) <- "RNA"
Immune_markers19 <- FindConservedMarkers(refquery, ident.1 = 19,  grouping.var = "Diagnosis", verbose = T)
head(Immune_markers19)
write.csv(Immune_markers19, 'Immune_markers19.csv', row.names = TRUE)

# markers20
DefaultAssay(refquery) <- "RNA"
Immune_markers20 <- FindConservedMarkers(refquery, ident.1 = 20,  grouping.var = "Diagnosis", verbose = T)
head(Immune_markers20)
write.csv(Immune_markers20, 'Immune_markers20.csv', row.names = TRUE)

# markers21
DefaultAssay(refquery) <- "RNA"
Immune_markers21 <- FindConservedMarkers(refquery, ident.1 = 21,  grouping.var = "Diagnosis", verbose = T)
head(Immune_markers21)
write.csv(Immune_markers21, 'Immune_markers21.csv', row.names = TRUE)

# markers22
DefaultAssay(refquery) <- "RNA"
Immune_markers21 <- FindConservedMarkers(refquery, ident.1 = 22,  grouping.var = "Diagnosis", verbose = T)
head(Immune_markers22)
write.csv(Immune_markers22, 'Immune_markers22.csv', row.names = TRUE)
