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

# subset to Endothelial cells 
refquery <- subset(refquery, cell_type_main == 'Endothelial cells')
refquery <- ScaleData(object = refquery)
refquery <- RunPCA(object = refquery)
refquery <- FindNeighbors(refquery, dims = 1:20)
refquery <- FindClusters(refquery, resolution = 0.6)
refquery <- RunUMAP(object = refquery, dims = 1:20)

# FindConservedEndothelial_markers
# markers0
DefaultAssay(refquery) <- "RNA"
Endothelial_markers0 <- FindConservedMarkers(refquery, ident.1 = 0,  grouping.var = "Diagnosis", verbose = T)
head(Endothelial_markers0)
write.csv(Endothelial_markers0, 'Endothelial_markers0.csv', row.names = TRUE)

# markers1
DefaultAssay(refquery) <- "RNA"
Endothelial_markers1 <- FindConservedMarkers(refquery, ident.1 = 1,  grouping.var = "Diagnosis", verbose = T)
head(Endothelial_markers1)
write.csv(Endothelial_markers1, 'Endothelial_markers1.csv', row.names = TRUE)

# markers2
DefaultAssay(refquery) <- "RNA"
Endothelial_markers2 <- FindConservedMarkers(refquery, ident.1 = 2,  grouping.var = "Diagnosis", verbose = T)
head(Endothelial_markers2)
write.csv(Endothelial_markers2, 'Endothelial_markers2.csv', row.names = TRUE)

# markers3
DefaultAssay(refquery) <- "RNA"
Endothelial_markers3 <- FindConservedMarkers(refquery, ident.1 = 3,  grouping.var = "Diagnosis", verbose = T)
head(Endothelial_markers3)
write.csv(Endothelial_markers3, 'Endothelial_markers3.csv', row.names = TRUE)

# markers4
DefaultAssay(refquery) <- "RNA"
Endothelial_markers4 <- FindConservedMarkers(refquery, ident.1 = 4,  grouping.var = "Diagnosis", verbose = T)
head(Endothelial_markers4)
write.csv(Endothelial_markers4, 'Endothelial_markers4.csv', row.names = TRUE)

# markers5
DefaultAssay(refquery) <- "RNA"
Endothelial_markers5 <- FindConservedMarkers(refquery, ident.1 = 5,  grouping.var = "Diagnosis", verbose = T)
head(Endothelial_markers5)
write.csv(Endothelial_markers5, 'Endothelial_markers5.csv', row.names = TRUE)

# markers6
DefaultAssay(refquery) <- "RNA"
Endothelial_markers6 <- FindConservedMarkers(refquery, ident.1 = 6,  grouping.var = "Diagnosis", verbose = T)
head(Endothelial_markers6)
write.csv(Endothelial_markers6, 'Endothelial_markers6.csv', row.names = TRUE)

# markers7
DefaultAssay(refquery) <- "RNA"
Endothelial_markers7 <- FindConservedMarkers(refquery, ident.1 = 7,  grouping.var = "Diagnosis", verbose = T)
head(Endothelial_markers7)
write.csv(Endothelial_markers7, 'Endothelial_markers7.csv', row.names = TRUE)

# markers8
DefaultAssay(refquery) <- "RNA"
Endothelial_markers8 <- FindConservedMarkers(refquery, ident.1 = 8,  grouping.var = "Diagnosis", verbose = T)
head(Endothelial_markers8)
write.csv(Endothelial_markers8, 'Endothelial_markers8.csv', row.names = TRUE)

# markers9
DefaultAssay(refquery) <- "RNA"
Endothelial_markers9 <- FindConservedMarkers(refquery, ident.1 = 9,  grouping.var = "Diagnosis", verbose = T)
head(Endothelial_markers9)
write.csv(Endothelial_markers9, 'Endothelial_markers9.csv', row.names = TRUE)

# markers10
DefaultAssay(refquery) <- "RNA"
Endothelial_markers10 <- FindConservedMarkers(refquery, ident.1 = 10,  grouping.var = "Diagnosis", verbose = T)
head(Endothelial_markers10)
write.csv(Endothelial_markers10, 'Endothelial_markers10.csv', row.names = TRUE)

# markers11
DefaultAssay(refquery) <- "RNA"
Endothelial_markers11 <- FindConservedMarkers(refquery, ident.1 = 11,  grouping.var = "Diagnosis", verbose = T)
head(Endothelial_markers11)
write.csv(Endothelial_markers11, 'Endothelial_markers11.csv', row.names = TRUE)

# markers12
DefaultAssay(refquery) <- "RNA"
Endothelial_markers12 <- FindConservedMarkers(refquery, ident.1 = 12,  grouping.var = "Diagnosis", verbose = T)
head(Endothelial_markers12)
write.csv(Endothelial_markers12, 'Endothelial_markers12.csv', row.names = TRUE)

# markers13
DefaultAssay(refquery) <- "RNA"
Endothelial_markers13 <- FindConservedMarkers(refquery, ident.1 = 13,  grouping.var = "Diagnosis", verbose = T)
head(Endothelial_markers13)
write.csv(Endothelial_markers13, 'Endothelial_markers13.csv', row.names = TRUE)

# markers14
DefaultAssay(refquery) <- "RNA"
Endothelial_markers14 <- FindConservedMarkers(refquery, ident.1 = 14,  grouping.var = "Diagnosis", verbose = T)
head(Endothelial_markers14)
write.csv(Endothelial_markers14, 'Endothelial_markers14.csv', row.names = TRUE)
