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

# subset to Mesenchymal cells 
refquery <- subset(refquery, cell_type_main == 'Mesenchymal cells')
refquery <- ScaleData(object = refquery)
refquery <- RunPCA(object = refquery)
refquery <- FindNeighbors(refquery, dims = 1:20)
refquery <- FindClusters(refquery, resolution = 0.6)
refquery <- RunUMAP(object = refquery, dims = 1:20)

# FindConservedMesenchymal_markers
# markers0
DefaultAssay(refquery) <- "RNA"
Mesenchymal_markers0 <- FindConservedMarkers(refquery, ident.1 = 0,  grouping.var = "Diagnosis", verbose = T)
head(Mesenchymal_markers0)
write.csv(Mesenchymal_markers0, 'Mesenchymal_markers0.csv', row.names = TRUE)

# markers1
DefaultAssay(refquery) <- "RNA"
Mesenchymal_markers1 <- FindConservedMarkers(refquery, ident.1 = 1,  grouping.var = "Diagnosis", verbose = T)
head(Mesenchymal_markers1)
write.csv(Mesenchymal_markers1, 'Mesenchymal_markers1.csv', row.names = TRUE)

# markers2
DefaultAssay(refquery) <- "RNA"
Mesenchymal_markers2 <- FindConservedMarkers(refquery, ident.1 = 2,  grouping.var = "Diagnosis", verbose = T)
head(Mesenchymal_markers2)
write.csv(Mesenchymal_markers2, 'Mesenchymal_markers2.csv', row.names = TRUE)

# markers3
DefaultAssay(refquery) <- "RNA"
Mesenchymal_markers3 <- FindConservedMarkers(refquery, ident.1 = 3,  grouping.var = "Diagnosis", verbose = T)
head(Mesenchymal_markers3)
write.csv(Mesenchymal_markers3, 'Mesenchymal_markers3.csv', row.names = TRUE)

# markers4
DefaultAssay(refquery) <- "RNA"
Mesenchymal_markers4 <- FindConservedMarkers(refquery, ident.1 = 4,  grouping.var = "Diagnosis", verbose = T)
head(Mesenchymal_markers4)
write.csv(Mesenchymal_markers4, 'Mesenchymal_markers4.csv', row.names = TRUE)

# markers5
DefaultAssay(refquery) <- "RNA"
Mesenchymal_markers5 <- FindConservedMarkers(refquery, ident.1 = 5,  grouping.var = "Diagnosis", verbose = T)
head(Mesenchymal_markers5)
write.csv(Mesenchymal_markers5, 'Mesenchymal_markers5.csv', row.names = TRUE)

# markers6
DefaultAssay(refquery) <- "RNA"
Mesenchymal_markers6 <- FindConservedMarkers(refquery, ident.1 = 6,  grouping.var = "Diagnosis", verbose = T)
head(Mesenchymal_markers6)
write.csv(Mesenchymal_markers6, 'Mesenchymal_markers6.csv', row.names = TRUE)

# markers7
DefaultAssay(refquery) <- "RNA"
Mesenchymal_markers7 <- FindConservedMarkers(refquery, ident.1 = 7,  grouping.var = "Diagnosis", verbose = T)
head(Mesenchymal_markers7)
write.csv(Mesenchymal_markers7, 'Mesenchymal_markers7.csv', row.names = TRUE)

# markers8
DefaultAssay(refquery) <- "RNA"
Mesenchymal_markers8 <- FindConservedMarkers(refquery, ident.1 = 8,  grouping.var = "Diagnosis", verbose = T)
head(Mesenchymal_markers8)
write.csv(Mesenchymal_markers8, 'Mesenchymal_markers8.csv', row.names = TRUE)

# markers9
DefaultAssay(refquery) <- "RNA"
Mesenchymal_markers9 <- FindConservedMarkers(refquery, ident.1 = 9,  grouping.var = "Diagnosis", verbose = T)
head(Mesenchymal_markers9)
write.csv(Mesenchymal_markers9, 'Mesenchymal_markers9.csv', row.names = TRUE)

# markers10
DefaultAssay(refquery) <- "RNA"
Mesenchymal_markers10 <- FindConservedMarkers(refquery, ident.1 = 10,  grouping.var = "Diagnosis", verbose = T)
head(Mesenchymal_markers10)
write.csv(Mesenchymal_markers10, 'Mesenchymal_markers10.csv', row.names = TRUE)

# markers11
DefaultAssay(refquery) <- "RNA"
Mesenchymal_markers11 <- FindConservedMarkers(refquery, ident.1 = 11,  grouping.var = "Diagnosis", verbose = T)
head(Mesenchymal_markers11)
write.csv(Mesenchymal_markers11, 'Mesenchymal_markers11.csv', row.names = TRUE)

# markers12
DefaultAssay(refquery) <- "RNA"
Mesenchymal_markers12 <- FindConservedMarkers(refquery, ident.1 = 12,  grouping.var = "Diagnosis", verbose = T)
head(Mesenchymal_markers12)
write.csv(Mesenchymal_markers12, 'Mesenchymal_markers12.csv', row.names = TRUE)

# markers13
DefaultAssay(refquery) <- "RNA"
Mesenchymal_markers13 <- FindConservedMarkers(refquery, ident.1 = 13,  grouping.var = "Diagnosis", verbose = T)
head(Mesenchymal_markers13)
write.csv(Mesenchymal_markers13, 'Mesenchymal_markers13.csv', row.names = TRUE)
