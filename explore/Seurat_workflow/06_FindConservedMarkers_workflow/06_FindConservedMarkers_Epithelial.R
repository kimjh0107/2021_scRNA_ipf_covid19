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

# subset to Epithelial cells 
refquery <- subset(refquery, cell_type_main == 'Epithelial cells')
refquery <- ScaleData(object = refquery)
refquery <- RunPCA(object = refquery)
refquery <- FindNeighbors(refquery, dims = 1:20)
refquery <- FindClusters(refquery, resolution = 0.6)
refquery <- RunUMAP(object = refquery, dims = 1:20)

# FindConservedepithelial_markers
# markers0
DefaultAssay(refquery) <- "RNA"
epithelial_markers0 <- FindConservedMarkers(refquery, ident.1 = 0,  grouping.var = "Diagnosis", verbose = T)
head(epithelial_markers0)
write.csv(epithelial_markers0, 'epithelial_markers0.csv', row.names = TRUE)

# markers1
DefaultAssay(refquery) <- "RNA"
epithelial_markers1 <- FindConservedMarkers(refquery, ident.1 = 1,  grouping.var = "Diagnosis", verbose = T)
head(epithelial_markers1)
write.csv(epithelial_markers1, 'epithelial_markers1.csv', row.names = TRUE)

# markers2
DefaultAssay(refquery) <- "RNA"
epithelial_markers2 <- FindConservedMarkers(refquery, ident.1 = 2,  grouping.var = "Diagnosis", verbose = T)
head(epithelial_markers2)
write.csv(epithelial_markers2, 'epithelial_markers2.csv', row.names = TRUE)

# markers3
DefaultAssay(refquery) <- "RNA"
epithelial_markers3 <- FindConservedMarkers(refquery, ident.1 = 3,  grouping.var = "Diagnosis", verbose = T)
head(epithelial_markers3)
write.csv(epithelial_markers3, 'epithelial_markers3.csv', row.names = TRUE)

# markers4
DefaultAssay(refquery) <- "RNA"
epithelial_markers4 <- FindConservedMarkers(refquery, ident.1 = 4,  grouping.var = "Diagnosis", verbose = T)
head(epithelial_markers4)
write.csv(epithelial_markers4, 'epithelial_markers4.csv', row.names = TRUE)

# markers5
DefaultAssay(refquery) <- "RNA"
epithelial_markers5 <- FindConservedMarkers(refquery, ident.1 = 5,  grouping.var = "Diagnosis", verbose = T)
head(epithelial_markers5)
write.csv(epithelial_markers5, 'epithelial_markers5.csv', row.names = TRUE)

# markers6
DefaultAssay(refquery) <- "RNA"
epithelial_markers6 <- FindConservedMarkers(refquery, ident.1 = 6,  grouping.var = "Diagnosis", verbose = T)
head(epithelial_markers6)
write.csv(epithelial_markers6, 'epithelial_markers6.csv', row.names = TRUE)

# markers7
DefaultAssay(refquery) <- "RNA"
epithelial_markers7 <- FindConservedMarkers(refquery, ident.1 = 7,  grouping.var = "Diagnosis", verbose = T)
head(epithelial_markers7)
write.csv(epithelial_markers7, 'epithelial_markers7.csv', row.names = TRUE)

# markers8
DefaultAssay(refquery) <- "RNA"
epithelial_markers8 <- FindConservedMarkers(refquery, ident.1 = 8,  grouping.var = "Diagnosis", verbose = T)
head(epithelial_markers8)
write.csv(epithelial_markers8, 'epithelial_markers8.csv', row.names = TRUE)

# markers9
DefaultAssay(refquery) <- "RNA"
epithelial_markers9 <- FindConservedMarkers(refquery, ident.1 = 9,  grouping.var = "Diagnosis", verbose = T)
head(epithelial_markers9)
write.csv(epithelial_markers9, 'epithelial_markers9.csv', row.names = TRUE)

# markers10
DefaultAssay(refquery) <- "RNA"
epithelial_markers10 <- FindConservedMarkers(refquery, ident.1 = 10,  grouping.var = "Diagnosis", verbose = T)
head(epithelial_markers10)
write.csv(epithelial_markers10, 'epithelial_markers10.csv', row.names = TRUE)

# markers11
DefaultAssay(refquery) <- "RNA"
epithelial_markers11 <- FindConservedMarkers(refquery, ident.1 = 11,  grouping.var = "Diagnosis", verbose = T)
head(epithelial_markers11)
write.csv(epithelial_markers11, 'epithelial_markers11.csv', row.names = TRUE)

# markers12
DefaultAssay(refquery) <- "RNA"
epithelial_markers12 <- FindConservedMarkers(refquery, ident.1 = 12,  grouping.var = "Diagnosis", verbose = T)
head(epithelial_markers12)
write.csv(epithelial_markers12, 'epithelial_markers12.csv', row.names = TRUE)

# markers13
DefaultAssay(refquery) <- "RNA"
epithelial_markers13 <- FindConservedMarkers(refquery, ident.1 = 13,  grouping.var = "Diagnosis", verbose = T)
head(epithelial_markers13)
write.csv(epithelial_markers13, 'epithelial_markers13.csv', row.names = TRUE)

# markers14
DefaultAssay(refquery) <- "RNA"
epithelial_markers14 <- FindConservedMarkers(refquery, ident.1 = 14,  grouping.var = "Diagnosis", verbose = T)
head(epithelial_markers14)
write.csv(epithelial_markers14, 'epithelial_markers14.csv', row.names = TRUE)

# markers15
DefaultAssay(refquery) <- "RNA"
epithelial_markers15 <- FindConservedMarkers(refquery, ident.1 = 15,  grouping.var = "Diagnosis", verbose = T)
head(epithelial_markers15)
write.csv(epithelial_markers15, 'epithelial_markers15.csv', row.names = TRUE)

# markers16
DefaultAssay(refquery) <- "RNA"
epithelial_markers16 <- FindConservedMarkers(refquery, ident.1 = 16,  grouping.var = "Diagnosis", verbose = T)
head(epithelial_markers16)
write.csv(epithelial_markers16, 'epithelial_markers16.csv', row.names = TRUE)

# markers17
DefaultAssay(refquery) <- "RNA"
epithelial_markers17 <- FindConservedMarkers(refquery, ident.1 = 17,  grouping.var = "Diagnosis", verbose = T)
head(epithelial_markers17)
write.csv(epithelial_markers17, 'epithelial_markers17.csv', row.names = TRUE)

# markers18
DefaultAssay(refquery) <- "RNA"
epithelial_markers18 <- FindConservedMarkers(refquery, ident.1 = 18,  grouping.var = "Diagnosis", verbose = T)
head(epithelial_markers18)
write.csv(epithelial_markers18, 'epithelial_markers18.csv', row.names = TRUE)

# markers19
DefaultAssay(refquery) <- "RNA"
epithelial_markers19 <- FindConservedMarkers(refquery, ident.1 = 19,  grouping.var = "Diagnosis", verbose = T)
head(epithelial_markers19)
write.csv(epithelial_markers19, 'epithelial_markers19.csv', row.names = TRUE)

# markers20
DefaultAssay(refquery) <- "RNA"
epithelial_markers20 <- FindConservedMarkers(refquery, ident.1 = 20,  grouping.var = "Diagnosis", verbose = T)
head(epithelial_markers20)
write.csv(epithelial_markers20, 'epithelial_markers20.csv', row.names = TRUE)

# markers21
DefaultAssay(refquery) <- "RNA"
epithelial_markers21 <- FindConservedMarkers(refquery, ident.1 = 21,  grouping.var = "Diagnosis", verbose = T)
head(epithelial_markers21)
write.csv(epithelial_markers21, 'epithelial_markers21.csv', row.names = TRUE)
