library(Seurat)
library(here)
library(tidyverse)
library(cowplot)
library(knitr)
library(sctransform)
library(Matrix)
library(future)
library(future.apply)

plan('multiprocess', workers = 12)
options(future.globals.maxSize = 2e+05 * 1024^2)


refquery <- readRDS(here('explore/RDS/04_Integrate_fill_null.RDS'))
#head(refquery)




DefaultAssay(refquery) <- "RNA"
markers0 <- FindConservedMarkers(refquery, ident.1 = 0,  grouping.var = "Diagnosis", verbose = T)
head(markers0)
write.csv(markers0, 'markers0.csv', row.names = TRUE)



DefaultAssay(refquery) <- "RNA"
markers1 <- FindConservedMarkers(refquery, ident.1 = 1,  grouping.var = "Diagnosis", verbose = T)
head(markers1)
write.csv(markers1, 'markers1.csv', row.names = TRUE)



DefaultAssay(refquery) <- "RNA"
markers2 <- FindConservedMarkers(refquery, ident.1 = 2,  grouping.var = "Diagnosis", verbose = T)
head(markers2)
write.csv(markers2, 'markers2.csv', row.names = TRUE)



DefaultAssay(refquery) <- "RNA"
markers3 <- FindConservedMarkers(refquery, ident.1 = 3,  grouping.var = "Diagnosis", verbose = T)
head(markers3)
write.csv(markers3, 'markers3.csv', row.names = TRUE)



DefaultAssay(refquery) <- "RNA"
markers4 <- FindConservedMarkers(refquery, ident.1 = 4,  grouping.var = "Diagnosis", verbose = T)
head(markers4)
write.csv(markers4, 'markers4.csv', row.names = TRUE)



DefaultAssay(refquery) <- "RNA"
markers5 <- FindConservedMarkers(refquery, ident.1 = 5,  grouping.var = "Diagnosis", verbose = T)
head(markers5)
write.csv(markers5, 'markers0.csv', row.names = TRUE)



DefaultAssay(refquery) <- "RNA"
markers6 <- FindConservedMarkers(refquery, ident.1 = 6,  grouping.var = "Diagnosis", verbose = T)
head(markers6)
write.csv(markers6, 'markers0.csv', row.names = TRUE)



DefaultAssay(refquery) <- "RNA"
markers7 <- FindConservedMarkers(refquery, ident.1 = 7,  grouping.var = "Diagnosis", verbose = T)
head(markers7)
write.csv(markers0, 'markers7.csv', row.names = TRUE)



DefaultAssay(refquery) <- "RNA"
markers8 <- FindConservedMarkers(refquery, ident.1 = 8,  grouping.var = "Diagnosis", verbose = T)
head(markers8)
write.csv(markers8, 'markers8.csv', row.names = TRUE)



DefaultAssay(refquery) <- "RNA"
markers9 <- FindConservedMarkers(refquery, ident.1 = 9,  grouping.var = "Diagnosis", verbose = T)
head(markers9)
write.csv(markers9, 'markers9.csv', row.names = TRUE)



DefaultAssay(refquery) <- "RNA"
markers10 <- FindConservedMarkers(refquery, ident.1 = 10,  grouping.var = "Diagnosis", verbose = T)
head(markers10)
write.csv(markers10, 'markers10.csv', row.names = TRUE)



DefaultAssay(refquery) <- "RNA"
markers11 <- FindConservedMarkers(refquery, ident.1 = 11,  grouping.var = "Diagnosis", verbose = T)
head(markers11)
write.csv(markers11, 'markers11.csv', row.names = TRUE)



DefaultAssay(refquery) <- "RNA"
markers12 <- FindConservedMarkers(refquery, ident.1 = 12,  grouping.var = "Diagnosis", verbose = T)
head(markers12)
write.csv(markers12, 'markers12.csv', row.names = TRUE)



DefaultAssay(refquery) <- "RNA"
markers13 <- FindConservedMarkers(refquery, ident.1 = 13,  grouping.var = "Diagnosis", verbose = T)
head(markers13)
write.csv(markers13, 'markers0.csv', row.names = TRUE)



DefaultAssay(refquery) <- "RNA"
markers14 <- FindConservedMarkers(refquery, ident.1 = 14,  grouping.var = "Diagnosis", verbose = T)
head(markers14)
write.csv(markers14, 'markers14.csv', row.names = TRUE)



DefaultAssay(refquery) <- "RNA"
markers15 <- FindConservedMarkers(refquery, ident.1 = 15,  grouping.var = "Diagnosis", verbose = T)
head(markers15)
write.csv(markers0, 'markers15.csv', row.names = TRUE)



DefaultAssay(refquery) <- "RNA"
markers16 <- FindConservedMarkers(refquery, ident.1 = 16,  grouping.var = "Diagnosis", verbose = T)
head(markers16)
write.csv(markers16, 'markers16.csv', row.names = TRUE)



DefaultAssay(refquery) <- "RNA"
markers17 <- FindConservedMarkers(refquery, ident.1 = 17,  grouping.var = "Diagnosis", verbose = T)
head(markers17)
write.csv(markers17, 'markers17.csv', row.names = TRUE)



DefaultAssay(refquery) <- "RNA"
markers18 <- FindConservedMarkers(refquery, ident.1 = 18,  grouping.var = "Diagnosis", verbose = T)
head(markers18)
write.csv(markers18, 'markers18.csv', row.names = TRUE)



DefaultAssay(refquery) <- "RNA"
markers19 <- FindConservedMarkers(refquery, ident.1 = 19,  grouping.var = "Diagnosis", verbose = T)
head(markers19)
write.csv(markers19, 'markers19.csv', row.names = TRUE)



DefaultAssay(refquery) <- "RNA"
markers20 <- FindConservedMarkers(refquery, ident.1 = 20,  grouping.var = "Diagnosis", verbose = T)
head(markers20)
write.csv(markers20, 'markers20.csv', row.names = TRUE)



DefaultAssay(refquery) <- "RNA"
markers21 <- FindConservedMarkers(refquery, ident.1 = 21,  grouping.var = "Diagnosis", verbose = T)
head(markers21)
write.csv(markers21, 'markers21.csv', row.names = TRUE)



DefaultAssay(refquery) <- "RNA"
markers22 <- FindConservedMarkers(refquery, ident.1 = 22,  grouping.var = "Diagnosis", verbose = T)
head(markers22)
write.csv(markers22, 'markers22.csv', row.names = TRUE)



DefaultAssay(refquery) <- "RNA"
markers23 <- FindConservedMarkers(refquery, ident.1 = 23,  grouping.var = "Diagnosis", verbose = T)
head(markers23)
write.csv(markers23, 'markers23.csv', row.names = TRUE)



DefaultAssay(refquery) <- "RNA"
markers24 <- FindConservedMarkers(refquery, ident.1 = 24,  grouping.var = "Diagnosis", verbose = T)
head(markers24)
write.csv(markers24, 'markers24.csv', row.names = TRUE)



DefaultAssay(refquery) <- "RNA"
markers25 <- FindConservedMarkers(refquery, ident.1 = 25,  grouping.var = "Diagnosis", verbose = T)
head(markers25)
write.csv(markers25, 'markers25.csv', row.names = TRUE)



DefaultAssay(refquery) <- "RNA"
markers26 <- FindConservedMarkers(refquery, ident.1 = 26,  grouping.var = "Diagnosis", verbose = T)
head(markers26)
write.csv(markers26, 'markers26.csv', row.names = TRUE)



DefaultAssay(refquery) <- "RNA"
markers27 <- FindConservedMarkers(refquery, ident.1 = 27,  grouping.var = "Diagnosis", verbose = T)
head(markers27)
write.csv(markers27, 'markers27.csv', row.names = TRUE)



DefaultAssay(refquery) <- "RNA"
markers28 <- FindConservedMarkers(refquery, ident.1 = 28,  grouping.var = "Diagnosis", verbose = T)
head(markers28)
write.csv(markers28, 'markers28.csv', row.names = TRUE)



DefaultAssay(refquery) <- "RNA"
markers29 <- FindConservedMarkers(refquery, ident.1 = 29,  grouping.var = "Diagnosis", verbose = T)
head(markers29)
write.csv(markers29, 'markers29.csv', row.names = TRUE)



DefaultAssay(refquery) <- "RNA"
markers30 <- FindConservedMarkers(refquery, ident.1 = 30,  grouping.var = "Diagnosis", verbose = T)
head(markers30)
write.csv(markers30, 'markers30.csv', row.names = TRUE)



DefaultAssay(refquery) <- "RNA"
markers31 <- FindConservedMarkers(refquery, ident.1 = 31,  grouping.var = "Diagnosis", verbose = T)
head(markers31)
write.csv(markers31, 'markers31.csv', row.names = TRUE)


DefaultAssay(refquery) <- "RNA"
markers33 <- FindConservedMarkers(refquery, ident.1 = 33,  grouping.var = "Diagnosis", verbose = T)
head(markers33)
write.csv(markers33, 'markers33.csv', row.names = TRUE)



DefaultAssay(refquery) <- "RNA"
markers34 <- FindConservedMarkers(refquery, ident.1 = 34,  grouping.var = "Diagnosis", verbose = T)
head(markers34)
write.csv(markers34, 'markers34.csv', row.names = TRUE)



