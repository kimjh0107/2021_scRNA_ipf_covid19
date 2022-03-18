library(Seurat)
library(here)
library(tidyverse)

data <- Read10X_h5(here('data/GSE149878_covid_lung/GSM4516279_C166_filtered_feature_bc_matrix.h5'), use.names = TRUE, unique.features = TRUE)
obj1 <- CreateSeuratObject(data)

# C168
data <- Read10X_h5(here('data/GSE149878_covid_lung/GSM4516280_C168_filtered_feature_bc_matrix.h5'), use.names = TRUE, unique.features = TRUE)
obj2 <- CreateSeuratObject(data)

# C170
data <- Read10X_h5(here('data/GSE149878_covid_lung/GSM4516281_C170_filtered_feature_bc_matrix.h5'), use.names = TRUE, unique.features = TRUE)
obj3 <- CreateSeuratObject(data)

# C172
data <- Read10X_h5(here('data/GSE149878_covid_lung/GSM4516282_C172_filtered_feature_bc_matrix.h5'), use.names = TRUE, unique.features = TRUE)
obj4 <- CreateSeuratObject(data)

# merge SeuratObject
covid <- merge(obj1, y = c(obj2, obj3, obj4), add.cell.ids = c('obj1', 'obj2', 'obj3', 'obj4'), project = 'scRNA_IPF_data')
