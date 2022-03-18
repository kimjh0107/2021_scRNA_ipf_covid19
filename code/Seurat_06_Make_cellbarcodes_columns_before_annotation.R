library(Seurat)
library(here)
library(tidyverse)

refquery <- readRDS(here('explore/RDS/05_Integrate_df_annotation.RDS'))

# make new cellbarcodes columns 
refquery@meta.data$cellbarcodes <- rownames(refquery@meta.data)

saveRDS(refquery, file = "05_New_Integrate_df.RDS")
