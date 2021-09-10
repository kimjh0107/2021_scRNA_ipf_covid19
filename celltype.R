library(Seurat)
library(here)
library(tidyverse)
library(cowplot)
library(knitr)






refquery <- readRDS(here('explore/RDS/05_Integrate_df_annotation.RDS'))

Immune <- read.csv(here('explore/RDS/celltype/Immune_celltype.csv'))
Endothelial <- read.csv(here('explore/RDS/celltype/Endothelial_celltype.csv'))
Epithelial <- read.csv(here('explore/RDS/celltype/Epithelial_celltype'))
Mesenchymal <- read.csv(here('explore/RDS/celltype/Mesenchymal_celltype.csv'))



celltypes <- rbind.data.frame(Immune, Endothelial,Epithelial ,Mesenchymal)
refquery@meta.data$barcode <- rownames(refquery@meta.data)
refquery@meta.data <- left_join(refquery@meta.data, celltypes, by = 'cellbarcodes')

rownames(refquery@meta.data) <- refquery@meta.data$cellbarcodes
refquery@meta.data$cell_type_submain <- ifelse(is.na(refquery@meta.data$cell_type_submain) == T, refquery@meta.data$cell_type_main, refquery@meta.data$cell_type_submain)


saveRDS(total, file = "total_celltype.RDS")


