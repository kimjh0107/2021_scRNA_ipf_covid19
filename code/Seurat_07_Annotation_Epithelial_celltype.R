library(Seurat)
library(here)
library(tidyverse)

refquery <- readRDS(here('explore/RDS/05_New_Integrate_df.RDS'))

# Epithelial cells 
refquery <- subset(refquery, cell_type_main == 'Epithelial cells')
refquery <- ScaleData(object = refquery)
refquery <- RunPCA(object = refquery)
refquery <- FindNeighbors(refquery, dims = 1:20)
refquery <- FindClusters(refquery, resolution = 0.6)
refquery <- RunUMAP(object = refquery, dims = 1:20)

refquery@meta.data$cell_type_submain <- ifelse(refquery@meta.data$integrated_snn_res.0.6 %in% c(11, 13), "AT1", "NA")
refquery@meta.data$cell_type_submain <- ifelse(refquery@meta.data$integrated_snn_res.0.6 %in% c(10, 5, 8, 7, 14, 15), "AT2", refquery@meta.data$cell_type_submain)
refquery@meta.data$cell_type_submain <- ifelse(refquery@meta.data$integrated_snn_res.0.6 %in% c(6), "Basal Cell", refquery@meta.data$cell_type_submain)
refquery@meta.data$cell_type_submain <- ifelse(refquery@meta.data$integrated_snn_res.0.6 %in% c(1, 2, 0, 20, 18, 12, 9, 21, 19), "Ciliated Cell", refquery@meta.data$cell_type_submain)
refquery@meta.data$cell_type_submain <- ifelse(refquery@meta.data$integrated_snn_res.0.6 %in% c(3, 17, 16), "Club Cell", refquery@meta.data$cell_type_submain)
refquery@meta.data$cell_type_submain <- ifelse(refquery@meta.data$integrated_snn_res.0.6 %in% c(4), "Gobelt Cel", refquery@meta.data$cell_type_submain)

saveRDS(refquery, file = "06_Epithelial_celltype.RDS" )
celltype <- refquery@meta.data %>% select("cellbarcodes", "cell_type_submain")
write.csv(celltype, "New_Epithelial_celltype")