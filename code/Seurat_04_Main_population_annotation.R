library(Seurat)
library(here)
library(tidyverse)

refquery <- readRDS(here('explore/RDS/04_Integrate_fill_null.RDS'))
head(refquery)

# epithelial cells = (EPCAM)
# endothelial cells = (PECAM-1 = CD31)
# immune cells = (PTPRC = CD45)
# <Mesenchymal Cells>
# 1) Smooth Muscle Cells - ACTA2, MYH11, PDGFRB
# 2) Mesothelial Cells - WT1, UPK3B
# 3) Myofibroblasts - LUM, PDGFRA, ACTA2, MYLK
# 4) HAS1 High Fibroblasts - LUM, PDGFRA, HAS1, TWIST
# 5) Fibroblasts - LUM, PDGFRA
# 6) PLIN2+ Fibroblasts  - LUM, PDGFRA, PLIN2

# make markers 
features_epithelial <- c("EPCAM")
features_endothelial <- c("PECAM1")
features_immune <- c("PTPRC")
features_mesenchymal_1 <- c("ACTA2", "MYH11", "PDGFRB")
features_mesenchymal_2 <- c("WT1", "UPK3B")
features_mesenchymal_3 <- c("LUM", "PDGFRA", "ACTA2", "MYLK")
features_mesenchymal_4 <- c("LUM", "PDGFRA", "HAS1", "TWIST")
features_mesenchymal_5 <- c("LUM", "PDGFRA")
features_mesenchymal_6 <- c("PLIN2")

# check out expression 
FeaturePlot(refquery, features = features_epithelial)
FeaturePlot(refquery, features = features_endothelial)
FeaturePlot(refquery, features = features_immune)
FeaturePlot(refquery, features = features_mesenchymal_1)
FeaturePlot(refquery, features = features_mesenchymal_2)
FeaturePlot(refquery, features = features_mesenchymal_3)
FeaturePlot(refquery, features = features_mesenchymal_4)
FeaturePlot(refquery, features = features_mesenchymal_5)
FeaturePlot(refquery, features = features_mesenchymal_6)

DimPlot(refquery, reduction = "umap", group.by = "seurat_clusters", label = TRUE, label.size = 3, repel = TRUE) 
DimPlot(refquery, reduction = "umap", group.by = "integrated_snn_res.0.8", label = TRUE, label.size = 3, repel = TRUE)

# population annotation
refquery@meta.data$cell_type_main <- ifelse(refquery$integrated_snn_res.0.8 %in% c(3, 29, 22, 10, 33, 17, 14, 21, 5, 2, 28), 'Epithelial cells', 'NA')
refquery@meta.data$cell_type_main <- ifelse(refquery$integrated_snn_res.0.8 %in% c(25, 15, 16, 20), 'Endothelial cells', refquery$cell_type_main)
refquery@meta.data$cell_type_main <- ifelse(refquery$integrated_snn_res.0.8 %in% c(0, 9, 1, 18, 6, 12, 4, 23, 32, 11, 8, 7, 30, 24, 26, 31, 19), 'Immune cells', refquery$cell_type_main)
refquery@meta.data$cell_type_main <- ifelse(refquery$integrated_snn_res.0.8 %in% c(13, 27, 34), 'Mesenchymal cells', refquery$cell_type_main)

saveRDS(refquery, file = "05_Integrate_df_annotation.RDS")


