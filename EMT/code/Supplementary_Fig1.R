library(here)
library(Seurat)
library(tidyverse)
library(ggpubr)


# Sup_Fig1a 
df <- readRDS(here("EMT/RDS/df.RDS"))
DefaultAssay(df) <- "RNA"

p4 <- DimPlot(df, group.by = "cell_type_main", label = T, repel = T) + NoLegend()
#ggsave(here("EMT/figure/Supplementary_Fig1//Fig1a_DimPlot_celltype_main.pdf"), width = 7, height = 5)

# Sup_Fig1a 
# FeaturePlot 
DefaultAssay(df) <- "RNA"

# FeaturePlot 
p1 <- FeaturePlot(df, features = c('EPCAM'), min.cutoff = 'q10', max.cutoff = 'q90', order = T, raster = T)
#ggsave(here("final/figure/FeaturePlot_macro_FABP4.pdf"), width = 7, height = 5)

p2 <- FeaturePlot(df, features = c('PECAM1'), min.cutoff = 'q10', max.cutoff = 'q90', order = T, raster = T)
#ggsave(here("final/figure/FeaturePlot_macro_SPP1.pdf"), width = 7, height = 5)

p3 <- FeaturePlot(df, features = c('PTPRC'), min.cutoff = 'q10', max.cutoff = 'q90', order = T, raster = T)
#ggsave(here("final/figure/FeaturePlot_macro_SPP1.pdf"), width = 7, height = 5)

#ggsave(here("final/figure/FeaturePlot_macro_SPP1.pdf"), width = 7, height = 5)


ggarrange(p4,p1,p2,p3, ncol = 4, nrow = 1)
ggsave(here("EMT/figure/Supplementary_Fig1/Fig1ab.pdf"), width = 20, height = 5)


DimPlot(df, group.by = "new_submain", split.by = "Diagnosis", label = T, repel = T) + NoLegend()
ggsave(here("EMT/figure/Supplementary_Fig1//Fig1c_split_dimplot.pdf"), width = 20, height = 5)




# EMT부분을 위한 새로운 figure 추가해주도록 하기 
df <- readRDS(here("EMT/RDS/EMT.RDS"))
DefaultAssay(df) <- "RNA"
p1 <- DimPlot(df, group.by = "new_submain", label = T, repel = T) + NoLegend()
p2 <- DimPlot(df, group.by = "Diagnosis", label = F, repel = T, shuffle = T, cols = c("#00A0FA", "red", "purple")) + ggplot2::theme(legend.position = "bottom")

# fibro markers 
p3 <- FeaturePlot(df, features = c('COL1A1'), min.cutoff = 'q10', max.cutoff = 'q90', order = T, raster = T)
p4 <- FeaturePlot(df, features = c('PDGFRA'), min.cutoff = 'q10', max.cutoff = 'q90', order = T, raster = T)

p5 <- FeaturePlot(df, features = c('SFTPA1'), min.cutoff = 'q10', max.cutoff = 'q90', order = T, raster = T)
p6 <- FeaturePlot(df, features = c('SFTPC'), min.cutoff = 'q10', max.cutoff = 'q90', order = T, raster = T)

p7 <- FeaturePlot(df, features = c('AGER'), min.cutoff = 'q10', max.cutoff = 'q90', order = T, raster = T)
p8 <- FeaturePlot(df, features = c('CLIC5'), min.cutoff = 'q10', max.cutoff = 'q90', order = T, raster = T)


ggarrange(p1,p2,p3,p4, ncol = 4, nrow = 1)
ggsave(here("EMT/figure/Supplementary_Fig1/DimPlot_EMT_and_featureplot.pdf"), width = 20, height = 5)

ggarrange(p5,p6,p7,p8, ncol = 4, nrow = 1)
ggsave(here("EMT/figure/Supplementary_Fig1/Featureplot_EMT_markers.pdf"), width = 20, height = 5)


df <- readRDS(here("EMT/RDS/fibro.RDS"))
DefaultAssay(df) <- "RNA"
DimPlot(df, group.by = "Diagnosis", label = F, repel = T, shuffle = T, cols = c("#00A0FA", "red", "purple")) + ggplot2::theme(legend.position = "bottom")
ggsave(here("EMT/figure/DimPlot_fibro_by_diagnosis.pdf"), width = 7, height = 5)
