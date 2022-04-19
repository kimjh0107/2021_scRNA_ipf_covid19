library(here)
library(Seurat)
library(tidyverse)
library(ggpubr)

df <- readRDS(here("EMT/RDS/snRNA.RDS"))
df <- subset(df, subset = Diagnosis %in% c("Control", "COVID-19"))

unique(df$Diagnosis_tag)
DefaultAssay(df) <- "RNA"
p4 <- DimPlot(df, group.by = "cell_type_main", label = T, repel = T) + NoLegend()
p1 <- FeaturePlot(df, features = c('EPCAM'), min.cutoff = 'q10', max.cutoff = 'q90', order = T, raster = T)
p2 <- FeaturePlot(df, features = c('PECAM1'), min.cutoff = 'q10', max.cutoff = 'q90', order = T, raster = T)
p3 <- FeaturePlot(df, features = c('PTPRC'), min.cutoff = 'q10', max.cutoff = 'q90', order = T, raster = T)

ggarrange(p4,p1,p2,p3, ncol = 4, nrow = 1)
ggsave(here("EMT/figure/snRNA_DimPlot_celltypemain_markers.pdf"), width = 20, height = 5)



# EMT부분을 위한 새로운 figure 추가해주도록 하기 
df <- readRDS(here("EMT/RDS/snRNA_EMT.RDS"))
DefaultAssay(df) <- "RNA"
p1 <- DimPlot(df, group.by = "cell_type_submain", label = T, repel = T) + NoLegend()
p2 <- DimPlot(df, group.by = "Diagnosis", label = F, repel = T, shuffle = T, cols = c("#00A0FA", "red", "purple")) + ggplot2::theme(legend.position = "bottom")

# fibro markers 
p3 <- FeaturePlot(df, features = c('COL1A1'), min.cutoff = 'q10', max.cutoff = 'q90', order = T, raster = T)
p4 <- FeaturePlot(df, features = c('PDGFRA'), min.cutoff = 'q10', max.cutoff = 'q90', order = T, raster = T)

p5 <- FeaturePlot(df, features = c('SFTPA1'), min.cutoff = 'q10', max.cutoff = 'q90', order = T, raster = T)
p6 <- FeaturePlot(df, features = c('SFTPC'), min.cutoff = 'q10', max.cutoff = 'q90', order = T, raster = T)

p7 <- FeaturePlot(df, features = c('AGER'), min.cutoff = 'q10', max.cutoff = 'q90', order = T, raster = T)
p8 <- FeaturePlot(df, features = c('CLIC5'), min.cutoff = 'q10', max.cutoff = 'q90', order = T, raster = T)


ggarrange(p1,p2,p3,p4, ncol = 4, nrow = 1)
ggsave(here("EMT/figure/snRNA_DimPlot_EMT_and_featureplot.pdf"), width = 20, height = 5)

ggarrange(p5,p6,p7,p8, ncol = 4, nrow = 1)
ggsave(here("EMT/figure/snRNA_Featureplot_EMT_markers.pdf"), width = 20, height = 5)


