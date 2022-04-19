library(Seurat)
library(here)
library(tidyverse)
library(biomaRt)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DoMultiBarHeatmap)

# EMT subset cell heatmap 
df <- readRDS(here("EMT/RDS/EMT.RDS"))
DefaultAssay(df) <- "RNA"
df <- ScaleData(df)
Idents(df) <- "new_submain"
levels(df) <- c("AT1", "AT2", "Fibroblasts")
levels(df)
celltype_markers <- read_csv(here("EMT/csv/EMT_celltype_findallmarkers.csv"))
cell_type <- arrange(celltype_markers, cluster, desc(avg_log2FC))
cell_type <- cell_type %>% filter(avg_log2FC > 1.5)
cell_type
DoHeatmap(df, features = cell_type$gene)
ggsave(here("EMT/figure/EMT_heatmap_rna.pdf"), width = 15, height = 15)






# Celltype heatmap 
df <- readRDS(here("EMT/RDS/df.RDS"))
DefaultAssay(df) <- "RNA"
df <- ScaleData(df)
Idents(df) <- "new_submain"
celltype_markers <- read_csv(here("EMT/csv/celltype_findallmarkers.csv"))
cell_type <- arrange(celltype_markers, cluster, desc(avg_log2FC))
cell_type <- cell_type %>% filter(avg_log2FC > 2)
cell_type
DoHeatmap(df, features = cell_type$gene)
ggsave(here("EMT/figure/Heatmap//celltype_heatmap_rna.pdf"), width = 20, height = 20)




# Fibro DGE heatmap 
fibro_up_merge <- read_csv(here("EMT/csv/fibro_upregulated_merge_min0.csv"))
fibro_up_merge <- fibro_up_merge %>% filter(p_val_adj.x < 0.05 & p_val_adj.y < 0.05 & avg_log2FC.x > 0.25 & avg_log2FC.y > 0.25)
fibro_down_merge <- read_csv(here("EMT/csv/fibroell_downregulated_merge_min0.csv"))
fibro_down_merge <- fibro_down_merge %>% filter(p_val_adj.x < 0.05 & p_val_adj.y < 0.05 & avg_log2FC.x > 0.25 & avg_log2FC.y > 0.25)

fibro <- readRDS(here("EMT/RDS/fibro.RDS"))
DefaultAssay(fibro) <- "RNA"
fibro <- ScaleData(fibro)

Idents(fibro) <- "Diagnosis"

# heatmap 
DoHeatmap(fibro, features = fibro_up_merge$gene)
ggsave(here("EMT/figure/Heatmap/heatmap_fibro_up_merge.pdf"), width = 20, height = 20)
DoHeatmap(fibro, features = fibro_down_merge$gene)
ggsave(here("EMT/figure/Heatmap/heatmap_fibro_down_merge.pdf"), width = 20, height = 20)



# AT2 & AT1 
# AT2
at2_up_merge <- read_csv(here("EMT/csv/at2_upregulated_merge_min0.csv"))
at2_up_merge <- at2_up_merge %>% filter(p_val_adj.x < 0.05 & p_val_adj.y < 0.05 & avg_log2FC.x > 0.25 & avg_log2FC.y > 0.25)
at2_down_merge <- read_csv(here("EMT/csv/at2_downregulated_merge_min0.csv"))
at2_down_merge <- at2_down_merge %>% filter(p_val_adj.x < 0.05 & p_val_adj.y < 0.05 & avg_log2FC.x > 0.25 & avg_log2FC.y > 0.25)

# AT1
at1_up_merge <- read_csv(here("EMT/csv/at1_upregulated_merge_min0.csv"))
at1_up_merge <- at1_up_merge %>% filter(p_val_adj.x < 0.05 & p_val_adj.y < 0.05 & avg_log2FC.x > 0.25 & avg_log2FC.y > 0.25)
at1_down_merge <- read_csv(here("EMT/csv/at1_downregulated_merge_min0.csv"))
at1_down_merge <- at1_down_merge %>% filter(p_val_adj.x < 0.05 & p_val_adj.y < 0.05 & avg_log2FC.x > 0.25 & avg_log2FC.y > 0.25)

at <- readRDS(here("EMT/RDS/AT.RDS"))
at2 <- subset(at, subset = new_submain == "AT2")
at1 <- subset(at, subset = new_submain == "AT1")
DefaultAssay(at2) <- "RNA"
DefaultAssay(at1) <- "RNA"
at2 <- ScaleData(at2)
at1 <- ScaleData(at1)

Idents(at2) <- "Diagnosis"
Idents(at1) <- "Diagnosis"

# heatmap 
# AT2 
DoHeatmap(at2, features = at2_up_merge$gene)
ggsave(here("EMT/figure/Heatmap/heatmap_at2_up_merge.pdf"), width = 20, height = 20)
DoHeatmap(at2, features = at2_down_merge$gene)
ggsave(here("EMT/figure/Heatmap/heatmap_at2_down_merge.pdf"), width = 20, height = 20)

# AT1
DoHeatmap(at1, features = at1_up_merge$gene)
ggsave(here("EMT/figure/Heatmap/heatmap_at1_up_merge.pdf"), width = 20, height = 20)
DoHeatmap(at1, features = at1_down_merge$gene)
ggsave(here("EMT/figure/Heatmap/heatmap_at1_down_merge.pdf"), width = 20, height = 20)

