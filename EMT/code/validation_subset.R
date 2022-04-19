library(here)
library(Seurat)
library(tidyverse)

snRNA <- readRDS(here("final_before/RDS/snRNA.RDS"))
unique(snRNA$Diagnosis_tag)
df <- subset(snRNA, subset = Diagnosis %in% c("Control", "COVID-19"))
# subset into EMT, AT, fibro cell types 

# EMT
subset <- subset(df, subset = cell_type_submain %in% c("AT1", "AT2", "Fibroblasts"))
subset <- ScaleData(object = subset)
subset <- RunPCA(object = subset)
subset <- FindNeighbors(subset, dims = 1:30)
subset <- FindClusters(subset, resolution = 0.6)
subset <- RunUMAP(object = subset, dims = 1:30)
saveRDS(subset, file = here("EMT/RDS/val_EMT.RDS"))

# AT
subset <- subset(df, subset = cell_type_submain %in% c("AT1", "AT2"))
subset <- ScaleData(object = subset)
subset <- RunPCA(object = subset)
subset <- FindNeighbors(subset, dims = 1:30)
subset <- FindClusters(subset, resolution = 0.6)
subset <- RunUMAP(object = subset, dims = 1:30)
saveRDS(subset, file = here("EMT/RDS/val_AT.RDS"))

### 1
subset <- subset(df, subset = cell_type_submain %in% c("Fibroblasts"))
subset <- ScaleData(object = subset)
subset <- RunPCA(object = subset)
subset <- FindNeighbors(subset, dims = 1:30)
subset <- FindClusters(subset, resolution = 0.6)
subset <- RunUMAP(object = subset, dims = 1:30)
saveRDS(subset, file = here("EMT/RDS/val_fibro.RDS"))




library(here)
library(tidyverse)
library(Seurat)
scrna <- read_csv(here("EMT/csv/fibro_upregulated_covid_min0.csv"))
snrna <- read_csv(here("EMT/csv/fibro_upregulated_val_covid_min0.csv"))


scrna <- scrna %>% filter(p_val_adj < 0.05, avg_log2FC > 0.25)
snrna <- snnra %>% filter(p_val_adj < 0.05, avg_log2FC > 0.25)
scrna
snrna
inner <- inner_join(scrna, snrna, by = "gene")
inner
write_csv(inner,here("EMT/csv/merge_fibro_scrna_snrna_result.csv"))
