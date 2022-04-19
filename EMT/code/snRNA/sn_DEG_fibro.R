library(here)
library(Seurat)
library(tidyverse)
# DGE - fibro 
emt <- readRDS(here("EMT/RDS/snRNA_EMT.RDS"))
df <- subset(emt, subset = cell_type_submain == "Fibroblasts")

# Sup_Fig1a 
DefaultAssay(df) <- "RNA"
Idents(df) <- "Diagnosis"
# Upregulated
df_covid_control <- FindMarkers(df, ident.1 = "COVID-19", ident.2 = "Control", only.pos = TRUE, min.pct = 0.25,logfc.threshold = 0, test.use = "wilcox")
df_covid_control$gene <- rownames(df_covid_control)
write.csv(df_covid_control, here("EMT/csv/snRNA_fibro_upregulated_covid_min0.csv"))

# Downregulated_min0
df_covid_control <- FindMarkers(df, ident.1 = "Control", ident.2 = "COVID-19", only.pos = TRUE, min.pct = 0.25,logfc.threshold = 0, test.use = "wilcox")
df_covid_control$gene <- rownames(df_covid_control)
write.csv(df_covid_control, here("EMT/csv/snRNA_fibro_downregulated_covid_min0.csv"))