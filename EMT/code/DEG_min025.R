library(here)
library(Seurat)
library(tidyverse)
library(ggpubr)


# Sup_Fig1a 
df <- readRDS(here("EMT/RDS/EMT.RDS"))

at1 <- subset(df, subset = new_submain == "AT1")
at2 <- subset(df, subset = new_submain == "AT2")
fibro <- subset(df, subset = new_submain == "Fibroblasts")

Idents(at1) <- "Diagnosis"
Idents(at2) <- "Diagnosis"
Idents(fibro) <- "Diagnosis"

DefaultAssay(at1) <- "RNA"
DefaultAssay(at2) <- "RNA"
DefaultAssay(fibro) <- "RNA"

#### fibro
# Upregulated
df_covid_control <- FindMarkers(fibro, ident.1 = "COVID", ident.2 = "Control", only.pos = TRUE, min.pct = 0.25,logfc.threshold = 0.25, test.use = "wilcox")
df_covid_control$gene <- rownames(df_covid_control)
write.csv(df_covid_control, here("EMT/csv/fibro_upregulated_covid_min0.csv"))

df_ipf_control <- FindMarkers(df, ident.1 = "IPF", ident.2 = "Control", only.pos = TRUE, min.pct = 0.25,logfc.threshold = 0, test.use = "wilcox")
df_ipf_control$gene <- rownames(df_ipf_control)
write.csv(df_ipf_control, here("EMT/csv/fibro_upregulated_ipf_min0.csv"))

df_upregulated_merge <- inner_join(df_covid_control, df_ipf_control, by = "gene")
write.csv(df_upregulated_merge, here("EMT/csv/fibro_upregulated_merge_min0.csv"))



# Downregulated_min0
df_covid_control <- FindMarkers(df, ident.1 = "Control", ident.2 = "COVID", only.pos = TRUE, min.pct = 0.25,logfc.threshold = 0, test.use = "wilcox")
df_covid_control$gene <- rownames(df_covid_control)
write.csv(df_covid_control, here("EMT/csv/fibro_downregulated_covid_min0.csv"))

df_ipf_control <- FindMarkers(df, ident.1 = "Control", ident.2 = "IPF", only.pos = TRUE, min.pct = 0.25,logfc.threshold = 0, test.use = "wilcox")
df_ipf_control$gene <- rownames(df_ipf_control)
write.csv(df_ipf_control, here("EMT/csv/fibro_downregulated_ipf_min0.csv"))

df_upregulated_merge <- inner_join(df_covid_control, df_ipf_control, by = "gene")
write.csv(df_upregulated_merge, here("EMT/csv/fibroell_downregulated_merge_min0.csv"))