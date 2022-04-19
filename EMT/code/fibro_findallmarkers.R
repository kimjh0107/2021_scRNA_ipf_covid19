library(here)
library(Seurat)
library(tidyverse)
library(ggpubr)


# Sup_Fig1a 
df <- readRDS(here("EMT/RDS/fibro.RDS"))
DefaultAssay(df) <- "RNA"
Idents(df) <- "seurat_clusters"

unique(df$new_submain)
all_markers <- FindAllMarkers(df, logfc.threshold = 0.5, min.pct = 0.25, test.use = "wilcox")
all_markers$gene <- rownames(all_markers)
write_csv(all_markers, here("EMT/csv/fibro_clusters_findallmarkers.csv"))





Idents(df) <- "Diagnosis"
DefaultAssay(df) <- "RNA"

# Upregulated
df_covid_control <- FindMarkers(df, ident.1 = "COVID", ident.2 = "Control", only.pos = TRUE, min.pct = 0.25,logfc.threshold = 0, test.use = "wilcox")
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