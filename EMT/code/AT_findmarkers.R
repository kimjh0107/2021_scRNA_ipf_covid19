library(here)
library(Seurat)
library(tidyverse)
library(ggpubr)

df <- readRDS(here("EMT/RDS/EMT.RDS"))
df <- subset(df, subset = new_submain %in% c("AT2"))


Idents(df) <- "Diagnosis"
DefaultAssay(df) <- "RNA"

# Upregulated
df_covid_control <- FindMarkers(df, ident.1 = "COVID", ident.2 = "Control", only.pos = TRUE, min.pct = 0.25,logfc.threshold = 0, test.use = "wilcox")
df_covid_control$gene <- rownames(df_covid_control)
write.csv(df_covid_control, here("EMT/csv/at2_upregulated_covid_min0.csv"))

df_ipf_control <- FindMarkers(df, ident.1 = "IPF", ident.2 = "Control", only.pos = TRUE, min.pct = 0.25,logfc.threshold = 0, test.use = "wilcox")
df_ipf_control$gene <- rownames(df_ipf_control)
write.csv(df_ipf_control, here("EMT/csv/at2_upregulated_ipf_min0.csv"))

df_upregulated_merge <- inner_join(df_covid_control, df_ipf_control, by = "gene")
write.csv(df_upregulated_merge, here("EMT/csv/at2_upregulated_merge_min0.csv"))



# Downregulated_min0
df_covid_control <- FindMarkers(df, ident.1 = "Control", ident.2 = "COVID", only.pos = TRUE, min.pct = 0.25,logfc.threshold = 0, test.use = "wilcox")
df_covid_control$gene <- rownames(df_covid_control)
write.csv(df_covid_control, here("EMT/csv/at2_downregulated_covid_min0.csv"))

df_ipf_control <- FindMarkers(df, ident.1 = "Control", ident.2 = "IPF", only.pos = TRUE, min.pct = 0.25,logfc.threshold = 0, test.use = "wilcox")
df_ipf_control$gene <- rownames(df_ipf_control)
write.csv(df_ipf_control, here("EMT/csv/at2_downregulated_ipf_min0.csv"))

df_upregulated_merge <- inner_join(df_covid_control, df_ipf_control, by = "gene")
write.csv(df_upregulated_merge, here("EMT/csv/at2_downregulated_merge_min0.csv"))








df <- readRDS(here("EMT/RDS/EMT.RDS"))
df <- subset(df, subset = new_submain %in% c("AT1"))


Idents(df) <- "Diagnosis"
DefaultAssay(df) <- "RNA"

# Upregulated
df_covid_control <- FindMarkers(df, ident.1 = "COVID", ident.2 = "Control", only.pos = TRUE, min.pct = 0.25,logfc.threshold = 0, test.use = "wilcox")
df_covid_control$gene <- rownames(df_covid_control)
write.csv(df_covid_control, here("EMT/csv/at1_upregulated_covid_min0.csv"))

df_ipf_control <- FindMarkers(df, ident.1 = "IPF", ident.2 = "Control", only.pos = TRUE, min.pct = 0.25,logfc.threshold = 0, test.use = "wilcox")
df_ipf_control$gene <- rownames(df_ipf_control)
write.csv(df_ipf_control, here("EMT/csv/at1_upregulated_ipf_min0.csv"))

df_upregulated_merge <- inner_join(df_covid_control, df_ipf_control, by = "gene")
write.csv(df_upregulated_merge, here("EMT/csv/at1_upregulated_merge_min0.csv"))



# Downregulated_min0
df_covid_control <- FindMarkers(df, ident.1 = "Control", ident.2 = "COVID", only.pos = TRUE, min.pct = 0.25,logfc.threshold = 0, test.use = "wilcox")
df_covid_control$gene <- rownames(df_covid_control)
write.csv(df_covid_control, here("EMT/csv/at1_downregulated_covid_min0.csv"))

df_ipf_control <- FindMarkers(df, ident.1 = "Control", ident.2 = "IPF", only.pos = TRUE, min.pct = 0.25,logfc.threshold = 0, test.use = "wilcox")
df_ipf_control$gene <- rownames(df_ipf_control)
write.csv(df_ipf_control, here("EMT/csv/at1_downregulated_ipf_min0.csv"))

df_upregulated_merge <- inner_join(df_covid_control, df_ipf_control, by = "gene")
write.csv(df_upregulated_merge, here("EMT/csv/at1_downregulated_merge_min0.csv"))