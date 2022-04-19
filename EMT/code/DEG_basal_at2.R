library(here)
library(Seurat)
library(tidyverse)
library(ggpubr)


# Sup_Fig1a 
df <- readRDS(here("EMT/RDS/df.RDS"))

basal <- subset(df, subset = new_submain == "Basal cells")
at2 <- subset(df, subset = new_submain == "AT2")
Idents(basal) <- "Diagnosis"
Idents(at2) <- "Diagnosis"

DefaultAssay(basal) <- "RNA"
DefaultAssay(at2) <- "RNA"


#### basal
# Upregulated
df_covid_control <- FindMarkers(basal, ident.1 = "COVID", ident.2 = "Control", only.pos = TRUE, min.pct = 0.1,logfc.threshold = 0.25, test.use = "wilcox")
df_covid_control$gene <- rownames(df_covid_control)
write.csv(df_covid_control, here("EMT/csv/basal_upregulated_covid_minpct01.csv"))

df_ipf_control <- FindMarkers(basal, ident.1 = "IPF", ident.2 = "Control", only.pos = TRUE, min.pct = 0.1,logfc.threshold = 0.25, test.use = "wilcox")
df_ipf_control$gene <- rownames(df_ipf_control)
write.csv(df_ipf_control, here("EMT/csv/basal_upregulated_ipf_minpct01.csv"))

df_upregulated_merge <- inner_join(df_covid_control, df_ipf_control, by = "gene")
write.csv(df_upregulated_merge, here("EMT/csv/basal_upregulated_merge_minpct01.csv"))



# Downregulated_min0
df_covid_control <- FindMarkers(basal, ident.1 = "Control", ident.2 = "COVID", only.pos = TRUE, min.pct = 0.1,logfc.threshold = 0.25, test.use = "wilcox")
df_covid_control$gene <- rownames(df_covid_control)
write.csv(df_covid_control, here("EMT/csv/basal_downregulated_covid_minpct01.csv"))

df_ipf_control <- FindMarkers(basal, ident.1 = "Control", ident.2 = "IPF", only.pos = TRUE, min.pct = 0.1,logfc.threshold = 0.25, test.use = "wilcox")
df_ipf_control$gene <- rownames(df_ipf_control)
write.csv(df_ipf_control, here("EMT/csv/basal_downregulated_ipf_minpct01.csv"))

df_upregulated_merge <- inner_join(df_covid_control, df_ipf_control, by = "gene")
write.csv(df_upregulated_merge, here("EMT/csv/basal_downregulated_merge_minpct01.csv"))






#### AT2
# Upregulated
df_covid_control <- FindMarkers(at2, ident.1 = "COVID", ident.2 = "Control", only.pos = TRUE, min.pct = 0.1,logfc.threshold = 0.25, test.use = "wilcox")
df_covid_control$gene <- rownames(df_covid_control)
write.csv(df_covid_control, here("EMT/csv/at2_upregulated_covid_minpct01.csv"))

df_ipf_control <- FindMarkers(at2, ident.1 = "IPF", ident.2 = "Control", only.pos = TRUE, min.pct = 0.1,logfc.threshold = 0.25, test.use = "wilcox")
df_ipf_control$gene <- rownames(df_ipf_control)
write.csv(df_ipf_control, here("EMT/csv/at2_upregulated_ipf_minpct01.csv"))

df_upregulated_merge <- inner_join(df_covid_control, df_ipf_control, by = "gene")
write.csv(df_upregulated_merge, here("EMT/csv/at2_upregulated_merge_minpct01.csv"))



# Downregulated_min0
df_covid_control <- FindMarkers(at2, ident.1 = "Control", ident.2 = "COVID", only.pos = TRUE, min.pct = 0.1,logfc.threshold = 0.25, test.use = "wilcox")
df_covid_control$gene <- rownames(df_covid_control)
write.csv(df_covid_control, here("EMT/csv/at2_downregulated_covid_minpct01.csv"))

df_ipf_control <- FindMarkers(at2, ident.1 = "Control", ident.2 = "IPF", only.pos = TRUE, min.pct = 0.1,logfc.threshold = 0.25, test.use = "wilcox")
df_ipf_control$gene <- rownames(df_ipf_control)
write.csv(df_ipf_control, here("EMT/csv/at2_downregulated_ipf_minpct01.csv"))

df_upregulated_merge <- inner_join(df_covid_control, df_ipf_control, by = "gene")
write.csv(df_upregulated_merge, here("EMT/csv/at2_downregulated_merge_minpct01.csv"))