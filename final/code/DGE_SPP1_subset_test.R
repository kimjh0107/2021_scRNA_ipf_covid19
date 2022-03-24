library(Seurat)
library(tidyverse)
library(future)
library(here)

Macrophages <- readRDS(here("final_before/RDS/macro_annotation.RDS"))

unique(Macrophages$new_submain)

Macrophages <- subset(Macrophages, subset = new_submain == "SPP1+")

 Macrophages
Idents(Macrophages) <- "Diagnosis"
DefaultAssay(Macrophages) <- "RNA"

# Upregulated
Macrophages_covid_control <- FindMarkers(Macrophages, ident.1 = "COVID-19", ident.2 = "Control", only.pos = TRUE, min.pct = 0.25)
Macrophages_covid_control$gene <- rownames(Macrophages_covid_control)
write.csv(Macrophages_covid_control, here("final/csv/Macrophages_SPP1_test_upregulated_covid.csv"))

Macrophages_ipf_control <- FindMarkers(Macrophages, ident.1 = "IPF", ident.2 = "Control", only.pos = TRUE, min.pct = 0.25)
Macrophages_ipf_control$gene <- rownames(Macrophages_ipf_control)
write.csv(Macrophages_ipf_control, here("final/csv/Macrophages_SPP1_test_upregulated_ipf.csv"))

Macrophages_upregulated_merge <- inner_join(Macrophages_covid_control, Macrophages_ipf_control, by = "gene")
write.csv(Macrophages_upregulated_merge, here("final/csv/Macrophages_SPP1_test_upregulated_merge.csv"))


# Downregulated
Macrophages_covid_control <- FindMarkers(Macrophages, ident.1 = "Control", ident.2 = "COVID-19", only.pos = TRUE, min.pct = 0.25)
Macrophages_covid_control$gene <- rownames(Macrophages_covid_control)
write.csv(Macrophages_covid_control, here("final/csv/Macrophages_SPP1_test_downregulated_covid.csv"))

Macrophages_ipf_control <- FindMarkers(Macrophages, ident.1 = "Control", ident.2 = "IPF", only.pos = TRUE, min.pct = 0.25)
Macrophages_ipf_control$gene <- rownames(Macrophages_ipf_control)
write.csv(Macrophages_ipf_control, here("final/csv/Macrophages_SPP1_test_downregulated_ipf.csv"))

Macrophages_upregulated_merge <- inner_join(Macrophages_covid_control, Macrophages_ipf_control, by = "gene")
write.csv(Macrophages_upregulated_merge, here("final/csv/Macrophages_SPP1_test_downregulated_merge.csv"))





Macrophages <- readRDS(here("final_before/RDS/macro_annotation.RDS"))
Macrophages <- subset(Macrophages, subset = new_submain == "SPP1+FCN1+")

 Macrophages
Idents(Macrophages) <- "Diagnosis"
DefaultAssay(Macrophages) <- "RNA"

# Upregulated
Macrophages_covid_control <- FindMarkers(Macrophages, ident.1 = "COVID-19", ident.2 = "Control", only.pos = TRUE, min.pct = 0.25)
Macrophages_covid_control$gene <- rownames(Macrophages_covid_control)
write.csv(Macrophages_covid_control, here("final/csv/Macrophages_SPP1_FCN1_test_upregulated_covid.csv"))

Macrophages_ipf_control <- FindMarkers(Macrophages, ident.1 = "IPF", ident.2 = "Control", only.pos = TRUE, min.pct = 0.25)
Macrophages_ipf_control$gene <- rownames(Macrophages_ipf_control)
write.csv(Macrophages_ipf_control, here("final/csv/Macrophages_SPP1_FCN1_test_upregulated_ipf.csv"))

Macrophages_upregulated_merge <- inner_join(Macrophages_covid_control, Macrophages_ipf_control, by = "gene")
write.csv(Macrophages_upregulated_merge, here("final/csv/Macrophages_SPP1_FCN1_test_upregulated_merge.csv"))


# Downregulated
Macrophages_covid_control <- FindMarkers(Macrophages, ident.1 = "Control", ident.2 = "COVID-19", only.pos = TRUE, min.pct = 0.25)
Macrophages_covid_control$gene <- rownames(Macrophages_covid_control)
write.csv(Macrophages_covid_control, here("final/csv/Macrophages_SPP1_FCN1_test_downregulated_covid.csv"))

Macrophages_ipf_control <- FindMarkers(Macrophages, ident.1 = "Control", ident.2 = "IPF", only.pos = TRUE, min.pct = 0.25)
Macrophages_ipf_control$gene <- rownames(Macrophages_ipf_control)
write.csv(Macrophages_ipf_control, here("final/csv/Macrophages_SPP1_FCN1_test_downregulated_ipf.csv"))

Macrophages_upregulated_merge <- inner_join(Macrophages_covid_control, Macrophages_ipf_control, by = "gene")
write.csv(Macrophages_upregulated_merge, here("final/csv/Macrophages_SPP1_FCN1_test_downregulated_merge.csv"))