library(Seurat)
library(tidyverse)
library(here)
library(ggpubr)

tcell <- readRDS(here("final/RDS/T_NK.RDS"))




Idents(tcell) <- "Diagnosis"
DefaultAssay(tcell) <- "RNA"

# Upregulated
tcell_covid_control <- FindMarkers(tcell, ident.1 = "COVID", ident.2 = "Control", only.pos = TRUE, min.pct = 0.25)
tcell_covid_control$gene <- rownames(tcell_covid_control)
write.csv(tcell_covid_control, here("final/csv/tcell_upregulated_covid.csv"))

tcell_ipf_control <- FindMarkers(tcell, ident.1 = "IPF", ident.2 = "Control", only.pos = TRUE, min.pct = 0.25)
tcell_ipf_control$gene <- rownames(tcell_ipf_control)
write.csv(tcell_ipf_control, here("final/csv/tcell_upregulated_ipf.csv"))

tcell_upregulated_merge <- inner_join(tcell_covid_control, tcell_ipf_control, by = "gene")
write.csv(tcell_upregulated_merge, here("final/csv/tcell_upregulated_merge.csv"))

# Downregulated
tcell_covid_control <- FindMarkers(tcell, ident.1 = "Control", ident.2 = "COVID", only.pos = TRUE, min.pct = 0.25)
tcell_covid_control$gene <- rownames(tcell_covid_control)
write.csv(tcell_covid_control, here("final/csv/tcell_downregulated_covid.csv"))

tcell_ipf_control <- FindMarkers(tcell, ident.1 = "Control", ident.2 = "IPF", only.pos = TRUE, min.pct = 0.25)
tcell_ipf_control$gene <- rownames(tcell_ipf_control)
write.csv(tcell_ipf_control, here("final/csv/tcell_downregulated_ipf.csv"))

tcell_upregulated_merge <- inner_join(tcell_covid_control, tcell_ipf_control, by = "gene")
write.csv(tcell_upregulated_merge, here("final/csv/tcell_downregulated_merge.csv"))
