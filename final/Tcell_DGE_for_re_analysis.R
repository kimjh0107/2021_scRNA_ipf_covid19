library(Seurat)
library(tidyverse)
library(here)

tcell <- readRDS(here("final/RDS/T_NK.RDS"))



# tcell
Idents(tcell) <- "Diagnosis"
DefaultAssay(tcell) <- "RNA"

# Upregulated
tcell_covid_control <- FindMarkers(tcell, ident.1 = "COVID-19", ident.2 = "Control", only.pos = TRUE, min.pct = 0.25, test.use = "MAST")
tcell_covid_control$gene <- rownames(tcell_covid_control)
write.csv(tcell_covid_control, here("final/csv/tcell_upregulated_covid.csv"))

tcell_ipf_control <- FindMarkers(tcell, ident.1 = "IPF", ident.2 = "Control", only.pos = TRUE, min.pct = 0.25, test.use = "MAST")
tcell_ipf_control$gene <- rownames(tcell_ipf_control)
write.csv(tcell_ipf_control, here("final/csv/tcell_upregulated_ipf.csv"))

tcell_upregulated_merge <- inner_join(tcell_covid_control, tcell_ipf_control, by = "gene")
write.csv(tcell_upregulated_merge, here("final/csv/tcell_upregulated_merge.csv"))

tcell_upregulated_merge <- tcell_upregulated_merge %>% filter(avg_log2FC > 0.5)
write_csv(covid_up_filter,here("final/csv/tcell_upregulated_merge_filter.csv"))


# Downregulated_min0
tcell_covid_control <- FindMarkers(tcell, ident.1 = "Control", ident.2 = "COVID-19", only.pos = TRUE, min.pct = 0.25, test.use = "MAST")
tcell_covid_control$gene <- rownames(tcell_covid_control)
write.csv(tcell_covid_control, here("final/csv/tcell_downregulated_covid.csv"))

tcell_ipf_control <- FindMarkers(tcell, ident.1 = "Control", ident.2 = "IPF", only.pos = TRUE, min.pct = 0.25, test.use = "MAST")
tcell_ipf_control$gene <- rownames(tcell_ipf_control)
write.csv(tcell_ipf_control, here("final/csv/tcell_downregulated_ipf.csv"))

tcell_upregulated_merge <- inner_join(tcell_covid_control, tcell_ipf_control, by = "gene")
write.csv(tcell_upregulated_merge, here("final/csv/tcell_downregulated_merge.csv"))

tcell_upregulated_merge <- tcell_upregulated_merge %>% filter(avg_log2FC > 0.5)
write_csv(tcell_upregulated_merge, here("final/csv/tcell_downregulated_merge_filter.csv"))
