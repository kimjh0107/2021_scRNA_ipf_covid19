library(Seurat)
library(tidyverse)
library(here)
library(presto)
library(msigdbr)
library(fgsea)
library(ggpubr)


#### 1.  HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION ####
mesenchymal <- readRDS(here("EMT/RDS/fibro.RDS"))

mesenchymal.genes <- wilcoxauc(mesenchymal, 'Diagnosis')

# Load MsigDB
m_df<- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP")
m_df@gs_name
m_df <- subset(m_df, subset = gs_name == "REACTOME_SIGNALING_BY_FGFR_IN_DISEASE")
m_df
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

# arrange geneset 
mesenchymal.genes %>% arrange(desc(logFC), desc(auc)) 

# select only the feature and auc columns for fgsea, which statistics to use is an open question
control_genes<- mesenchymal.genes %>%
  dplyr::filter(group == "Control") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)

covid_genes<- mesenchymal.genes %>%
  dplyr::filter(group == "COVID") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)

ipf_genes<- mesenchymal.genes %>%
  dplyr::filter(group == "IPF") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)

# rank gene sets 
ranks_control <- deframe(control_genes)
ranks_covid <- deframe(covid_genes)
ranks_ipf <- deframe(ipf_genes)

# plot by each diagnosis 
# Control
fgseaRes <- fgsea(fgsea_sets, stats = ranks_control, nperm = 1000) %>%
  as_tibble() %>%
  arrange(desc(NES)) %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  head()

p1 <- plotEnrichment(fgsea_sets[["REACTOME_SIGNALING_BY_FGFR_IN_DISEASE"]],
               ranks_control) + labs(title="Control REACTOME_SIGNALING_BY_FGFR_IN_DISEASE") + theme_classic()

# COVID
fgseaRes <- fgsea(fgsea_sets, stats = ranks_covid, nperm = 1000) %>%
  as_tibble() %>%
  arrange(desc(NES)) %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  head()

p2 <- plotEnrichment(fgsea_sets[["REACTOME_SIGNALING_BY_FGFR_IN_DISEASE"]],
               ranks_covid) + labs(title="COVID-19 REACTOME_SIGNALING_BY_FGFR_IN_DISEASE") + theme_classic()

# COVID
fgseaRes <- fgsea(fgsea_sets, stats = ranks_ipf, nperm = 1000) %>%
  as_tibble() %>%
  arrange(desc(NES)) %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  head()

p3 <- plotEnrichment(fgsea_sets[["REACTOME_SIGNALING_BY_FGFR_IN_DISEASE"]],
               ranks_ipf) + labs(title="IPF REACTOME_SIGNALING_BY_FGFR_IN_DISEASE") + theme_classic()

ggarrange(p1,p2,p3, ncol = 3, nrow = 1)
ggsave(here("EMT/figure/GSEA_fibro_FGFR_fibro.pdf"), width = 20, height = 5)
