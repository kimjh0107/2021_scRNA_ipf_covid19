library(Seurat)
library(tidyverse)
library(here)
library(presto)
library(msigdbr)
library(fgsea)
library(ggpubr)


#### 1.  HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION ####
mesenchymal <- readRDS(here("EMT/RDS/snRNA_EMT.RDS"))

mesenchymal.genes <- wilcoxauc(mesenchymal, 'Diagnosis')

# Load MsigDB
m_df<- msigdbr(species = "Homo sapiens", category = "H")
m_df <- subset(m_df, subset = gs_name == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

# arrange geneset 
mesenchymal.genes %>% arrange(desc(logFC), desc(auc)) 

# select only the feature and auc columns for fgsea, which statistics to use is an open question
control_genes<- mesenchymal.genes %>%
  dplyr::filter(group == "Control") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)

covid_genes<- mesenchymal.genes %>%
  dplyr::filter(group == "COVID-19") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)


# rank gene sets 
ranks_control <- deframe(control_genes)
ranks_covid <- deframe(covid_genes)

# plot by each diagnosis 
# Control
fgseaRes <- fgsea(fgsea_sets, stats = ranks_control, nperm = 1000) %>%
  as_tibble() %>%
  arrange(desc(NES)) %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  head()

p1 <- plotEnrichment(fgsea_sets[["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"]],
               ranks_control) + labs(title="Control HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION") + theme_classic()

# COVID
fgseaRes <- fgsea(fgsea_sets, stats = ranks_covid, nperm = 1000) %>%
  as_tibble() %>%
  arrange(desc(NES)) %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  head()

p2 <- plotEnrichment(fgsea_sets[["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"]],
               ranks_covid) + labs(title="COVID-19 HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION") + theme_classic()

ggarrange(p1,p2,ncol = 2, nrow = 1)
ggsave(here("EMT/figure/snRNA_GSEA_bulk_EMT_fibro.pdf"), width = 20, height = 5)





#### 2, HALLMARK_TGF_BETA_SIGNALING ####
mesenchymal <- readRDS(here("EMT/RDS/snRNA_EMT.RDS"))


mesenchymal.genes <- wilcoxauc(mesenchymal, 'Diagnosis')

# Load MsigDB
m_df<- msigdbr(species = "Homo sapiens", category = "H")
m_df <- subset(m_df, subset = gs_name == "HALLMARK_TGF_BETA_SIGNALING")
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

# arrange geneset 
mesenchymal.genes %>% arrange(desc(logFC), desc(auc)) 

# select only the feature and auc columns for fgsea, which statistics to use is an open question
control_genes<- mesenchymal.genes %>%
  dplyr::filter(group == "Control") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)

covid_genes<- mesenchymal.genes %>%
  dplyr::filter(group == "COVID-19") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)


# rank gene sets 
ranks_control <- deframe(control_genes)
ranks_covid <- deframe(covid_genes)

# plot by each diagnosis 
# Control
fgseaRes <- fgsea(fgsea_sets, stats = ranks_control, nperm = 1000) %>%
  as_tibble() %>%
  arrange(desc(NES)) %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  head()

p1 <- plotEnrichment(fgsea_sets[["HALLMARK_TGF_BETA_SIGNALING"]],
               ranks_control) + labs(title="Control HALLMARK_TGF_BETA_SIGNALING") + theme_classic()

# COVID
fgseaRes <- fgsea(fgsea_sets, stats = ranks_covid, nperm = 1000) %>%
  as_tibble() %>%
  arrange(desc(NES)) %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  head()

p2 <- plotEnrichment(fgsea_sets[["HALLMARK_TGF_BETA_SIGNALING"]],
               ranks_covid) + labs(title="COVID-19 HALLMARK_TGF_BETA_SIGNALING") + theme_classic()

ggarrange(p1,p2, ncol = 2, nrow = 1)
ggsave(here("EMT/figure/snRNA_GSEA_bulk_TGFB_fibro.pdf"), width = 20, height = 5)






#### 3, HALLMARK_NOTCH_SIGNALING ####
mesenchymal <- readRDS(here("EMT/RDS/snRNA_EMT.RDS"))


mesenchymal.genes <- wilcoxauc(mesenchymal, 'Diagnosis')

# Load MsigDB
m_df<- msigdbr(species = "Homo sapiens", category = "H")
m_df <- subset(m_df, subset = gs_name == "HALLMARK_NOTCH_SIGNALING")
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

# arrange geneset 
mesenchymal.genes %>% arrange(desc(logFC), desc(auc)) 

# select only the feature and auc columns for fgsea, which statistics to use is an open question
control_genes<- mesenchymal.genes %>%
  dplyr::filter(group == "Control") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)

covid_genes<- mesenchymal.genes %>%
  dplyr::filter(group == "COVID-19") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)


# rank gene sets 
ranks_control <- deframe(control_genes)
ranks_covid <- deframe(covid_genes)

# plot by each diagnosis 
# Control
fgseaRes <- fgsea(fgsea_sets, stats = ranks_control, nperm = 1000) %>%
  as_tibble() %>%
  arrange(desc(NES)) %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  head()

p1 <- plotEnrichment(fgsea_sets[["HALLMARK_NOTCH_SIGNALING"]],
               ranks_control) + labs(title="Control HALLMARK_NOTCH_SIGNALING") + theme_classic()

# COVID
fgseaRes <- fgsea(fgsea_sets, stats = ranks_covid, nperm = 1000) %>%
  as_tibble() %>%
  arrange(desc(NES)) %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  head()

p2 <- plotEnrichment(fgsea_sets[["HALLMARK_NOTCH_SIGNALING"]],
               ranks_covid) + labs(title="COVID-19 HALLMARK_NOTCH_SIGNALING") + theme_classic()

ggarrange(p1,p2, ncol = 2, nrow = 1)
ggsave(here("EMT/figure/snRNA_GSEA_bulk_NOTCH_fibro.pdf"), width = 20, height = 5)





#### 4, HALLMARK_TNFA_SIGNALING_VIA_NFKB ####
mesenchymal <- readRDS(here("EMT/RDS/snRNA_EMT.RDS"))


mesenchymal.genes <- wilcoxauc(mesenchymal, 'Diagnosis')

# Load MsigDB
m_df<- msigdbr(species = "Homo sapiens", category = "H")
m_df <- subset(m_df, subset = gs_name == "HALLMARK_TNFA_SIGNALING_VIA_NFKB")
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

# arrange geneset 
mesenchymal.genes %>% arrange(desc(logFC), desc(auc)) 

# select only the feature and auc columns for fgsea, which statistics to use is an open question
control_genes<- mesenchymal.genes %>%
  dplyr::filter(group == "Control") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)

covid_genes<- mesenchymal.genes %>%
  dplyr::filter(group == "COVID-19") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)


# rank gene sets 
ranks_control <- deframe(control_genes)
ranks_covid <- deframe(covid_genes)

# plot by each diagnosis 
# Control
fgseaRes <- fgsea(fgsea_sets, stats = ranks_control, nperm = 1000) %>%
  as_tibble() %>%
  arrange(desc(NES)) %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  head()

p1 <- plotEnrichment(fgsea_sets[["HALLMARK_TNFA_SIGNALING_VIA_NFKB"]],
               ranks_control) + labs(title="Control HALLMARK_TNFA_SIGNALING_VIA_NFKB") + theme_classic()

# COVID
fgseaRes <- fgsea(fgsea_sets, stats = ranks_covid, nperm = 1000) %>%
  as_tibble() %>%
  arrange(desc(NES)) %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  head()

p2 <- plotEnrichment(fgsea_sets[["HALLMARK_TNFA_SIGNALING_VIA_NFKB"]],
               ranks_covid) + labs(title="COVID-19 HALLMARK_TNFA_SIGNALING_VIA_NFKB") + theme_classic()

ggarrange(p1,p2, ncol = 2, nrow = 1)
ggsave(here("EMT/figure/snRNA_GSEA_bulk_TNFA_fibro.pdf"), width = 20, height = 5)






#### 5, HALLMARK_WNT_BETA_CATENIN_SIGNALING ####
mesenchymal <- readRDS(here("EMT/RDS/snRNA_EMT.RDS"))


mesenchymal.genes <- wilcoxauc(mesenchymal, 'Diagnosis')

# Load MsigDB
m_df<- msigdbr(species = "Homo sapiens", category = "H")
m_df <- subset(m_df, subset = gs_name == "HALLMARK_WNT_BETA_CATENIN_SIGNALING")
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

# arrange geneset 
mesenchymal.genes %>% arrange(desc(logFC), desc(auc)) 

# select only the feature and auc columns for fgsea, which statistics to use is an open question
control_genes<- mesenchymal.genes %>%
  dplyr::filter(group == "Control") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)

covid_genes<- mesenchymal.genes %>%
  dplyr::filter(group == "COVID-19") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)


# rank gene sets 
ranks_control <- deframe(control_genes)
ranks_covid <- deframe(covid_genes)

# plot by each diagnosis 
# Control
fgseaRes <- fgsea(fgsea_sets, stats = ranks_control, nperm = 1000) %>%
  as_tibble() %>%
  arrange(desc(NES)) %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  head()

p1 <- plotEnrichment(fgsea_sets[["HALLMARK_WNT_BETA_CATENIN_SIGNALING"]],
               ranks_control) + labs(title="Control HALLMARK_WNT_BETA_CATENIN_SIGNALING") + theme_classic()

# COVID
fgseaRes <- fgsea(fgsea_sets, stats = ranks_covid, nperm = 1000) %>%
  as_tibble() %>%
  arrange(desc(NES)) %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  head()

p2 <- plotEnrichment(fgsea_sets[["HALLMARK_WNT_BETA_CATENIN_SIGNALING"]],
               ranks_covid) + labs(title="COVID-19 HALLMARK_WNT_BETA_CATENIN_SIGNALING") + theme_classic()

ggarrange(p1,p2, ncol = 2, nrow = 1)
ggsave(here("EMT/figure/snRNA_GSEA_bulk_WNT_fibro.pdf"), width = 20, height = 5)
