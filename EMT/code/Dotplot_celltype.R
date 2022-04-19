library(here)
library(Seurat)
library(tidyverse)

df <- readRDS(here("EMT/RDS/df.RDS"))
DefaultAssay(df) <- "RNA"
Idents(df) <- "new_submain"
levels(df) <- df$new_submain
genesets <- c(
    "CLIC5", "AGER", "PDPN",
    "SFTPC", "SFTPB", "SFTPD", "ETV5", "PDPN", "MUC1",
    "CD79A", "CD24", "MS4A1", "CD19",
    "KRT5", "KRT14", "TP63", "DAPL1",
    "CD3E", "CD8A", "COTL1", 
    "CD3E", "CD8A", "GZMK", "DUSP2",
    "FOXJ1", "TUBB1", "TP73", "CCDC78",
    "SCGB3A2", "CCKAR",
    "PECAM1",
    "COL1A1", "PDGFRA",
    "MUC5B", "MUC5AC", "SPDEF",
    "MARCO", "MSR1", "MRC1",
     "CD1C", "PLD4",
    "CD14", "S100A8", "	FCGR3A", 
    "S100A8", "S100A9", "IFITM2", "FCGR3B",
    "KLRD1", "NKG7", "TYROBP", "NCAM1",
    "LILRB4", "IRF8", "LILRA4",
    "CD79A", "CD27", "SLAMF7",
    "CNN1", "ACTA2", "TAGLN", "RGS5"
)

DotPlot(df, features = genesets)






#### cell 논문 signature geneset list #### 

library(here)
library(Seurat)
library(tidyverse)

fibro <- readRDS(here("EMT/RDS/fibro.RDS"))

fibro
DefaultAssay(fibro) <- "RNA"
Idents(fibro) <- "Diagnosis"

fibro_genesets <- c(
    "COL1A1", "COL1A2", "COL3A1", "COL4A1", "COL4A2", "COL5A1", "COL6A1",
    "COL6A2", "COL6A3", "PCOLCE", "FN1", "LAMB1", "LAMP1", "LUM", "MMP2", "MMP14",
    "SPARC", "PLIN2", "FBLIM1", "PTX3"
)

myofibro_genesets <- c(
    "COL1A1", "COL1A2", "COL3A1", "COL4A1", "COL4A2", "COL5A1", "COL6A1",
    "COL6A3", "COL18A1", "COL27A1", "COL4A3BP", "PCOLCE", "FN1", "FBLIM1",
    "LAMA2", "LAMB1", "LAMC1", "LUM", "MMP2", "MMP14", "KIAA1199", "PLIN2",
    "MYLK", "MYO10", "MTO16", "MOXD1", "CCBE1"
)

DotPlot(fibro, features = fibro_genesets, cols = c("lightgrey", "red", "purple"))
ggsave(here("EMT/figure/Dotplot.pdf"))

