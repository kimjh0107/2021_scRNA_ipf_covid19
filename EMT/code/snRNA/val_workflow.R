library(here)
library(Seurat)
library(tidyverse)

emt <- readRDS(here("EMT/RDS/snRNA_EMT.RDS"))
DefaultAssay(df) <- ranks_control

at1 <- subset(emt, subset = new_submain == "AT1")
at2 <- subset(emt, subset = new_submain == "AT2")
fibro <- subset(emt, subset = new_submain == "Fibroblasts")


#### fibro부터 우선적으로 확인 
DefaultAssay(fibro) <- "RNA"
Idents(fibro) <- "Diagnosis"
fibro_geneset <- read_csv(here("EMT/csv/lung_fibrosis_genesets.txt"))
fibro_geneset

DotPlot(fibro, features = fibro_geneset$genes, cols = c("lightgrey", "red"), col.min = -2.5, col.max = 2.5,
dot.min =0, dot.scale= 6) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
ggsave(here("EMT/figure/snRNA_dotplot_fibro_by_lungfibro_gene.pdf"), width = 12, height = 4)


#### at2
DefaultAssay(at2) <- "RNA"
Idents(at2) <- "Diagnosis"
fibro_geneset <- read_csv(here("EMT/csv/lung_fibrosis_genesets.txt"))
at2_geneset

DotPlot(at2, features = fibro_geneset$genes, cols = c("lightgrey", "red"), col.min = -2.5, col.max = 2.5,
dot.min =0, dot.scale= 6) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
ggsave(here("EMT/figure/snRNA_dotplot_at2_by_lungfibro_gene.pdf"), width = 12, height = 4)

#### at1
DefaultAssay(at1) <- "RNA"
Idents(at1) <- "Diagnosis"
fibro_geneset <- read_csv(here("EMT/csv/lung_fibrosis_genesets.txt"))
at1_geneset

DotPlot(at1, features = fibro_geneset$genes, cols = c("lightgrey", "red"), col.min = -2.5, col.max = 2.5,
dot.min =0, dot.scale= 6) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
ggsave(here("EMT/figure/snRNA_dotplot_at1_by_lungfibro_gene.pdf"), width = 12, height = 4)





#### EMT genesets #### 
emt_genesets <- read_csv(here("EMT/csv/EMT_genesets_GO.txt"))
emt_genesets <- emt_genesets %>% filter(! duplicated(genes))
change_gene <- toupper(emt_genesets$genes)  # 대문자로 변경 
emt_genesets$genes <- change_gene

#### fibro부터 우선적으로 확인 
DefaultAssay(fibro) <- "RNA"
Idents(fibro) <- "Diagnosis"

DotPlot(fibro, features = emt_genesets$genes, cols = c("lightgrey", "red"), col.min = -2.5, col.max = 2.5,
dot.min =0, dot.scale= 6) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6))
ggsave(here("EMT/figure/snRNA_dotplot_fibro_by_EMT_gene.pdf"), width = 18, height = 4)

#### AT2  
DefaultAssay(at2) <- "RNA"
Idents(at2) <- "Diagnosis"

DotPlot(at2, features = emt_genesets$genes, cols = c("lightgrey", "red"), col.min = -2.5, col.max = 2.5,
dot.min =0, dot.scale= 6) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6))
ggsave(here("EMT/figure/snRNA_dotplot_at2_by_EMT_gene.pdf"), width = 18, height = 4)

#### AT1
DefaultAssay(at1) <- "RNA"
Idents(at1) <- "Diagnosis"

DotPlot(at1, features = emt_genesets$genes, cols = c("lightgrey", "red"), col.min = -2.5, col.max = 2.5,
dot.min =0, dot.scale= 6) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6))
ggsave(here("EMT/figure/snRNA_dotplot_at1_by_EMT_gene.pdf"), width = 18, height = 4)




#### MET genesets #### 
met_genesets <- read_csv(here("EMT/csv/MET_genesets.txt"))
met_genesets

#### fibro부터 우선적으로 확인 
DefaultAssay(fibro) <- "RNA"
Idents(fibro) <- "Diagnosis"

DotPlot(fibro, features = met_genesets$genes, cols = c("lightgrey", "red"), col.min = -2.5, col.max = 2.5,
dot.min =0, dot.scale= 6) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6))
ggsave(here("EMT/figure/snRNA_dotplot_fibro_by_MET_gene.pdf"), width = 6, height = 4)

#### AT2  
DefaultAssay(at2) <- "RNA"
Idents(at2) <- "Diagnosis"

DotPlot(at2, features = met_genesets$genes, cols = c("lightgrey", "red"), col.min = -2.5, col.max = 2.5,
dot.min =0, dot.scale= 6) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6))
ggsave(here("EMT/figure/snRNA_dotplot_at2_by_MET_gene.pdf"), width = 6, height = 4)

#### AT1
DefaultAssay(at1) <- "RNA"
Idents(at1) <- "Diagnosis"

DotPlot(at1, features = met_genesets$genes, cols = c("lightgrey", "red"), col.min = -2.5, col.max = 2.5,
dot.min =0, dot.scale= 6) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6))
ggsave(here("EMT/figure/snRNA_dotplot_at1_by_MET_gene.pdf"), width = 6, height = 4)

