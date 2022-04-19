library(here)
library(Seurat)
library(tidyverse)
library(CellChat)
library(ggpubr)
#### Lung fibrosis genesets #### 
emt <- readRDS(here("EMT/RDS/EMT.RDS"))

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
ggsave(here("EMT/figure/dotplot_fibro_by_lungfibro_gene.pdf"), width = 12, height = 4)


#### at2
DefaultAssay(at2) <- "RNA"
Idents(at2) <- "Diagnosis"
fibro_geneset <- read_csv(here("EMT/csv/lung_fibrosis_genesets.txt"))
at2_geneset

DotPlot(at2, features = fibro_geneset$genes, cols = c("lightgrey", "red"), col.min = -2.5, col.max = 2.5,
dot.min =0, dot.scale= 6) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
ggsave(here("EMT/figure/dotplot_at2_by_lungfibro_gene.pdf"), width = 12, height = 4)

#### at1
DefaultAssay(at1) <- "RNA"
Idents(at1) <- "Diagnosis"
fibro_geneset <- read_csv(here("EMT/csv/lung_fibrosis_genesets.txt"))
at1_geneset

DotPlot(at1, features = fibro_geneset$genes, cols = c("lightgrey", "red"), col.min = -2.5, col.max = 2.5,
dot.min =0, dot.scale= 6) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
ggsave(here("EMT/figure/dotplot_at1_by_lungfibro_gene.pdf"), width = 12, height = 4)




#### EMT genesets #### 
emt_genesets <- read_csv(here("EMT/csv/EMT_genesets_GO.txt"))
emt_genesets
emt_genesets <- emt_genesets %>% filter(! duplicated(genes))
change_gene <- toupper(emt_genesets$genes)  # 대문자로 변경 
emt_genesets$genes <- change_gene
emt_genesets

#### fibro부터 우선적으로 확인 
DefaultAssay(fibro) <- "RNA"
Idents(fibro) <- "Diagnosis"

DotPlot(fibro, features = emt_genesets$genes, cols = c("lightgrey", "red"), col.min = -2.5, col.max = 2.5,
dot.min =0, dot.scale= 6) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6))
ggsave(here("EMT/figure/dotplot_fibro_by_EMT_gene.pdf"), width = 18, height = 4)

#### AT2  
DefaultAssay(at2) <- "RNA"
Idents(at2) <- "Diagnosis"

DotPlot(at2, features = emt_genesets$genes, cols = c("lightgrey", "red"), col.min = -2.5, col.max = 2.5,
dot.min =0, dot.scale= 6) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6))
ggsave(here("EMT/figure/dotplot_at2_by_EMT_gene.pdf"), width = 18, height = 4)

#### AT1
DefaultAssay(at1) <- "RNA"
Idents(at1) <- "Diagnosis"

DotPlot(at1, features = emt_genesets$genes, cols = c("lightgrey", "red"), col.min = -2.5, col.max = 2.5,
dot.min =0, dot.scale= 6) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6))
ggsave(here("EMT/figure/dotplot_at1_by_EMT_gene.pdf"), width = 18, height = 4)






#### MET genesets #### 
met_genesets <- read_csv(here("EMT/csv/MET_genesets.txt"))
met_genesets

#### fibro부터 우선적으로 확인 
DefaultAssay(fibro) <- "RNA"
Idents(fibro) <- "Diagnosis"

DotPlot(fibro, features = met_genesets$genes, cols = c("lightgrey", "red"), col.min = -2.5, col.max = 2.5,
dot.min =0, dot.scale= 6) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6))
ggsave(here("EMT/figure/dotplot_fibro_by_MET_gene.pdf"), width = 6, height = 4)

#### AT2  
DefaultAssay(at2) <- "RNA"
Idents(at2) <- "Diagnosis"

DotPlot(at2, features = met_genesets$genes, cols = c("lightgrey", "red"), col.min = -2.5, col.max = 2.5,
dot.min =0, dot.scale= 6) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6))
ggsave(here("EMT/figure/dotplot_at2_by_MET_gene.pdf"), width = 6, height = 4)

#### AT1
DefaultAssay(at1) <- "RNA"
Idents(at1) <- "Diagnosis"

DotPlot(at1, features = met_genesets$genes, cols = c("lightgrey", "red"), col.min = -2.5, col.max = 2.5,
dot.min =0, dot.scale= 6) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6))
ggsave(here("EMT/figure/dotplot_at1_by_MET_gene.pdf"), width = 6, height = 4)




p1 <- VlnPlot(fibro, features = c("HGF"))
p2 <- VlnPlot(fibro, features = c("TNC"))
p3 <- VlnPlot(fibro, features = c("FGF7"))
p4 <- VlnPlot(fibro, features = c("FGFR1"))

p5 <- VlnPlot(fibro, features = c("WNT5A"))
p6 <- VlnPlot(fibro, features = c("FGF10"))
p7 <- VlnPlot(fibro, features = c("WNT2B"))
p8 <- VlnPlot(fibro, features = c("HOXA5"))
p9 <- VlnPlot(fibro, features = c("MET"))
p10 <- VlnPlot(fibro, features = c("KLF4"))

plot_arrange <- ggarrange(p1,p2,p3,p4,p5,
                          p6,p7,p8,p9,p10,
                          ncol =5, nrow = 2)
ggsave(here("EMT/figure/vlnplot_MET.pdf"), width = 20, height = 10)
plot_arrange




#### Barplot #### 
# function for emt scoring 
get_gene_emt_score <- function(gene){
  mean_exp <- exp_data %>% dplyr::select({{gene}}) %>% pull() %>% mean()
  return (tibble(gene = gene, emt_score = mean_exp))
}

final_score <- function(exp_data){
  final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
  final_score <- final$emt_score %>% mean()
  return (final_score)
}

get_plot <- function(celltype){
  merge_data_frame %>% 
    pivot_longer(cols = (-Diagnosis), names_to="cell", values_to="exp") %>%
    filter(cell == celltype) %>% 
    ggplot(aes(x = Diagnosis, y = exp)) + geom_bar(stat = 'identity', fill = c("skyblue", "red", "purple"), alpha= 0.6 ) + theme_classic() +
      labs(title= celltype) 
}

fibro <- readRDS(here("final/RDS/EMT.RDS"))
fibro <- subset(fibro, subset = new_submain == "Fibroblasts")
control <- subset(fibro, subset = Diagnosis == "Control")
covid <- subset(fibro, subset = Diagnosis == "COVID")
ipf <- subset(fibro, subset = Diagnosis == "IPF")

DefaultAssay(control) <- "RNA"
DefaultAssay(covid) <- "RNA"
DefaultAssay(ipf) <- "RNA"

met_genesets <- read_csv(here("EMT/csv/MET_genesets.txt"))
gene_list <- met_genesets$genes
gene_list

# control
exp_data <- FetchData(object = control, vars = gene_list)
MET_control_score <- final_score(exp_data)
# covid 
exp_data <- FetchData(object = covid, vars = gene_list)
MET_covid_score <- final_score(exp_data)
# ipf 
exp_data <- FetchData(object = ipf, vars = gene_list)
MET_ipf_score <- final_score(exp_data)

MET_score <- c(MET_control_score, MET_covid_score, MET_ipf_score)
MET_score
Diagnosis <- c("Control", "COVID-19", "IPF")
merge_data_frame <- data.frame(Diagnosis, MET_score)
merge_data_frame

figure <- get_plot("MET_score")
figure
ggsave(plot = figure, here("EMT/figure//Barplot_MET_scoring.pdf"), height = 4, width = 3)


met_genesets <- read_csv(here("EMT/csv/MET_genesets2.txt"))
gene_list <- met_genesets$genes
gene_list

# control
exp_data <- FetchData(object = control, vars = gene_list)
MET_control_score <- final_score(exp_data)
# covid 
exp_data <- FetchData(object = covid, vars = gene_list)
MET_covid_score <- final_score(exp_data)
# ipf 
exp_data <- FetchData(object = ipf, vars = gene_list)
MET_ipf_score <- final_score(exp_data)

MET_score <- c(MET_control_score, MET_covid_score, MET_ipf_score)
MET_score
Diagnosis <- c("Control", "COVID-19", "IPF")
merge_data_frame <- data.frame(Diagnosis, MET_score)
merge_data_frame

figure <- get_plot("MET_score")
figure
ggsave(plot = figure, here("EMT/figure//Barplot_MET_scoring_by_only_expressed_in_covid.pdf"), height = 4, width = 3)
