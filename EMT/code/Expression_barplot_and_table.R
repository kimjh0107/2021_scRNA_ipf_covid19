library(Seurat)
library(tidyverse)
library(here)


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


at <- readRDS(here("EMT/RDS/AT.RDS"))
at1 <- subset(at, subset = new_submain %in% c("AT1"))
at2 <- subset(at, subset = new_submain %in% c("AT2"))

incomple_transition_marker_gene <- c("ETV5")
late_AT1_maturation <- c("CAV1")


#### incomple_transition_marker_gene #### 
# at 전체 데이터셋에 대해서 분석 진행 
control <- subset(at2, subset = Diagnosis == "Control")
covid <- subset(at2, subset = Diagnosis == "COVID")
ipf <- subset(at2, subset = Diagnosis == "IPF")
DefaultAssay(control) <- "RNA"
DefaultAssay(covid) <- "RNA"
DefaultAssay(ipf) <- "RNA"

# Control 
exp_data <- FetchData(object = control, vars = incomple_transition_marker_gene)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
control_final_emt_score <- final$emt_score %>% mean()

# COVID-19 
exp_data <- FetchData(object = covid, vars = incomple_transition_marker_gene)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
covid_final_emt_score <- final$emt_score %>% mean()

# IPF
exp_data <- FetchData(object = ipf, vars = incomple_transition_marker_gene)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
ipf_final_emt_score <- final$emt_score %>% mean()

Incomple_transition_marker_expression <- c(control_final_emt_score, covid_final_emt_score, ipf_final_emt_score)


control <- subset(at1, subset = Diagnosis == "Control")
covid <- subset(at1, subset = Diagnosis == "COVID")
ipf <- subset(at1, subset = Diagnosis == "IPF")
DefaultAssay(control) <- "RNA"
DefaultAssay(covid) <- "RNA"
DefaultAssay(ipf) <- "RNA"

# Control 
exp_data <- FetchData(object = control, vars = late_AT1_maturation)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
control_final_emt_score <- final$emt_score %>% mean()

# COVID-19 
exp_data <- FetchData(object = covid, vars = late_AT1_maturation)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
covid_final_emt_score <- final$emt_score %>% mean()

# IPF
exp_data <- FetchData(object = ipf, vars = late_AT1_maturation)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
ipf_final_emt_score <- final$emt_score %>% mean()

late_AT1_maturation_marker_expression <- c(control_final_emt_score, covid_final_emt_score, ipf_final_emt_score)


Diagnosis <- c("Control", "COVID", "IPF")
merge_data_frame <- data.frame(Diagnosis,
                               Incomple_transition_marker_expression,
                               late_AT1_maturation_marker_expression
)
merge_data_frame    
write_csv(merge_data_frame, here("EMT/csv/AT_late_maturation_and_incomplete_transition_expression_table.csv"))


figure <- get_plot("late_AT1_maturation_marker_expression")
ggsave(plot = figure, here("EMT/figure/Figure3/pre_at1_barplot_expression.pdf"), height = 4, width = 3)

figure <- get_plot("Incomple_transition_marker_expression")
ggsave(plot = figure, here("EMT/figure/Figure3/incomplete_transition_barplot_expression.pdf"), height = 4, width = 3)






#### pathogenic fibroblasts gene sets  ####
fibro <- readRDS(here("EMT/RDS/fibro.RDS"))
control <- subset(fibro, subset = Diagnosis == "Control")
covid <- subset(fibro, subset = Diagnosis == "COVID")
ipf <- subset(fibro, subset = Diagnosis == "IPF")
DefaultAssay(control) <- "RNA"
DefaultAssay(covid) <- "RNA"
DefaultAssay(ipf) <- "RNA"

pathogenic_gene <- c("CTHRC1", "COL1A1", "COL3A1")


#### 4 ####
# Control 
exp_data <- FetchData(object = control, vars = pathogenic_gene)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
control_final_emt_score <- final$emt_score %>% mean()

# COVID-19 
exp_data <- FetchData(object = covid, vars = pathogenic_gene)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
covid_final_emt_score <- final$emt_score %>% mean()

# IPF
exp_data <- FetchData(object = ipf, vars = pathogenic_gene)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
ipf_final_emt_score <- final$emt_score %>% mean()

Pathogenic_geneic_fibrosis_expression <- c(control_final_emt_score, covid_final_emt_score, ipf_final_emt_score)
Diagnosis <- c("Control", "COVID-19", "IPF")
merge_data_frame <- data.frame(Diagnosis, 
                               Pathogenic_geneic_fibrosis_expression
)
write_csv(merge_data_frame, here("EMT/csv/pathogenic_fibrosis_scRNA_expression_table.csv"))

figure <- get_plot("Pathogenic_geneic_fibrosis_expression")
ggsave(plot = figure, here("EMT/figure/Figure3/pathogenic_fibrosis_barplot.pdf"), height = 4, width = 3)


# ECM composition check 
reactome_2 <- read_csv(here("final/RDS/ECM_geneset/REACTOME_ECM_PROTEOGLYCANS.csv"))
reactome_2_list <- reactome_2$Gene.Symbol

#### Reactome 2 #### 
# Control 
exp_data <- FetchData(object = control, vars = reactome_2_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
control_final_emt_score <- final$emt_score %>% mean()

# COVID-19 
exp_data <- FetchData(object = covid, vars = reactome_2_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
covid_final_emt_score <- final$emt_score %>% mean()

# IPF
exp_data <- FetchData(object = ipf, vars = reactome_2_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
ipf_final_emt_score <- final$emt_score %>% mean()

ECM_score <- c(control_final_emt_score, covid_final_emt_score, ipf_final_emt_score)
Diagnosis <- c("Control", "COVID-19", "IPF")
merge_data_frame <- data.frame(
                                Diagnosis, ECM_score)
write_csv(merge_data_frame, here("EMT/csv/ECM_expression_table.csv"))

figure <- get_plot("ECM_score")
ggsave(plot = figure, here("EMT/figure/Figure3/ECM_barplot.pdf"), height = 4, width = 3)







#### MET  ####
fibro <- readRDS(here("EMT/RDS/fibro.RDS"))
control <- subset(fibro, subset = Diagnosis == "Control")
covid <- subset(fibro, subset = Diagnosis == "COVID")
ipf <- subset(fibro, subset = Diagnosis == "IPF")
DefaultAssay(control) <- "RNA"
DefaultAssay(covid) <- "RNA"
DefaultAssay(ipf) <- "RNA"

# MET markers check out 
MET_markers <- c("TNC",
                 "FGFR1",
                 "WNT5A",
                 "FGF7",
                 "HGF",
                 "FGF10",
                 "WNT2B",
                 "HOXA5")

# Control 
exp_data <- FetchData(object = control, vars = MET_markers)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
control_final_emt_score <- final$emt_score %>% mean()
# COVID-19 
exp_data <- FetchData(object = covid, vars = MET_markers)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
covid_final_emt_score <- final$emt_score %>% mean()

# IPF
exp_data <- FetchData(object = ipf, vars = MET_markers)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
ipf_final_emt_score <- final$emt_score %>% mean()

MET_markers_fibrosis_expression <- c(control_final_emt_score, covid_final_emt_score, ipf_final_emt_score)
Diagnosis <- c("Control", "COVID-19", "IPF")
merge_data_frame <- data.frame(Diagnosis, 
                               MET_markers_fibrosis_expression
)

write_csv(merge_data_frame, here("EMT/csv/MET_scRNA_expression_table.csv"))

figure <- get_plot("MET_expression")
ggsave(plot = figure, here("EMT/figure/Figure3/MET_expression_barplot.pdf"), height = 4, width = 3)












#### FGFR  ####
fibro <- readRDS(here("EMT/RDS/fibro.RDS"))
control <- subset(fibro, subset = Diagnosis == "Control")
covid <- subset(fibro, subset = Diagnosis == "COVID")
ipf <- subset(fibro, subset = Diagnosis == "IPF")
DefaultAssay(control) <- "RNA"
DefaultAssay(covid) <- "RNA"
DefaultAssay(ipf) <- "RNA"

# MET markers check out 
FGFR_markers <- c("BAG4","BCR","CEP43","CNTRL","CPSF6","CUX1","ERLIN2","FGF1","FGF10","FGF16","FGF17","FGF18","FGF2","FGF20","FGF22","FGF23","FGF3","FGF4","FGF5","FGF6","FGF7","FGF8","FGF9","FGFR1","FGFR1OP2","FGFR2","FGFR3","FGFR4","FRS2","GAB1","GAB2","GRB2","GTF2F1","GTF2F2","HRAS","KRAS","LRRFIP1","MYO18A","NCBP1","NCBP2","NRAS","PIK3CA","PIK3R1","PLCG1","POLR2A","POLR2B","POLR2C","POLR2D","POLR2E","POLR2F","POLR2G","POLR2H","POLR2I","POLR2J","POLR2K","POLR2L","SOS1","STAT1","STAT3","STAT5A","STAT5B","TRIM24","ZMYM2")

# Control 
exp_data <- FetchData(object = control, vars = FGFR_markers)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
control_final_emt_score <- final$emt_score %>% mean()
# COVID-19 
exp_data <- FetchData(object = covid, vars = FGFR_markers)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
covid_final_emt_score <- final$emt_score %>% mean()

# IPF
exp_data <- FetchData(object = ipf, vars = FGFR_markers)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
ipf_final_emt_score <- final$emt_score %>% mean()

FGFR_markers_fibrosis_expression <- c(control_final_emt_score, covid_final_emt_score, ipf_final_emt_score)
Diagnosis <- c("Control", "COVID-19", "IPF")
merge_data_frame <- data.frame(Diagnosis, 
                               FGFR_markers_fibrosis_expression
)

write_csv(merge_data_frame, here("EMT/csv/FGFR_markers_scRNA_expression_table.csv"))

figure <- get_plot("FGFR_markers_fibrosis_expression")
ggsave(plot = figure, here("EMT/figure/FGFR_expression_barplot.pdf"), height = 4, width = 3)





# MET markers check out 
FGFR_markers <- c("ANOS1","BRAF","CBL","ESRP1","ESRP2","FGF1","FGF10","FGF16","FGF17","FGF18","FGF19","FGF2","FGF20","FGF22","FGF23","FGF3","FGF4","FGF5","FGF6","FGF7","FGF8","FGF9","FGFBP1","FGFBP2","FGFBP3","FGFR1","FGFR2","FGFR3","FGFR4","FGFRL1","FLRT1","FLRT2","FLRT3","FRS2","FRS3","GAB1","GALNT3","GRB2","GTF2F1","GTF2F2","HNRNPA1","HNRNPF","HNRNPH1","HNRNPM","HRAS","KL","KLB","KRAS","MAPK1","MAPK3","MKNK1","NCBP1","NCBP2","NRAS","PIK3CA","PIK3R1","PLCG1","POLR2A","POLR2B","POLR2C","POLR2D","POLR2E","POLR2F","POLR2G","POLR2H","POLR2I","POLR2J","POLR2K","POLR2L","PPP2CA","PPP2CB","PPP2R1A","PTBP1","PTPN11","RBFOX2","RPS27A","SHC1","SOS1","SPRED1","SPRED2","SPRY2","SRC","TIA1","TIAL1","UBA52","UBB","UBC")

# Control 
exp_data <- FetchData(object = control, vars = FGFR_markers)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
control_final_emt_score <- final$emt_score %>% mean()
# COVID-19 
exp_data <- FetchData(object = covid, vars = FGFR_markers)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
covid_final_emt_score <- final$emt_score %>% mean()

# IPF
exp_data <- FetchData(object = ipf, vars = FGFR_markers)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
ipf_final_emt_score <- final$emt_score %>% mean()

FGFR_markers_fibrosis_expression <- c(control_final_emt_score, covid_final_emt_score, ipf_final_emt_score)
Diagnosis <- c("Control", "COVID-19", "IPF")
merge_data_frame <- data.frame(Diagnosis, 
                               FGFR_markers_fibrosis_expression
)

figure <- get_plot("FGFR_markers_fibrosis_expression")
ggsave(plot = figure, here("EMT/figure/FGFR_expression_barplot_2.pdf"), height = 4, width = 3)
