library(Seurat)
library(tidyverse)
library(here)
library(ggpubr)

df <- readRDS(here("final/RDS/df.RDS"))


AT2 <- subset(df, subset = new_submain == "AT2")
Monocytes <- subset(df, subset = new_submain == "Monocytes")
Endothelial <- subset(df, subset = new_submain == "Endothelial cells")
Macrophages <- subset(df, subset = new_submain == "Macrophages")
NK <- subset(df, subset = new_submain == "NK cells")
Ciliated <- subset(df, subset = new_submain == "Ciliated cells")
Club <- subset(df, subset = new_submain == "Club cells")
mDCs <- subset(df, subset = new_submain == "mDCs")
Bcell <- subset(df, subset = new_submain == "B cells")
AT1 <- subset(df, subset = new_submain == "AT1")
Fibroblasts <- subset(df, subset = new_submain == "Fibroblasts")
Basal <- subset(df, subset = new_submain == "Basal cells")
Plasma <- subset(df, subset = new_submain == "Plasma cells")
Gobelt <- subset(df, subset = new_submain == "Gobelt cells")
CD4 <- subset(df, subset = new_submain == "CD4 Tcells")
Mast <- subset(df, subset = new_submain == "Mast cells")
Smooth <- subset(df, subset = new_submain == "Smooth Muscle cells")
CD8 <- subset(df, subset = new_submain == "CD8 Tcells")
pDCs <- subset(df, subset = new_submain == "pDCs")
Neutrophils <- subset(df, subset = new_submain == "Neutrophils")



# fetch data 
gene_list <- c("COL1A1")

# function for emt scoring 
get_gene_emt_score <- function(gene){
  mean_exp <- exp_data %>% select({{gene}}) %>% pull() %>% mean()
  return (tibble(gene = gene, emt_score = mean_exp))
}

# AT2 
control <- subset(AT2, subset = Diagnosis == "Control")
covid <- subset(AT2, subset = Diagnosis == "COVID")
ipf <- subset(AT2, subset = Diagnosis == "IPF")
DefaultAssay(control) <- "RNA"
DefaultAssay(covid) <- "RNA"
DefaultAssay(ipf) <- "RNA"

# Control 
exp_data <- FetchData(object = control, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
control_final_emt_score <- final$emt_score %>% mean()

# COVID-19 
exp_data <- FetchData(object = covid, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
covid_final_emt_score <- final$emt_score %>% mean()

# IPF
exp_data <- FetchData(object = ipf, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
ipf_final_emt_score <- final$emt_score %>% mean()

AT2_Score <- c(control_final_emt_score, covid_final_emt_score, ipf_final_emt_score)
Diagnosis <- c("Control", "COVID-19", "IPF")
AT2_Score <- data.frame( Diagnosis,AT2_Score)
AT2_Score



# Monocytes 
control <- subset(Monocytes, subset = Diagnosis == "Control")
covid <- subset(Monocytes, subset = Diagnosis == "COVID")
ipf <- subset(Monocytes, subset = Diagnosis == "IPF")
DefaultAssay(control) <- "RNA"
DefaultAssay(covid) <- "RNA"
DefaultAssay(ipf) <- "RNA"

# Control 
exp_data <- FetchData(object = control, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
control_final_emt_score <- final$emt_score %>% mean()

# COVID-19 
exp_data <- FetchData(object = covid, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
covid_final_emt_score <- final$emt_score %>% mean()

# IPF
exp_data <- FetchData(object = ipf, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
ipf_final_emt_score <- final$emt_score %>% mean()

Monocytes_Score <- c(control_final_emt_score, covid_final_emt_score, ipf_final_emt_score)
Diagnosis <- c("Control", "COVID-19", "IPF")
Monocytes_Score <- data.frame( Diagnosis,Monocytes_Score)
Monocytes_Score

# Endothelial 
control <- subset(Endothelial, subset = Diagnosis == "Control")
covid <- subset(Endothelial, subset = Diagnosis == "COVID")
ipf <- subset(Endothelial, subset = Diagnosis == "IPF")
DefaultAssay(control) <- "RNA"
DefaultAssay(covid) <- "RNA"
DefaultAssay(ipf) <- "RNA"

# Control 
exp_data <- FetchData(object = control, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
control_final_emt_score <- final$emt_score %>% mean()

# COVID-19 
exp_data <- FetchData(object = covid, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
covid_final_emt_score <- final$emt_score %>% mean()

# IPF
exp_data <- FetchData(object = ipf, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
ipf_final_emt_score <- final$emt_score %>% mean()

Endothelial_Score <- c(control_final_emt_score, covid_final_emt_score, ipf_final_emt_score)
Diagnosis <- c("Control", "COVID-19", "IPF")
Endothelial_Score <- data.frame( Diagnosis,Endothelial_Score)
Endothelial_Score

# Macrophages 
control <- subset(Macrophages, subset = Diagnosis == "Control")
covid <- subset(Macrophages, subset = Diagnosis == "COVID")
ipf <- subset(Macrophages, subset = Diagnosis == "IPF")
DefaultAssay(control) <- "RNA"
DefaultAssay(covid) <- "RNA"
DefaultAssay(ipf) <- "RNA"

# Control 
exp_data <- FetchData(object = control, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
control_final_emt_score <- final$emt_score %>% mean()

# COVID-19 
exp_data <- FetchData(object = covid, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
covid_final_emt_score <- final$emt_score %>% mean()

# IPF
exp_data <- FetchData(object = ipf, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
ipf_final_emt_score <- final$emt_score %>% mean()

Macrophages_Score <- c(control_final_emt_score, covid_final_emt_score, ipf_final_emt_score)
Diagnosis <- c("Control", "COVID-19", "IPF")
Macrophages_Score <- data.frame( Diagnosis,Macrophages_Score)
Macrophages_Score

# NK 
control <- subset(NK, subset = Diagnosis == "Control")
covid <- subset(NK, subset = Diagnosis == "COVID")
ipf <- subset(NK, subset = Diagnosis == "IPF")
DefaultAssay(control) <- "RNA"
DefaultAssay(covid) <- "RNA"
DefaultAssay(ipf) <- "RNA"

# Control 
exp_data <- FetchData(object = control, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
control_final_emt_score <- final$emt_score %>% mean()

# COVID-19 
exp_data <- FetchData(object = covid, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
covid_final_emt_score <- final$emt_score %>% mean()

# IPF
exp_data <- FetchData(object = ipf, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
ipf_final_emt_score <- final$emt_score %>% mean()

NK_Score <- c(control_final_emt_score, covid_final_emt_score, ipf_final_emt_score)
Diagnosis <- c("Control", "COVID-19", "IPF")
NK_Score <- data.frame( Diagnosis,NK_Score)
NK_Score


# Ciliated 
exp_data <- FetchData(object = Ciliated, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
Ciliated_score <- final$emt_score %>% mean()

control <- subset(Ciliated, subset = Diagnosis == "Control")
covid <- subset(Ciliated, subset = Diagnosis == "COVID")
ipf <- subset(Ciliated, subset = Diagnosis == "IPF")
DefaultAssay(control) <- "RNA"
DefaultAssay(covid) <- "RNA"
DefaultAssay(ipf) <- "RNA"

# Control 
exp_data <- FetchData(object = control, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
control_final_emt_score <- final$emt_score %>% mean()

# COVID-19 
exp_data <- FetchData(object = covid, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
covid_final_emt_score <- final$emt_score %>% mean()

# IPF
exp_data <- FetchData(object = ipf, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
ipf_final_emt_score <- final$emt_score %>% mean()

Ciliated_Score <- c(control_final_emt_score, covid_final_emt_score, ipf_final_emt_score)
Diagnosis <- c("Control", "COVID-19", "IPF")
Ciliated_Score <- data.frame( Diagnosis,Ciliated_Score)
Ciliated_Score


# Club 
control <- subset(Club, subset = Diagnosis == "Control")
covid <- subset(Club, subset = Diagnosis == "COVID")
ipf <- subset(Club, subset = Diagnosis == "IPF")
DefaultAssay(control) <- "RNA"
DefaultAssay(covid) <- "RNA"
DefaultAssay(ipf) <- "RNA"

# Control 
exp_data <- FetchData(object = control, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
control_final_emt_score <- final$emt_score %>% mean()

# COVID-19 
exp_data <- FetchData(object = covid, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
covid_final_emt_score <- final$emt_score %>% mean()

# IPF
exp_data <- FetchData(object = ipf, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
ipf_final_emt_score <- final$emt_score %>% mean()

Club_Score <- c(control_final_emt_score, covid_final_emt_score, ipf_final_emt_score)
Diagnosis <- c("Control", "COVID-19", "IPF")
Club_Score <- data.frame( Diagnosis,Club_Score)
Club_Score

# mDCs 
control <- subset(mDCs, subset = Diagnosis == "Control")
covid <- subset(mDCs, subset = Diagnosis == "COVID")
ipf <- subset(mDCs, subset = Diagnosis == "IPF")
DefaultAssay(control) <- "RNA"
DefaultAssay(covid) <- "RNA"
DefaultAssay(ipf) <- "RNA"

# Control 
exp_data <- FetchData(object = control, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
control_final_emt_score <- final$emt_score %>% mean()

# COVID-19 
exp_data <- FetchData(object = covid, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
covid_final_emt_score <- final$emt_score %>% mean()

# IPF
exp_data <- FetchData(object = ipf, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
ipf_final_emt_score <- final$emt_score %>% mean()

mDCs_Score <- c(control_final_emt_score, covid_final_emt_score, ipf_final_emt_score)
Diagnosis <- c("Control", "COVID-19", "IPF")
mDCs_Score <- data.frame( Diagnosis,mDCs_Score)
mDCs_Score


# Bcell 
control <- subset(Bcell, subset = Diagnosis == "Control")
covid <- subset(Bcell, subset = Diagnosis == "COVID")
ipf <- subset(Bcell, subset = Diagnosis == "IPF")
DefaultAssay(control) <- "RNA"
DefaultAssay(covid) <- "RNA"
DefaultAssay(ipf) <- "RNA"

# Control 
exp_data <- FetchData(object = control, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
control_final_emt_score <- final$emt_score %>% mean()

# COVID-19 
exp_data <- FetchData(object = covid, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
covid_final_emt_score <- final$emt_score %>% mean()

# IPF
exp_data <- FetchData(object = ipf, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
ipf_final_emt_score <- final$emt_score %>% mean()

Bcell_Score <- c(control_final_emt_score, covid_final_emt_score, ipf_final_emt_score)
Diagnosis <- c("Control", "COVID-19", "IPF")
Bcell_Score <- data.frame( Diagnosis,Bcell_Score)
Bcell_Score


# AT1 
control <- subset(AT1, subset = Diagnosis == "Control")
covid <- subset(AT1, subset = Diagnosis == "COVID")
ipf <- subset(AT1, subset = Diagnosis == "IPF")
DefaultAssay(control) <- "RNA"
DefaultAssay(covid) <- "RNA"
DefaultAssay(ipf) <- "RNA"

# Control 
exp_data <- FetchData(object = control, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
control_final_emt_score <- final$emt_score %>% mean()

# COVID-19 
exp_data <- FetchData(object = covid, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
covid_final_emt_score <- final$emt_score %>% mean()

# IPF
exp_data <- FetchData(object = ipf, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
ipf_final_emt_score <- final$emt_score %>% mean()

AT1_Score <- c(control_final_emt_score, covid_final_emt_score, ipf_final_emt_score)
Diagnosis <- c("Control", "COVID-19", "IPF")
AT1_Score <- data.frame( Diagnosis,AT1_Score)
AT1_Score


# Fibroblasts 
control <- subset(Fibroblasts, subset = Diagnosis == "Control")
covid <- subset(Fibroblasts, subset = Diagnosis == "COVID")
ipf <- subset(Fibroblasts, subset = Diagnosis == "IPF")
DefaultAssay(control) <- "RNA"
DefaultAssay(covid) <- "RNA"
DefaultAssay(ipf) <- "RNA"

# Control 
exp_data <- FetchData(object = control, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
control_final_emt_score <- final$emt_score %>% mean()

# COVID-19 
exp_data <- FetchData(object = covid, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
covid_final_emt_score <- final$emt_score %>% mean()

# IPF
exp_data <- FetchData(object = ipf, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
ipf_final_emt_score <- final$emt_score %>% mean()

Fibroblasts_Score <- c(control_final_emt_score, covid_final_emt_score, ipf_final_emt_score)
Diagnosis <- c("Control", "COVID-19", "IPF")
Fibroblasts_Score <- data.frame( Diagnosis,Fibroblasts_Score)
Fibroblasts_Score



# Basal 
control <- subset(Basal, subset = Diagnosis == "Control")
covid <- subset(Basal, subset = Diagnosis == "COVID")
ipf <- subset(Basal, subset = Diagnosis == "IPF")
DefaultAssay(control) <- "RNA"
DefaultAssay(covid) <- "RNA"
DefaultAssay(ipf) <- "RNA"

# Control 
exp_data <- FetchData(object = control, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
control_final_emt_score <- final$emt_score %>% mean()

# COVID-19 
exp_data <- FetchData(object = covid, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
covid_final_emt_score <- final$emt_score %>% mean()

# IPF
exp_data <- FetchData(object = ipf, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
ipf_final_emt_score <- final$emt_score %>% mean()

Basal_Score <- c(control_final_emt_score, covid_final_emt_score, ipf_final_emt_score)
Diagnosis <- c("Control", "COVID-19", "IPF")
Basal_Score <- data.frame( Diagnosis,Basal_Score)
Basal_Score


# Plasma 
control <- subset(Plasma, subset = Diagnosis == "Control")
covid <- subset(Plasma, subset = Diagnosis == "COVID")
ipf <- subset(Plasma, subset = Diagnosis == "IPF")
DefaultAssay(control) <- "RNA"
DefaultAssay(covid) <- "RNA"
DefaultAssay(ipf) <- "RNA"

# Control 
exp_data <- FetchData(object = control, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
control_final_emt_score <- final$emt_score %>% mean()

# COVID-19 
exp_data <- FetchData(object = covid, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
covid_final_emt_score <- final$emt_score %>% mean()

# IPF
exp_data <- FetchData(object = ipf, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
ipf_final_emt_score <- final$emt_score %>% mean()

Plasma_Score <- c(control_final_emt_score, covid_final_emt_score, ipf_final_emt_score)
Diagnosis <- c("Control", "COVID-19", "IPF")
Plasma_Score <- data.frame( Diagnosis,Plasma_Score)
Plasma_Score



# Gobelt 
control <- subset(Gobelt, subset = Diagnosis == "Control")
covid <- subset(Gobelt, subset = Diagnosis == "COVID")
ipf <- subset(Gobelt, subset = Diagnosis == "IPF")
DefaultAssay(control) <- "RNA"
DefaultAssay(covid) <- "RNA"
DefaultAssay(ipf) <- "RNA"

# Control 
exp_data <- FetchData(object = control, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
control_final_emt_score <- final$emt_score %>% mean()

# COVID-19 
exp_data <- FetchData(object = covid, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
covid_final_emt_score <- final$emt_score %>% mean()

# IPF
exp_data <- FetchData(object = ipf, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
ipf_final_emt_score <- final$emt_score %>% mean()

Gobelt_Score <- c(control_final_emt_score, covid_final_emt_score, ipf_final_emt_score)
Diagnosis <- c("Control", "COVID-19", "IPF")
Gobelt_Score <- data.frame( Diagnosis,Gobelt_Score)
Gobelt_Score

# CD4 
control <- subset(CD4, subset = Diagnosis == "Control")
covid <- subset(CD4, subset = Diagnosis == "COVID")
ipf <- subset(CD4, subset = Diagnosis == "IPF")
DefaultAssay(control) <- "RNA"
DefaultAssay(covid) <- "RNA"
DefaultAssay(ipf) <- "RNA"

# Control 
exp_data <- FetchData(object = control, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
control_final_emt_score <- final$emt_score %>% mean()

# COVID-19 
exp_data <- FetchData(object = covid, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
covid_final_emt_score <- final$emt_score %>% mean()

# IPF
exp_data <- FetchData(object = ipf, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
ipf_final_emt_score <- final$emt_score %>% mean()

CD4_Score <- c(control_final_emt_score, covid_final_emt_score, ipf_final_emt_score)
Diagnosis <- c("Control", "COVID-19", "IPF")
CD4_Score <- data.frame( Diagnosis,CD4_Score)
CD4_Score


# Mast 
control <- subset(Mast, subset = Diagnosis == "Control")
covid <- subset(Mast, subset = Diagnosis == "COVID")
ipf <- subset(Mast, subset = Diagnosis == "IPF")
DefaultAssay(control) <- "RNA"
DefaultAssay(covid) <- "RNA"
DefaultAssay(ipf) <- "RNA"

# Control 
exp_data <- FetchData(object = control, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
control_final_emt_score <- final$emt_score %>% mean()

# COVID-19 
exp_data <- FetchData(object = covid, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
covid_final_emt_score <- final$emt_score %>% mean()

# IPF
exp_data <- FetchData(object = ipf, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
ipf_final_emt_score <- final$emt_score %>% mean()

Mast_Score <- c(control_final_emt_score, covid_final_emt_score, ipf_final_emt_score)
Diagnosis <- c("Control", "COVID-19", "IPF")
Mast_Score <- data.frame( Diagnosis,Mast_Score)
Mast_Score

# Smooth 
control <- subset(Smooth, subset = Diagnosis == "Control")
covid <- subset(Smooth, subset = Diagnosis == "COVID")
ipf <- subset(Smooth, subset = Diagnosis == "IPF")
DefaultAssay(control) <- "RNA"
DefaultAssay(covid) <- "RNA"
DefaultAssay(ipf) <- "RNA"

# Control 
exp_data <- FetchData(object = control, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
control_final_emt_score <- final$emt_score %>% mean()

# COVID-19 
exp_data <- FetchData(object = covid, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
covid_final_emt_score <- final$emt_score %>% mean()

# IPF
exp_data <- FetchData(object = ipf, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
ipf_final_emt_score <- final$emt_score %>% mean()

Smooth_Score <- c(control_final_emt_score, covid_final_emt_score, ipf_final_emt_score)
Diagnosis <- c("Control", "COVID-19", "IPF")
Smooth_Score <- data.frame( Diagnosis,Smooth_Score)
Smooth_Score


# CD8 
control <- subset(CD8, subset = Diagnosis == "Control")
covid <- subset(CD8, subset = Diagnosis == "COVID")
ipf <- subset(CD8, subset = Diagnosis == "IPF")
DefaultAssay(control) <- "RNA"
DefaultAssay(covid) <- "RNA"
DefaultAssay(ipf) <- "RNA"

# Control 
exp_data <- FetchData(object = control, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
control_final_emt_score <- final$emt_score %>% mean()

# COVID-19 
exp_data <- FetchData(object = covid, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
covid_final_emt_score <- final$emt_score %>% mean()

# IPF
exp_data <- FetchData(object = ipf, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
ipf_final_emt_score <- final$emt_score %>% mean()

CD8_Score <- c(control_final_emt_score, covid_final_emt_score, ipf_final_emt_score)
Diagnosis <- c("Control", "COVID-19", "IPF")
CD8_Score <- data.frame( Diagnosis,CD8_Score)
CD8_Score


# pDCs 
control <- subset(pDCs, subset = Diagnosis == "Control")
covid <- subset(pDCs, subset = Diagnosis == "COVID")
ipf <- subset(pDCs, subset = Diagnosis == "IPF")
DefaultAssay(control) <- "RNA"
DefaultAssay(covid) <- "RNA"
DefaultAssay(ipf) <- "RNA"

# Control 
exp_data <- FetchData(object = control, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
control_final_emt_score <- final$emt_score %>% mean()

# COVID-19 
exp_data <- FetchData(object = covid, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
covid_final_emt_score <- final$emt_score %>% mean()

# IPF
exp_data <- FetchData(object = ipf, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
ipf_final_emt_score <- final$emt_score %>% mean()

pDCs_Score <- c(control_final_emt_score, covid_final_emt_score, ipf_final_emt_score)
Diagnosis <- c("Control", "COVID-19", "IPF")
pDCs_Score <- data.frame( Diagnosis,pDCs_Score)
pDCs_Score


# Neutrophils 
control <- subset(Neutrophils, subset = Diagnosis == "Control")
covid <- subset(Neutrophils, subset = Diagnosis == "COVID")
ipf <- subset(Neutrophils, subset = Diagnosis == "IPF")
DefaultAssay(control) <- "RNA"
DefaultAssay(covid) <- "RNA"
DefaultAssay(ipf) <- "RNA"

# Control 
exp_data <- FetchData(object = control, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
control_final_emt_score <- final$emt_score %>% mean()

# COVID-19 
exp_data <- FetchData(object = covid, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
covid_final_emt_score <- final$emt_score %>% mean()

# IPF
exp_data <- FetchData(object = ipf, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
ipf_final_emt_score <- final$emt_score %>% mean()

Neutrophils_Score <- c(control_final_emt_score, covid_final_emt_score, ipf_final_emt_score)
Neutrophils_Score
Diagnosis <- c("Control", "COVID-19", "IPF")
Neutrophils_Score <- data.frame( Diagnosis,Neutrophils_Score)
Neutrophils_Score


test <- left_join(AT2_Score, Monocytes_Score, by = "Diagnosis")
test <- left_join(test, Endothelial_Score, by = "Diagnosis")
test <- left_join(test, Macrophages_Score, by = "Diagnosis")
test <- left_join(test, NK_Score, by = "Diagnosis")
test <- left_join(test, Ciliated_Score, by = "Diagnosis")
test <- left_join(test, Club_Score, by = "Diagnosis")
test <- left_join(test, mDCs_Score, by = "Diagnosis")
test <- left_join(test, Bcell_Score, by = "Diagnosis")
test <- left_join(test, AT1_Score, by = "Diagnosis")
test <- left_join(test, Fibroblasts_Score, by = "Diagnosis")
test <- left_join(test, Basal_Score, by = "Diagnosis")
test <- left_join(test, Plasma_Score, by = "Diagnosis")
test <- left_join(test, Gobelt_Score, by = "Diagnosis")
test <- left_join(test, CD4_Score, by = "Diagnosis")
test <- left_join(test, Mast_Score, by = "Diagnosis")
test <- left_join(test, Smooth_Score, by = "Diagnosis")
test <- left_join(test, CD8_Score, by = "Diagnosis")
test <- left_join(test, pDCs_Score, by = "Diagnosis")
test <- left_join(test, Neutrophils_Score, by = "Diagnosis")


write_csv(test, here("final/csv/fibro_related_expression_check_test2.csv"))

test
Score <- c(AT2_Score, Monocytes_Score, Endothelial_Score, Macrophages_Score, NK_Score,
           Ciliated_Score, Club_Score, mDCs_Score, Bcell_Score, AT1_Score, Fibroblasts_Score, Basal_Score,
           Plasma_Score, Gobelt_Score, CD4_Score, Mast_Score, Smooth_Score, CD8_Score, pDCs_Score, Neutrophils_Score)
Score

celltype <- c("AT2_score", "Monocytes_score", "Endothelial_score", "Macrophages_score", "NK_score",
           "Ciliated_score", "Club_score", "mDCs_score", "Bcell_score", "AT1_score", "Fibroblasts_score", "Basal_score",
           "Plasma_score", "Gobelt_score", "CD4_score", "Mast_score", "Smooth_score", "CD8_score", "pDCs_score", "Neutrophils_score")
#Diagnosis <- c("Control", "COVID-19", "IPF")

merge <- data.frame(Score, celltype)
merge
write_csv(merge, here("final/csv/fibro_related_expression_check2.csv"))

Score
merge
