library(Seurat)
library(tidyverse)
library(here)
library(ggpubr)

df <- readRDS(here("final/RDS/df.RDS"))
DefaultAssay(df) <- "RNA"

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
gene_list <- c("COL1A1", "PDGFRA")

# function for emt scoring 
get_gene_emt_score <- function(gene){
  mean_exp <- exp_data %>% select({{gene}}) %>% pull() %>% mean()
  return (tibble(gene = gene, emt_score = mean_exp))
}

# AT2 
exp_data <- FetchData(object = AT2, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
AT2_score <- final$emt_score %>% mean()

# Monocytes 
exp_data <- FetchData(object = Monocytes, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
Monocytes_score <- final$emt_score %>% mean()

# Endothelial 
exp_data <- FetchData(object = Endothelial, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
Endothelial_score <- final$emt_score %>% mean()

# Macrophages 
exp_data <- FetchData(object = Macrophages, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
Macrophages_score <- final$emt_score %>% mean()

# NK 
exp_data <- FetchData(object = NK, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
NK_score <- final$emt_score %>% mean()

# Ciliated 
exp_data <- FetchData(object = Ciliated, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
Ciliated_score <- final$emt_score %>% mean()

# Club 
exp_data <- FetchData(object = Club, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
Club_score <- final$emt_score %>% mean()

# mDCs 
exp_data <- FetchData(object = mDCs, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
mDCs_score <- final$emt_score %>% mean()

# Bcell 
exp_data <- FetchData(object = Bcell, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
Bcell_score <- final$emt_score %>% mean()

# AT1 
exp_data <- FetchData(object = AT1, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
AT1_score <- final$emt_score %>% mean()

# Fibroblasts 
exp_data <- FetchData(object = Fibroblasts, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
Fibroblasts_score <- final$emt_score %>% mean()

# Basal 
exp_data <- FetchData(object = Basal, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
Basal_score <- final$emt_score %>% mean()

# Plasma 
exp_data <- FetchData(object = Plasma, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
Plasma_score <- final$emt_score %>% mean()

# Gobelt 
exp_data <- FetchData(object = Gobelt, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
Gobelt_score <- final$emt_score %>% mean()

# CD4 
exp_data <- FetchData(object = CD4, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
CD4_score <- final$emt_score %>% mean()

# Mast 
exp_data <- FetchData(object = Mast, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
Mast_score <- final$emt_score %>% mean()

# Smooth 
exp_data <- FetchData(object = Smooth, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
Smooth_score <- final$emt_score %>% mean()

# CD8 
exp_data <- FetchData(object = CD8, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
CD8_score <- final$emt_score %>% mean()

# pDCs 
exp_data <- FetchData(object = pDCs, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
pDCs_score <- final$emt_score %>% mean()

# Neutrophils 
exp_data <- FetchData(object = Neutrophils, vars = gene_list)
final <- colnames(exp_data) %>% map_df(get_gene_emt_score) 
Neutrophils_score <- final$emt_score %>% mean()


Score <- c(AT2_score, Monocytes_score, Endothelial_score, Macrophages_score, NK_score,
           Ciliated_score, Club_score, mDCs_score, Bcell_score, AT1_score, Fibroblasts_score, Basal_score,
           Plasma_score, Gobelt_score, CD4_score, Mast_score, Smooth_score, CD8_score, pDCs_score, Neutrophils_score)

celltype <- c("AT2_score", "Monocytes_score", "Endothelial_score", "Macrophages_score", "NK_score",
           "Ciliated_score", "Club_score", "mDCs_score", "Bcell_score", "AT1_score", "Fibroblasts_score", "Basal_score",
           "Plasma_score", "Gobelt_score", "CD4_score", "Mast_score", "Smooth_score", "CD8_score", "pDCs_score", "Neutrophils_score")
merge <- data.frame(celltype, Score)
write_csv(merge, here("final/csv/fibro_related_expression_check_2.csv"))

Score
merge
