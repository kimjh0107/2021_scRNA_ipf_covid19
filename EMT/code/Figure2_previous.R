library(here)
library(Seurat)
library(tidyverse)

fibro <- readRDS(here("EMT/RDS/fibro.RDS"))
DefaultAssay(fibro) <- "RNA"


DimPlot(fibro, group.by = "new_submain", label = T, repel = T) + NoLegend()
DimPlot(fibro, group.by = "seurat_clusters", label = T, repel = T) + NoLegend()

#ggsave(here("EMT/figure/Figure2/Fig1a_DimPlot_diagnosis.pdf"), width = 7, height = 5)
FeaturePlot(fibro, features = c("COL1A1"))
FeaturePlot(fibro, features = c("PDGFRA"))
FeaturePlot(fibro, features = c("ELN"))
FeaturePlot(fibro, features = c("ACTA2"))
VlnPlot(fibro, features = c("ACTA2") )


fibro@meta.data$new_submain2 <- ifelse(fibro@meta.data$seurat_clusters %in% c(8,7,9,6,11), "MyoFibroblasts", "NA")
fibro@meta.data$new_submain2 <- ifelse(fibro@meta.data$seurat_clusters %in% c(0,1,2,3,4,5,10,12,13), "Fibroblasts", fibro@meta.data$new_submain2)
DimPlot(fibro, group.by = "new_submain2", label = T, repel = T) + NoLegend()


control <- subset(fibro, subset = Diagnosis == "Control")
covid <- subset(fibro, subset = Diagnosis == "COVID")
ipf <- subset(fibro, subset = Diagnosis == "IPF")


### function for subset celltypes###
get_control_proportion <- function(diagnosis){
  diagnosis_cells <- diagnosis@meta.data
  cell_count_control <- diagnosis_cells %>% 
    select(new_submain2) %>% 
    group_by(new_submain2) %>% 
    summarise(count_cell = n(), .groups = "drop") %>% 
    mutate(total = sum(count_cell)) %>% 
    mutate(Control_Proportion = count_cell/total *100) %>% 
    select(new_submain2, Control_Proportion)
}

get_covid_proportion <- function(diagnosis){
  diagnosis_cells <- diagnosis@meta.data
  cell_count_control <- diagnosis_cells %>% 
    select(new_submain2) %>% 
    group_by(new_submain2) %>% 
    summarise(count_cell = n(), .groups = "drop") %>% 
    mutate(total = sum(count_cell)) %>% 
    mutate(COVID_Proportion = count_cell/total *100) %>% 
    select(new_submain2, COVID_Proportion)
}

get_ipf_proportion <- function(diagnosis){
  diagnosis_cells <- diagnosis@meta.data
  cell_count_control <- diagnosis_cells %>% 
    select(new_submain2) %>% 
    group_by(new_submain2) %>% 
    summarise(count_cell = n(), .groups = "drop") %>% 
    mutate(total = sum(count_cell)) %>% 
    mutate(IPF_Proportion = count_cell/total *100) %>% 
    select(new_submain2, IPF_Proportion)
}


# scRNA
control_proportion <- get_control_proportion(control)
covid_proportion <- get_covid_proportion(covid)
ipf_proportion <- get_ipf_proportion(ipf)

merge <- left_join(control_proportion, covid_proportion, by = "new_submain2")
merge <- left_join(merge, ipf_proportion, by = "new_submain2")
write_csv(merge, here("EMT/csv/Proportion_fibro_table.csv"))


