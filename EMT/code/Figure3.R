library(Seurat)
library(tidyverse)
library(here)
library(presto)
library(msigdbr)
library(fgsea)
library(ggpubr)

at <- readRDS(here("EMT/RDS/AT.RDS"))
DefaultAssay(at) <- "RNA"
DimPlot(at, group.by = "new_submain", label = T, repel = T) + NoLegend()
ggsave(here("EMT/figure/Figure3/Fig4a_DimPlot_celltype.pdf"), width = 7, height = 5)

DimPlot(at, group.by = "Diagnosis", label = F, repel = T, shuffle = T, cols = c("#00A0FA", "red", "purple")) + ggplot2::theme(legend.position = "bottom")
ggsave(here("EMT/figure/Figure3/Fig4b_DimPlot_diagnosis.pdf"), width = 7, height = 5)


# VlnPlot 
at1 <- subset(at, subset = new_submain %in% c("AT1"))
at2 <- subset(at, subset = new_submain %in% c("AT2"))

incomple_transition_marker_gene <- c("ETV5")
late_AT1_maturation <- c("CAV1")

DefaultAssay(at1) <- "RNA"
DefaultAssay(at2) <- "RNA"
Idents(at1) <- "Diagnosis"
Idents(at2) <- "Diagnosis"

VlnPlot(at2, features = c("ETV5"), pt.size = 0.0) + NoLegend()
ggsave(here("EMT/figure/Figure3/at_vlnplot_ETV5.pdf"), width = 5, height = 7)

VlnPlot(at1, features = c("CAV1"), pt.size = 0.0) + NoLegend()
ggsave(here("EMT/figure/Figure3/at_vlnplot_CAV1.pdf"), width = 5, height = 7)


# Proportion  - AT1, AT2 
at_control <- subset(at, subset = Diagnosis == "Control")
at_covid <- subset(at, subset = Diagnosis == "COVID")
at_ipf <- subset(at, subset = Diagnosis == "IPF")
DefaultAssay(at_control) <- "RNA"
DefaultAssay(at_covid) <- "RNA"
DefaultAssay(at_ipf) <- "RNA"



### function for subset celltypes###
get_control_proportion <- function(diagnosis){
  diagnosis_cells <- diagnosis@meta.data
  cell_count_control <- diagnosis_cells %>% 
    select(new_submain) %>% 
    group_by(new_submain) %>% 
    summarise(count_cell = n(), .groups = "drop") %>% 
    mutate(total = sum(count_cell)) %>% 
    mutate(Control_Proportion = count_cell/total *100) %>% 
    select(new_submain, Control_Proportion)
}

get_covid_proportion <- function(diagnosis){
  diagnosis_cells <- diagnosis@meta.data
  cell_count_control <- diagnosis_cells %>% 
    select(new_submain) %>% 
    group_by(new_submain) %>% 
    summarise(count_cell = n(), .groups = "drop") %>% 
    mutate(total = sum(count_cell)) %>% 
    mutate(COVID19_Proportion = count_cell/total *100) %>% 
    select(new_submain, COVID19_Proportion)
}

get_ipf_proportion <- function(diagnosis){
  diagnosis_cells <- diagnosis@meta.data
  cell_count_control <- diagnosis_cells %>% 
    select(new_submain) %>% 
    group_by(new_submain) %>% 
    summarise(count_cell = n(), .groups = "drop") %>% 
    mutate(total = sum(count_cell)) %>% 
    mutate(IPF_Proportion = count_cell/total *100) %>% 
    select(new_submain, IPF_Proportion)
}


control_proportion <- get_control_proportion(at_control)
covid_proportion <- get_covid_proportion(at_covid)
ipf_proportion <- get_ipf_proportion(at_ipf)

merge <- left_join(control_proportion, covid_proportion, by = "new_submain")
merge <- left_join(merge, ipf_proportion, by = "new_submain")
write_csv(merge, here("EMT/csv/at_proportion_table.csv"))




#### Proportion ####
#### scRNA ####
# lung proportion 
cells <- at@meta.data

cell_counts <- cells %>% 
  dplyr::select(Diagnosis, Diagnosis_tag2, new_submain) %>%
  group_by(Diagnosis, Diagnosis_tag2, new_submain) %>%
  summarize(total_cell_counts = n(), .groups = 'drop') 

# Control 
control_cell_count_1 <- cell_counts %>% 
  filter(Diagnosis_tag2 == 'HC1') %>% 
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

control_cell_count_2 <- cell_counts %>% 
  filter(Diagnosis_tag2 == 'HC2') %>% 
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

control_cell_count_3 <- cell_counts %>% 
  filter(Diagnosis_tag2 == 'HC3') %>% 
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

control_cell_count_4 <- cell_counts %>% 
  filter(Diagnosis_tag2 == 'HC4') %>% 
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

control_cell_count_5 <- cell_counts %>% 
  filter(Diagnosis_tag2 == 'HC5') %>% 
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

control_cell_count_6 <- cell_counts %>% 
  filter(Diagnosis_tag2 == 'HC6') %>% 
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

control_cell_count_7 <- cell_counts %>% 
  filter(Diagnosis_tag2 == 'HC7') %>% 
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

control_cell_count_8 <- cell_counts %>% 
  filter(Diagnosis_tag2 == 'HC8') %>% 
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

control_cell_count_9 <- cell_counts %>% 
  filter(Diagnosis_tag2 == 'HC9') %>% 
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

control_cell_count_10 <- cell_counts %>% 
  filter(Diagnosis_tag2 == 'HC10') %>% 
  mutate(proportion = total_cell_counts / sum(total_cell_counts))


# COVID
covid_cell_count_1 <- cell_counts %>% 
  filter(Diagnosis_tag2 == 'C1_1') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

covid_cell_count_2 <- cell_counts %>% 
  filter(Diagnosis_tag2 == 'C1_2') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

covid_cell_count_3 <- cell_counts %>% 
  filter(Diagnosis_tag2 == 'C1_3') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

covid_cell_count_4 <- cell_counts %>% 
  filter(Diagnosis_tag2 == 'C1_4') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))


# IPF
ipf_cell_count_1 <- cell_counts %>% 
  filter(Diagnosis_tag2 == 'I1') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))    

ipf_cell_count_2 <- cell_counts %>% 
  filter(Diagnosis_tag2 == 'I2') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

ipf_cell_count_3 <- cell_counts %>% 
  filter(Diagnosis_tag2 == 'I3') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

ipf_cell_count_4 <- cell_counts %>% 
  filter(Diagnosis_tag2 == 'I4') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

ipf_cell_count_5 <- cell_counts %>% 
  filter(Diagnosis_tag2 == 'I5') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

ipf_cell_count_6 <- cell_counts %>% 
  filter(Diagnosis_tag2 == 'I6') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

ipf_cell_count_7 <- cell_counts %>% 
  filter(Diagnosis_tag2 == 'I7') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

ipf_cell_count_8 <- cell_counts %>% 
  filter(Diagnosis_tag2 == 'I8') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

ipf_cell_count_9 <- cell_counts %>% 
  filter(Diagnosis_tag2 == 'I9') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

ipf_cell_count_10 <- cell_counts %>% 
  filter(Diagnosis_tag2 == 'I10') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

ipf_cell_count_11 <- cell_counts %>% 
  filter(Diagnosis_tag2 == 'I11') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

ipf_cell_count_12 <- cell_counts %>% 
  filter(Diagnosis_tag2 == 'I12') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))



cell_counts <- bind_rows(control_cell_count_1, 
                         control_cell_count_2, 
                         control_cell_count_3, 
                         control_cell_count_4, 
                         control_cell_count_5, 
                         control_cell_count_6, 
                         control_cell_count_7, 
                         control_cell_count_8, 
                         control_cell_count_9, 
                         control_cell_count_10, 
                         covid_cell_count_1, 
                         covid_cell_count_2, 
                         covid_cell_count_3, 
                         covid_cell_count_4,
                         ipf_cell_count_1, 
                         ipf_cell_count_2, 
                         ipf_cell_count_3, 
                         ipf_cell_count_4, 
                         ipf_cell_count_5, 
                         ipf_cell_count_6, 
                         ipf_cell_count_7, 
                         ipf_cell_count_8, 
                         ipf_cell_count_9, 
                         ipf_cell_count_10, 
                         ipf_cell_count_11, 
                         ipf_cell_count_12)


# ggplot - celltype
ggplot(at@meta.data, aes(Diagnosis_tag2, fill = new_submain)) + 
  geom_bar(position='fill') + 
  #scale_y_log10() + 
  #coord_flip() + 
  labs(y = '') + 
  theme_classic() + 
  RotatedAxis() + 
  scale_x_discrete(limits = c("HC1","HC2","HC3","HC4","HC5","HC6","HC7","HC8","HC9","HC10",
                              "C1_1","C1_2","C1_3","C1_4",
                              "I1","I2","I3","I4","I5","I6","I7","I8","I9","I10",
                              "I11","I12"
                              )) + 
  theme(legend.position = "right",axis.title.y = element_blank()) + 
  scale_fill_brewer(palette = "Paired")

ggsave(here("EMT/figure/Figure3/at_proportion_barplot.pdf"), width = 7, height = 5)


