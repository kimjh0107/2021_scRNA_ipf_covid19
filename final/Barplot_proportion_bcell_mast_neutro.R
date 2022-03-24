library(here)
library(tidyverse)
library(Seurat)

macro <- readRDS(here("final/RDS/bcell.RDS"))


#### scRNA ####
# lung proportion 
cells <- macro@meta.data

cell_counts <- cells %>% 
  dplyr::select(group, Diagnosis_tag, seurat_clusters) %>%
  group_by(group, Diagnosis_tag, seurat_clusters) %>%
  summarize(total_cell_counts = n(), .groups = 'drop') 

# Control 
control_cell_count_1 <- cell_counts %>% 
  filter(Diagnosis_tag == 'HC1') %>% 
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

control_cell_count_2 <- cell_counts %>% 
  filter(Diagnosis_tag == 'HC2') %>% 
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

control_cell_count_3 <- cell_counts %>% 
  filter(Diagnosis_tag == 'HC3') %>% 
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

control_cell_count_4 <- cell_counts %>% 
  filter(Diagnosis_tag == 'HC4') %>% 
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

control_cell_count_5 <- cell_counts %>% 
  filter(Diagnosis_tag == 'HC5') %>% 
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

control_cell_count_6 <- cell_counts %>% 
  filter(Diagnosis_tag == 'HC6') %>% 
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

control_cell_count_7 <- cell_counts %>% 
  filter(Diagnosis_tag == 'HC7') %>% 
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

control_cell_count_8 <- cell_counts %>% 
  filter(Diagnosis_tag == 'HC8') %>% 
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

control_cell_count_9 <- cell_counts %>% 
  filter(Diagnosis_tag == 'HC9') %>% 
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

control_cell_count_10 <- cell_counts %>% 
  filter(Diagnosis_tag == 'HC10') %>% 
  mutate(proportion = total_cell_counts / sum(total_cell_counts))


# COVID
covid_cell_count_1 <- cell_counts %>% 
  filter(Diagnosis_tag == 'C1_1') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

covid_cell_count_2 <- cell_counts %>% 
  filter(Diagnosis_tag == 'C1_2') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

covid_cell_count_3 <- cell_counts %>% 
  filter(Diagnosis_tag == 'C1_3') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

covid_cell_count_4 <- cell_counts %>% 
  filter(Diagnosis_tag == 'C1_4') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))


# IPF
ipf_cell_count_1 <- cell_counts %>% 
  filter(Diagnosis_tag == 'I1') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))    

ipf_cell_count_2 <- cell_counts %>% 
  filter(Diagnosis_tag == 'I2') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

ipf_cell_count_3 <- cell_counts %>% 
  filter(Diagnosis_tag == 'I3') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

ipf_cell_count_4 <- cell_counts %>% 
  filter(Diagnosis_tag == 'I4') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

ipf_cell_count_5 <- cell_counts %>% 
  filter(Diagnosis_tag == 'I5') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

ipf_cell_count_6 <- cell_counts %>% 
  filter(Diagnosis_tag == 'I6') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

ipf_cell_count_7 <- cell_counts %>% 
  filter(Diagnosis_tag == 'I7') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

ipf_cell_count_8 <- cell_counts %>% 
  filter(Diagnosis_tag == 'I8') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

ipf_cell_count_9 <- cell_counts %>% 
  filter(Diagnosis_tag == 'I9') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

ipf_cell_count_10 <- cell_counts %>% 
  filter(Diagnosis_tag == 'I10') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

ipf_cell_count_11 <- cell_counts %>% 
  filter(Diagnosis_tag == 'I11') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

ipf_cell_count_12 <- cell_counts %>% 
  filter(Diagnosis_tag == 'I12') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

ipf_cell_count_13 <- cell_counts %>% 
  filter(Diagnosis_tag == 'I13') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

ipf_cell_count_14 <- cell_counts %>% 
  filter(Diagnosis_tag == 'I14') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

ipf_cell_count_15 <- cell_counts %>% 
  filter(Diagnosis_tag == 'I15') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

ipf_cell_count_16 <- cell_counts %>% 
  filter(Diagnosis_tag == 'I16') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

ipf_cell_count_17 <- cell_counts %>% 
  filter(Diagnosis_tag == 'I17') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

ipf_cell_count_18 <- cell_counts %>% 
  filter(Diagnosis_tag == 'I18') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

ipf_cell_count_19 <- cell_counts %>% 
  filter(Diagnosis_tag == 'I19') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

ipf_cell_count_20 <- cell_counts %>% 
  filter(Diagnosis_tag == 'I20') %>%
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
                         ipf_cell_count_12, 
                         ipf_cell_count_13, 
                         ipf_cell_count_14, 
                         ipf_cell_count_15, 
                         ipf_cell_count_16, 
                         ipf_cell_count_17, 
                         ipf_cell_count_18, 
                         ipf_cell_count_19, 
                         ipf_cell_count_20)


# ggplot - celltype
ggplot(macro@meta.data, aes(Diagnosis_tag, fill = seurat_clusters)) + 
  geom_bar(position='fill') + 
  #scale_y_log10() + 
  #coord_flip() + 
  labs(y = '') + 
  theme_classic() + 
  RotatedAxis() + 
  scale_x_discrete(limits = c("HC1","HC2","HC3","HC4","HC5","HC6","HC7","HC8","HC9","HC10",
                              "C1_1","C1_2","C1_3","C1_4",
                              "I1","I2","I3","I4","I5","I6","I7","I8","I9","I10",
                              "I11","I12","I13","I14","I15","I16","I17","I18","I19","I20"
  )) + 
  theme(legend.position = "right",axis.title.y = element_blank())
ggsave(here("final/figure/Barplot_Bcell_Proportion.pdf"), width = 7, height = 5)







library(here)
library(tidyverse)
library(Seurat)

macro <- readRDS(here("final/RDS/mast.RDS"))

#### scRNA ####
# lung proportion 
cells <- macro@meta.data

cell_counts <- cells %>% 
  dplyr::select(group, Diagnosis_tag, seurat_clusters) %>%
  group_by(group, Diagnosis_tag, seurat_clusters) %>%
  summarize(total_cell_counts = n(), .groups = 'drop') 

# Control 
control_cell_count_1 <- cell_counts %>% 
  filter(Diagnosis_tag == 'HC1') %>% 
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

control_cell_count_2 <- cell_counts %>% 
  filter(Diagnosis_tag == 'HC2') %>% 
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

control_cell_count_3 <- cell_counts %>% 
  filter(Diagnosis_tag == 'HC3') %>% 
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

control_cell_count_4 <- cell_counts %>% 
  filter(Diagnosis_tag == 'HC4') %>% 
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

control_cell_count_5 <- cell_counts %>% 
  filter(Diagnosis_tag == 'HC5') %>% 
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

control_cell_count_6 <- cell_counts %>% 
  filter(Diagnosis_tag == 'HC6') %>% 
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

control_cell_count_7 <- cell_counts %>% 
  filter(Diagnosis_tag == 'HC7') %>% 
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

control_cell_count_8 <- cell_counts %>% 
  filter(Diagnosis_tag == 'HC8') %>% 
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

control_cell_count_9 <- cell_counts %>% 
  filter(Diagnosis_tag == 'HC9') %>% 
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

control_cell_count_10 <- cell_counts %>% 
  filter(Diagnosis_tag == 'HC10') %>% 
  mutate(proportion = total_cell_counts / sum(total_cell_counts))


# COVID
covid_cell_count_1 <- cell_counts %>% 
  filter(Diagnosis_tag == 'C1_1') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

covid_cell_count_2 <- cell_counts %>% 
  filter(Diagnosis_tag == 'C1_2') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

covid_cell_count_3 <- cell_counts %>% 
  filter(Diagnosis_tag == 'C1_3') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

covid_cell_count_4 <- cell_counts %>% 
  filter(Diagnosis_tag == 'C1_4') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))


# IPF
ipf_cell_count_1 <- cell_counts %>% 
  filter(Diagnosis_tag == 'I1') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))    

ipf_cell_count_2 <- cell_counts %>% 
  filter(Diagnosis_tag == 'I2') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

ipf_cell_count_3 <- cell_counts %>% 
  filter(Diagnosis_tag == 'I3') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

ipf_cell_count_4 <- cell_counts %>% 
  filter(Diagnosis_tag == 'I4') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

ipf_cell_count_5 <- cell_counts %>% 
  filter(Diagnosis_tag == 'I5') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

ipf_cell_count_6 <- cell_counts %>% 
  filter(Diagnosis_tag == 'I6') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

ipf_cell_count_7 <- cell_counts %>% 
  filter(Diagnosis_tag == 'I7') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

ipf_cell_count_8 <- cell_counts %>% 
  filter(Diagnosis_tag == 'I8') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

ipf_cell_count_9 <- cell_counts %>% 
  filter(Diagnosis_tag == 'I9') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

ipf_cell_count_10 <- cell_counts %>% 
  filter(Diagnosis_tag == 'I10') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

ipf_cell_count_11 <- cell_counts %>% 
  filter(Diagnosis_tag == 'I11') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

ipf_cell_count_12 <- cell_counts %>% 
  filter(Diagnosis_tag == 'I12') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

ipf_cell_count_13 <- cell_counts %>% 
  filter(Diagnosis_tag == 'I13') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

ipf_cell_count_14 <- cell_counts %>% 
  filter(Diagnosis_tag == 'I14') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

ipf_cell_count_15 <- cell_counts %>% 
  filter(Diagnosis_tag == 'I15') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

ipf_cell_count_16 <- cell_counts %>% 
  filter(Diagnosis_tag == 'I16') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

ipf_cell_count_17 <- cell_counts %>% 
  filter(Diagnosis_tag == 'I17') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

ipf_cell_count_18 <- cell_counts %>% 
  filter(Diagnosis_tag == 'I18') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

ipf_cell_count_19 <- cell_counts %>% 
  filter(Diagnosis_tag == 'I19') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

ipf_cell_count_20 <- cell_counts %>% 
  filter(Diagnosis_tag == 'I20') %>%
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
                         ipf_cell_count_12, 
                         ipf_cell_count_13, 
                         ipf_cell_count_14, 
                         ipf_cell_count_15, 
                         ipf_cell_count_16, 
                         ipf_cell_count_17, 
                         ipf_cell_count_18, 
                         ipf_cell_count_19, 
                         ipf_cell_count_20)


# ggplot - celltype
ggplot(macro@meta.data, aes(Diagnosis_tag, fill = seurat_clusters)) + 
  geom_bar(position='fill') + 
  #scale_y_log10() + 
  #coord_flip() + 
  labs(y = '') + 
  theme_classic() + 
  RotatedAxis() + 
  scale_x_discrete(limits = c("HC1","HC2","HC3","HC4","HC5","HC6","HC7","HC8","HC9","HC10",
                              "C1_1","C1_2","C1_3","C1_4",
                              "I1","I2","I3","I4","I5","I6","I7","I8","I9","I10",
                              "I11","I12","I13","I14","I15","I16","I17","I18","I19","I20"
  )) + 
  theme(legend.position = "right",axis.title.y = element_blank())
ggsave(here("final/figure/Barplot_Mast_Proportion.pdf"), width = 7, height = 5)







library(here)
library(tidyverse)
library(Seurat)

macro <- readRDS(here("final/RDS/neutrophil.RDS"))

#### scRNA ####
# lung proportion 
cells <- macro@meta.data

cell_counts <- cells %>% 
  dplyr::select(group, Diagnosis_tag, seurat_clusters) %>%
  group_by(group, Diagnosis_tag, seurat_clusters) %>%
  summarize(total_cell_counts = n(), .groups = 'drop') 

# Control 
control_cell_count_1 <- cell_counts %>% 
  filter(Diagnosis_tag == 'HC1') %>% 
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

control_cell_count_2 <- cell_counts %>% 
  filter(Diagnosis_tag == 'HC2') %>% 
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

control_cell_count_3 <- cell_counts %>% 
  filter(Diagnosis_tag == 'HC3') %>% 
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

control_cell_count_4 <- cell_counts %>% 
  filter(Diagnosis_tag == 'HC4') %>% 
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

control_cell_count_5 <- cell_counts %>% 
  filter(Diagnosis_tag == 'HC5') %>% 
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

control_cell_count_6 <- cell_counts %>% 
  filter(Diagnosis_tag == 'HC6') %>% 
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

control_cell_count_7 <- cell_counts %>% 
  filter(Diagnosis_tag == 'HC7') %>% 
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

control_cell_count_8 <- cell_counts %>% 
  filter(Diagnosis_tag == 'HC8') %>% 
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

control_cell_count_9 <- cell_counts %>% 
  filter(Diagnosis_tag == 'HC9') %>% 
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

control_cell_count_10 <- cell_counts %>% 
  filter(Diagnosis_tag == 'HC10') %>% 
  mutate(proportion = total_cell_counts / sum(total_cell_counts))


# COVID
covid_cell_count_1 <- cell_counts %>% 
  filter(Diagnosis_tag == 'C1_1') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

covid_cell_count_2 <- cell_counts %>% 
  filter(Diagnosis_tag == 'C1_2') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

covid_cell_count_3 <- cell_counts %>% 
  filter(Diagnosis_tag == 'C1_3') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

covid_cell_count_4 <- cell_counts %>% 
  filter(Diagnosis_tag == 'C1_4') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))


# IPF
ipf_cell_count_1 <- cell_counts %>% 
  filter(Diagnosis_tag == 'I1') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))    

ipf_cell_count_2 <- cell_counts %>% 
  filter(Diagnosis_tag == 'I2') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

ipf_cell_count_3 <- cell_counts %>% 
  filter(Diagnosis_tag == 'I3') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

ipf_cell_count_4 <- cell_counts %>% 
  filter(Diagnosis_tag == 'I4') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

ipf_cell_count_5 <- cell_counts %>% 
  filter(Diagnosis_tag == 'I5') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

ipf_cell_count_6 <- cell_counts %>% 
  filter(Diagnosis_tag == 'I6') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

ipf_cell_count_7 <- cell_counts %>% 
  filter(Diagnosis_tag == 'I7') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

ipf_cell_count_8 <- cell_counts %>% 
  filter(Diagnosis_tag == 'I8') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

ipf_cell_count_9 <- cell_counts %>% 
  filter(Diagnosis_tag == 'I9') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

ipf_cell_count_10 <- cell_counts %>% 
  filter(Diagnosis_tag == 'I10') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

ipf_cell_count_11 <- cell_counts %>% 
  filter(Diagnosis_tag == 'I11') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

ipf_cell_count_12 <- cell_counts %>% 
  filter(Diagnosis_tag == 'I12') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

ipf_cell_count_13 <- cell_counts %>% 
  filter(Diagnosis_tag == 'I13') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

ipf_cell_count_14 <- cell_counts %>% 
  filter(Diagnosis_tag == 'I14') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

ipf_cell_count_15 <- cell_counts %>% 
  filter(Diagnosis_tag == 'I15') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

ipf_cell_count_16 <- cell_counts %>% 
  filter(Diagnosis_tag == 'I16') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

ipf_cell_count_17 <- cell_counts %>% 
  filter(Diagnosis_tag == 'I17') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

ipf_cell_count_18 <- cell_counts %>% 
  filter(Diagnosis_tag == 'I18') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

ipf_cell_count_19 <- cell_counts %>% 
  filter(Diagnosis_tag == 'I19') %>%
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

ipf_cell_count_20 <- cell_counts %>% 
  filter(Diagnosis_tag == 'I20') %>%
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
                         ipf_cell_count_12, 
                         ipf_cell_count_13, 
                         ipf_cell_count_14, 
                         ipf_cell_count_15, 
                         ipf_cell_count_16, 
                         ipf_cell_count_17, 
                         ipf_cell_count_18, 
                         ipf_cell_count_19, 
                         ipf_cell_count_20)


# ggplot - celltype
ggplot(macro@meta.data, aes(Diagnosis_tag, fill = seurat_clusters)) + 
  geom_bar(position='fill') + 
  #scale_y_log10() + 
  #coord_flip() + 
  labs(y = '') + 
  theme_classic() + 
  RotatedAxis() + 
  scale_x_discrete(limits = c("HC1","HC2","HC3","HC4","HC5","HC6","HC7","HC8","HC9","HC10",
                              "C1_1","C1_2","C1_3","C1_4",
                              "I1","I2","I3","I4","I5","I6","I7","I8","I9","I10",
                              "I11","I12","I13","I14","I15","I16","I17","I18","I19","I20"
  )) + 
  theme(legend.position = "right",axis.title.y = element_blank())
ggsave(here("final/figure/Barplot_Neutrophils_Proportion.pdf"), width = 7, height = 5)

