library(here)
library(Seurat)
library(tidyverse)

df <- readRDS(here("EMT/RDS/EMT.RDS"))
cells <- df@meta.data 
unique(cells$Diagnosis_tag2)

# lung proportion 
cells <- df@meta.data

cell_counts <- cells %>% 
  dplyr::select(Diagnosis, Diagnosis_tag2 ,new_submain) %>%
  group_by(Diagnosis, Diagnosis_tag2 ,new_submain) %>%
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
                         ipf_cell_count_10
                   )

# ggplot
ggplot(df@meta.data, aes(Diagnosis_tag2, fill = new_submain)) + 
  geom_bar(position='fill') + 
  #  scale_y_log10() + 
  #coord_flip() + 
  labs(y = '') + 
  theme_classic() + 
  RotatedAxis() + 
  scale_x_discrete(limits = c("HC1","HC2","HC3","HC4","HC5","HC6","HC7","HC8","HC9","HC10",
                              "C1_1","C1_2","C1_3","C1_4",
                              "I1","I2","I3","I4","I5","I6","I7","I8","I9","I10"
                              )) + 
  theme(legend.position = "right",axis.title.y = element_blank())
ggsave(here("EMT/figure/barplot_proportion_by_diagnosis_tag.pdf"), width = 7, height = 5)
