library(tidyverse)
library(Seurat)
library(here)


emt <- readRDS(here("EMT/RDS/fibro_annotation.RDS"))
emt_sn <- readRDS(here("EMT/RDS/snRNA_fibro_annotation.RDS"))



cells <- emt_sn@meta.data 

cell_counts <- cells %>% 
  dplyr::select(Diagnosis, Diagnosis_tag ,new_submain3) %>%
  group_by(Diagnosis, Diagnosis_tag ,new_submain3) %>%
  summarize(total_cell_counts = n(), .groups = 'drop') 

# Control 
control_cell_count <- cell_counts %>% 
  filter(Diagnosis == 'Control') %>% 
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

# Control 
covid_cell_count <- cell_counts %>% 
  filter(Diagnosis == 'COVID-19') %>% 
  mutate(proportion = total_cell_counts / sum(total_cell_counts))



cell_counts <- bind_rows(control_cell_count, covid_cell_count)
cell_counts

# ggplot
ggplot(emt_sn@meta.data, aes(Diagnosis, fill = new_submain3)) + 
  geom_bar(position='fill') + 
  scale_fill_manual(values = c("#56B4E9", "purple", "red")) + 
  labs(y = '') + 
  theme_classic() + 
  RotatedAxis() + 
  scale_x_discrete(limits = c("Control", "COVID-19")) + 
  theme(legend.position = "right",axis.title.y = element_blank())
ggsave(here("EMT/figure/barplot_proportion_fibro_annotate_by_diagnosis2.pdf"), width = 7, height = 5)








emt <- readRDS(here("EMT/RDS/fibro_annotation.RDS"))

cells <- emt@meta.data 

cell_counts <- cells %>% 
  dplyr::select(Diagnosis, Diagnosis_tag ,new_submain2) %>%
  group_by(Diagnosis, Diagnosis_tag ,new_submain2) %>%
  summarize(total_cell_counts = n(), .groups = 'drop') 

# Control 
control_cell_count <- cell_counts %>% 
  filter(Diagnosis == 'Control') %>% 
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

# Control 
covid_cell_count <- cell_counts %>% 
  filter(Diagnosis == 'COVID') %>% 
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

# Control 
ipf_cell_count <- cell_counts %>% 
  filter(Diagnosis == 'IPF') %>% 
  mutate(proportion = total_cell_counts / sum(total_cell_counts))

cell_counts <- bind_rows(control_cell_count, covid_cell_count, ipf_cell_count)
cell_counts

# ggplot
ggplot(emt@meta.data, aes(Diagnosis, fill = new_submain2)) + 
  geom_bar(position='fill') + 
  scale_fill_manual(values = c("#56B4E9", "purple", "red")) + 
  labs(y = '') + 
  theme_classic() + 
  RotatedAxis() + 
  scale_x_discrete(limits = c("Control", "COVID", "IPF")) + 
  theme(legend.position = "right",axis.title.y = element_blank())
ggsave(here("EMT/figure/barplot_proportion_fibro_annotate_by_diagnosis.pdf"), width = 7, height = 5)
