library(here)
library(Seurat)
library(tidyverse)

df <- readRDS(here("EMT/RDS/EMT.RDS"))
cells <- df@meta.data 

cell_counts <- cells %>% 
  dplyr::select(Diagnosis, Diagnosis_tag ,new_submain) %>%
  group_by(Diagnosis, Diagnosis_tag ,new_submain) %>%
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
ggplot(df@meta.data, aes(Diagnosis, fill = new_submain)) + 
  geom_bar(position='fill') + 
  scale_fill_manual(values = c("#56B4E9", "purple", "red")) + 
  labs(y = '') + 
  theme_classic() + 
  RotatedAxis() + 
  scale_x_discrete(limits = c("Control", "COVID", "IPF")) + 
  theme(legend.position = "right",axis.title.y = element_blank())
ggsave(here("EMT/figure/barplot_proportion_by_diagnosis.pdf"), width = 7, height = 5)


# proportion table 
sc_control <- subset(df, subset = Diagnosis == "Control")
sc_covid <- subset(df, subset = Diagnosis == "COVID")
sc_ipf <- subset(df, subset = Diagnosis == "IPF")


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


control_proportion <- get_control_proportion(sc_control)
covid_proportion <- get_covid_proportion(sc_covid)
ipf_proportion <- get_ipf_proportion(sc_ipf)
control_proportion
merge <- left_join(control_proportion, covid_proportion, by = "new_submain")
merge <- left_join(merge, ipf_proportion, by = "new_submain")
write_csv(merge, here("EMT/csv/Proportion_EMT_table_by_diagnosis.csv"))
