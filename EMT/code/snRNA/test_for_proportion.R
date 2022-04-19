library(here)
library(Seurat)
library(tidyverse)

emt <- readRDS(here("EMT/RDS/EMT.RDS"))

emt@meta.data$new_submain2 <- ifelse(emt@meta.data$seurat_clusters %in% c(10,8,17,6), "AT1", "NA")
emt@meta.data$new_submain2 <- ifelse(emt@meta.data$seurat_clusters %in% c(1,0,11,4,5,2,9,13,16), "AT2", emt@meta.data$new_submain2)
emt@meta.data$new_submain2 <- ifelse(emt@meta.data$seurat_clusters %in% c(15,14,7), "Myofibroblasts", emt@meta.data$new_submain2)
emt@meta.data$new_submain2 <- ifelse(emt@meta.data$seurat_clusters %in% c(3,12), "Fibroblasts", emt@meta.data$new_submain2)
#saveRDS(emt, here("EMT/RDS/emt_new_annotation2.RDS"))

#emt <- readRDS(here("EMT/RDS/emt_new_annotation.RDS"))
subset <- subset(emt, subset = new_submain2 %in% c("Fibroblasts", "Myofibroblasts"))
subset <- ScaleData(object = subset)
subset <- RunPCA(object = subset)
subset <- FindNeighbors(subset, dims = 1:30)
subset <- FindClusters(subset, resolution = 0.6)
subset <- RunUMAP(object = subset, dims = 1:30)
#saveRDS(subset, file = here("EMT/RDS/fibro_annotation.RDS"))


cells <- subset@meta.data 

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
ggplot(subset@meta.data, aes(Diagnosis, fill = new_submain2)) + 
  geom_bar(position='fill') + 
  scale_fill_manual(values = c("#56B4E9", "purple", "red")) + 
  labs(y = '') + 
  theme_classic() + 
  RotatedAxis() + 
  scale_x_discrete(limits = c("Control", "COVID", "IPF")) + 
  theme(legend.position = "right",axis.title.y = element_blank())
ggsave(here("EMT/figure/barplot_proportion_fibro_annotate_by_diagnosis_test.pdf"), width = 7, height = 5)
