library(Seurat)
library(here)
library(tidyverse)
library(biomaRt)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DoMultiBarHeatmap)

####################### fibro - merge ####################### 
upregulated_total <- read_csv(here("EMT/csv/Reactome/Reactome_fibro_merge_up.csv"))


figure <- upregulated_total %>% 
  arrange('Entities FDR') %>% 
  head(25) %>% 
  dplyr::select('Pathway name', 'Entities FDR') %>%
  dplyr::mutate(logFDR = - log(`Entities FDR`)) %>%
  ggplot(aes(reorder(`Pathway name`, logFDR), `logFDR`)) + 
  geom_bar(stat='identity', fill = "#CC3333", alpha=0.5) + 
  coord_flip() + 
  #  geom_text(aes(y=1, label = 'Pathway name'))
  labs(y = 'logFDR', title = 'Fibroblasts common upregulated pathway') + 
  theme_minimal() + theme_classic() + 
  theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 8, face = 'bold'))
figure
ggsave(plot = figure,here('EMT/figure/Reactome/Reactome_fibro_merge_upregulated.pdf'), width = 15, height = 10)

# total downregulated macro 
upregulated_total <- read_csv(here("EMT/csv/Reactome/Reactome_fibro_merge_down.csv"))

figure <- upregulated_total %>% 
  arrange('Entities FDR') %>% 
  head(25) %>% 
  dplyr::select('Pathway name', 'Entities FDR') %>%
  dplyr::mutate(logFDR = - log(`Entities FDR`)) %>%
  ggplot(aes(reorder(`Pathway name`, logFDR), `logFDR`)) + 
  geom_bar(stat='identity', fill = "#3359cc", alpha=0.5) + 
  coord_flip() + 
  #   geom_text(aes(y=1, label = 'Pathway name'))
  labs(y = 'logFDR', title = 'Fibroblasts common downregulated pathway') + 
  theme_minimal() + theme_classic() + 
  theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 8, face = 'bold'))
figure
ggsave(plot = figure,here('EMT/figure/Reactome/Reactome_fibro_merge_downregulated.pdf'), width = 15, height = 10)








####################### fibro - covid ####################### 
upregulated_total <- read_csv(here("EMT/csv/Reactome/Reactome_fibro_covid_up.csv"))


figure <- upregulated_total %>% 
  arrange('Entities FDR') %>% 
  head(25) %>% 
  dplyr::select('Pathway name', 'Entities FDR') %>%
  dplyr::mutate(logFDR = - log(`Entities FDR`)) %>%
  ggplot(aes(reorder(`Pathway name`, logFDR), `logFDR`)) + 
  geom_bar(stat='identity', fill = "#CC3333", alpha=0.5) + 
  coord_flip() + 
  #  geom_text(aes(y=1, label = 'Pathway name'))
  labs(y = 'logFDR', title = 'Fibroblasts COVID-19 upregulated pathway') + 
  theme_minimal() +theme_classic() + 
  theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 8, face = 'bold'))
figure
ggsave(plot = figure,here('EMT/figure/Reactome/Reactome_fibro_covid_upregulated.pdf'), width = 15, height = 10)

# total downregulated macro 
upregulated_total <- read_csv(here("EMT/csv/Reactome/Reactome_fibro_covid_down.csv"))

figure <- upregulated_total %>% 
  arrange('Entities FDR') %>% 
  head(25) %>% 
  dplyr::select('Pathway name', 'Entities FDR') %>%
  dplyr::mutate(logFDR = - log(`Entities FDR`)) %>%
  ggplot(aes(reorder(`Pathway name`, logFDR), `logFDR`)) + 
  geom_bar(stat='identity', fill = "#3359cc", alpha=0.5) + 
  coord_flip() + 
  #   geom_text(aes(y=1, label = 'Pathway name'))
  labs(y = 'logFDR', title = 'Fibroblasts COVID-19 downregulated pathway') + 
  theme_minimal() +theme_classic() + 
  theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 8, face = 'bold'))
figure
ggsave(plot = figure,here('EMT/figure/Reactome/Reactome_fibro_covid_downregulated.pdf'), width = 15, height = 10)







####################### fibro - ipf ####################### 
upregulated_total <- read_csv(here("EMT/csv/Reactome/Reactome_fibro_ipf_up.csv"))



figure <- upregulated_total %>% 
  arrange('Entities FDR') %>% 
  head(25) %>% 
  dplyr::select('Pathway name', 'Entities FDR') %>%
  dplyr::mutate(logFDR = - log(`Entities FDR`)) %>%
  ggplot(aes(reorder(`Pathway name`, logFDR), `logFDR`)) + 
  geom_bar(stat='identity', fill = "#CC3333", alpha=0.5) + 
  coord_flip() + 
  #  geom_text(aes(y=1, label = 'Pathway name'))
  labs(y = 'logFDR', title = 'Fibroblasts IPF upregulated pathway') + 
  theme_minimal() +theme_classic() + 
  theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 8, face = 'bold'))
figure
ggsave(plot = figure,here('EMT/figure/Reactome/Reactome_fibro_ipf_upregulated.pdf'), width = 15, height = 10)

# total downregulated macro 
upregulated_total <- read_csv(here("EMT/csv/Reactome/Reactome_fibro_ipf_down.csv"))
upregulated_total <- upregulated_total %>% filter(p_val_adj < 0.05 & avg_log2FC > 0.5 )
write_csv(upregulated_total, here("EMT/csv/Reactome/fibro_ipf_down_merge.csv"))

figure <- upregulated_total %>% 
  arrange('Entities FDR') %>% 
  head(25) %>% 
  dplyr::select('Pathway name', 'Entities FDR') %>%
  dplyr::mutate(logFDR = - log(`Entities FDR`)) %>%
  ggplot(aes(reorder(`Pathway name`, logFDR), `logFDR`)) + 
  geom_bar(stat='identity', fill = "#3359cc", alpha=0.5) + 
  coord_flip() + 
  #   geom_text(aes(y=1, label = 'Pathway name'))
  labs(y = 'logFDR', title = 'Fibroblasts IPF downregulated pathway') + 
  theme_minimal() +theme_classic() + 
  theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 8, face = 'bold'))
figure
ggsave(plot = figure,here('EMT/figure/Reactome/Reactome_fibro_ipf_downregulated.pdf'), width = 15, height = 10)











library(Seurat)
library(here)
library(tidyverse)
library(biomaRt)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DoMultiBarHeatmap)

####################### fibro - merge ####################### 
upregulated_total <- read_csv(here("EMT/csv/Reactome/Reactome_at2_merge_up.csv"))


figure <- upregulated_total %>% 
  arrange('Entities FDR') %>% 
  head(25) %>% 
  dplyr::select('Pathway name', 'Entities FDR') %>%
  dplyr::mutate(logFDR = - log(`Entities FDR`)) %>%
  ggplot(aes(reorder(`Pathway name`, logFDR), `logFDR`)) + 
  geom_bar(stat='identity', fill = "#CC3333", alpha=0.5) + 
  coord_flip() + 
  #  geom_text(aes(y=1, label = 'Pathway name'))
  labs(y = 'logFDR', title = 'AT2 common upregulated pathway') + 
  theme_minimal() + theme_classic() + 
  theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 8, face = 'bold'))
figure
ggsave(plot = figure,here('EMT/figure/Reactome/Reactome_at2_merge_upregulated.pdf'), width = 15, height = 10)

# total downregulated macro 
upregulated_total <- read_csv(here("EMT/csv/Reactome/Reactome_at2_merge_down.csv"))

figure <- upregulated_total %>% 
  arrange('Entities FDR') %>% 
  head(25) %>% 
  dplyr::select('Pathway name', 'Entities FDR') %>%
  dplyr::mutate(logFDR = - log(`Entities FDR`)) %>%
  ggplot(aes(reorder(`Pathway name`, logFDR), `logFDR`)) + 
  geom_bar(stat='identity', fill = "#3359cc", alpha=0.5) + 
  coord_flip() + 
  #   geom_text(aes(y=1, label = 'Pathway name'))
  labs(y = 'logFDR', title = 'AT2 common downregulated pathway') + 
  theme_minimal() + theme_classic() + 
  theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 8, face = 'bold'))
figure
ggsave(plot = figure,here('EMT/figure/Reactome/Reactome_at2_merge_downregulated.pdf'), width = 15, height = 10)








####################### at2 - covid ####################### 
upregulated_total <- read_csv(here("EMT/csv/Reactome/Reactome_at2_covid_up.csv"))


figure <- upregulated_total %>% 
  arrange('Entities FDR') %>% 
  head(25) %>% 
  dplyr::select('Pathway name', 'Entities FDR') %>%
  dplyr::mutate(logFDR = - log(`Entities FDR`)) %>%
  ggplot(aes(reorder(`Pathway name`, logFDR), `logFDR`)) + 
  geom_bar(stat='identity', fill = "#CC3333", alpha=0.5) + 
  coord_flip() + 
  #  geom_text(aes(y=1, label = 'Pathway name'))
  labs(y = 'logFDR', title = 'AT2 COVID-19 upregulated pathway') + 
  theme_minimal() +theme_classic() + 
  theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 8, face = 'bold'))
figure
ggsave(plot = figure,here('EMT/figure/Reactome/Reactome_at2_covid_upregulated.pdf'), width = 15, height = 10)

# total downregulated macro 
upregulated_total <- read_csv(here("EMT/csv/Reactome/Reactome_at2_covid_down.csv"))

figure <- upregulated_total %>% 
  arrange('Entities FDR') %>% 
  head(25) %>% 
  dplyr::select('Pathway name', 'Entities FDR') %>%
  dplyr::mutate(logFDR = - log(`Entities FDR`)) %>%
  ggplot(aes(reorder(`Pathway name`, logFDR), `logFDR`)) + 
  geom_bar(stat='identity', fill = "#3359cc", alpha=0.5) + 
  coord_flip() + 
  #   geom_text(aes(y=1, label = 'Pathway name'))
  labs(y = 'logFDR', title = 'AT2 COVID-19 downregulated pathway') + 
  theme_minimal() +theme_classic() + 
  theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 8, face = 'bold'))
figure
ggsave(plot = figure,here('EMT/figure/Reactome/Reactome_at2_covid_downregulated.pdf'), width = 15, height = 10)







####################### at2 - ipf ####################### 
upregulated_total <- read_csv(here("EMT/csv/Reactome/Reactome_at2_ipf_up.csv"))



figure <- upregulated_total %>% 
  arrange('Entities FDR') %>% 
  head(25) %>% 
  dplyr::select('Pathway name', 'Entities FDR') %>%
  dplyr::mutate(logFDR = - log(`Entities FDR`)) %>%
  ggplot(aes(reorder(`Pathway name`, logFDR), `logFDR`)) + 
  geom_bar(stat='identity', fill = "#CC3333", alpha=0.5) + 
  coord_flip() + 
  #  geom_text(aes(y=1, label = 'Pathway name'))
  labs(y = 'logFDR', title = 'AT2 IPF upregulated pathway') + 
  theme_minimal() +theme_classic() + 
  theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 8, face = 'bold'))
figure
ggsave(plot = figure,here('EMT/figure/Reactome/Reactome_at2_ipf_upregulated.pdf'), width = 15, height = 10)

# total downregulated macro 
upregulated_total <- read_csv(here("EMT/csv/Reactome/Reactome_at2_ipf_down.csv"))
upregulated_total <- upregulated_total %>% filter(p_val_adj < 0.05 & avg_log2FC > 0.5 )
write_csv(upregulated_total, here("EMT/csv/Reactome/at2_ipf_down_merge.csv"))

figure <- upregulated_total %>% 
  arrange('Entities FDR') %>% 
  head(25) %>% 
  dplyr::select('Pathway name', 'Entities FDR') %>%
  dplyr::mutate(logFDR = - log(`Entities FDR`)) %>%
  ggplot(aes(reorder(`Pathway name`, logFDR), `logFDR`)) + 
  geom_bar(stat='identity', fill = "#3359cc", alpha=0.5) + 
  coord_flip() + 
  #   geom_text(aes(y=1, label = 'Pathway name'))
  labs(y = 'logFDR', title = 'AT2 IPF downregulated pathway') + 
  theme_minimal() +theme_classic() + 
  theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 8, face = 'bold'))
figure
ggsave(plot = figure,here('EMT/figure/Reactome/Reactome_at2_ipf_downregulated.pdf'), width = 15, height = 10)
