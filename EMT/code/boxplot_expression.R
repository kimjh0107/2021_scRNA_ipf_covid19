library(here)
library(Seurat)
library(tidyverse)
library(ggpubr)

#### EMT  fibro ####
df <- readRDS(here("EMT/RDS/fibro.RDS"))
DefaultAssay(df) <- "RNA"
# Boxplot - MET gene expression level 
genes <- read_csv(here("EMT/csv/EMT_geneset.txt"))
gene_list <- genes %>% pull()
genes

my_comparisons <- list(c("Control", "COVID"), c("COVID", "IPF"), c("Control", "IPF"))
exp_data <- FetchData(object = df, vars = c(gene_list,"Diagnosis","Diagnosis_tag2", "new_submain"))

gene_data <- exp_data %>% 
  pivot_longer(SERP1:HIF1A, names_to='genes', values_to = 'value') %>% 
  group_by(Diagnosis, Diagnosis_tag2, new_submain) %>%
  summarize(exp_avg = mean(value)) 

plot <- gene_data %>%
  ggplot(aes(x = Diagnosis, y = exp_avg), color = "Diagnosis") + 
  geom_boxplot() + 
  geom_jitter() + 
  facet_wrap(new_submain ~. ) + 
  theme_bw() + 
  labs(title = 'EMT gene expressions') + 
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") +theme_classic()
plot


#### MET  fibro ####
df <- readRDS(here("EMT/RDS/fibro.RDS"))
DefaultAssay(df) <- "RNA"
# Boxplot - MET gene expression level 
genes <- read_csv(here("EMT/csv/MET_genesets2.txt"))
gene_list <- genes %>% pull()

my_comparisons <- list(c("Control", "COVID"), c("COVID", "IPF"), c("Control", "IPF"))
exp_data <- FetchData(object = df, vars = c(gene_list,"Diagnosis","Diagnosis_tag2", "new_submain"))

gene_data <- exp_data %>% 
  pivot_longer(TNC:HGF, names_to='genes', values_to = 'value') %>% 
  group_by(Diagnosis, Diagnosis_tag2, new_submain) %>%
  summarize(exp_avg = mean(value)) 

plot <- gene_data %>%
  ggplot(aes(x = Diagnosis, y = exp_avg), color = "Diagnosis") + 
  geom_boxplot() + 
  geom_jitter() + 
  facet_wrap(new_submain ~. ) + 
  theme_bw() + 
  labs(title = 'MET gene expressions') + 
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") +theme_classic()
plot





#### MET  fibro ####
df <- readRDS(here("EMT/RDS/EMT.RDS"))
DefaultAssay(df) <- "RNA"
# Boxplot - MET gene expression level 
genes <- read_csv(here("EMT/csv/MET_genesets2.txt"))
gene_list <- genes %>% pull()

my_comparisons <- list(c("Control", "COVID"), c("COVID", "IPF"), c("Control", "IPF"))
exp_data <- FetchData(object = df, vars = c(gene_list,"Diagnosis","Diagnosis_tag2", "new_submain"))

gene_data <- exp_data %>% 
  pivot_longer(TNC:HOXA5, names_to='genes', values_to = 'value') %>% 
  group_by(Diagnosis, Diagnosis_tag2, new_submain) %>%
  summarize(exp_avg = mean(value)) 

plot <- gene_data %>%
  ggplot(aes(x = Diagnosis, y = exp_avg), color = "Diagnosis") + 
  geom_boxplot() + 
  geom_jitter() + 
  facet_wrap(new_submain ~. ) + 
  theme_bw() + 
  labs(title = 'MET gene expressions') + 
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") +theme_classic()
plot





# AT 
at <- readRDS(here("EMT/RDS/AT.RDS"))
DefaultAssay(at) <- "RNA"
at1 <- subset(at, subset = new_submain %in% c("AT1"))
at2 <- subset(at, subset = new_submain %in% c("AT2"))


## incomple_transition_marker_gene 
gene_list <- c("ETV5")

my_comparisons <- list(c("Control", "COVID"), c("COVID", "IPF"), c("Control", "IPF"))
exp_data <- FetchData(object = at2, vars = c(gene_list,"Diagnosis","Diagnosis_tag2", "new_submain"))

gene_data <- exp_data %>% 
  pivot_longer(ETV5, names_to='genes', values_to = 'value') %>% 
  group_by(Diagnosis, Diagnosis_tag2, new_submain) %>%
  summarize(exp_avg = mean(value)) 

plot <- gene_data %>%
  ggplot(aes(x = Diagnosis, y = exp_avg), color = "Diagnosis") + 
  geom_boxplot() + 
  geom_jitter() + 
  facet_wrap(new_submain ~. ) + 
  theme_bw() + 
  labs(title = 'Incomple_transition_marker_gene ') + 
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") +theme_classic()
plot




## late_AT1_maturation 
gene_list <- c("CAV1")

my_comparisons <- list(c("Control", "COVID"), c("COVID", "IPF"), c("Control", "IPF"))
exp_data <- FetchData(object = at1, vars = c(gene_list,"Diagnosis","Diagnosis_tag2", "new_submain"))

gene_data <- exp_data %>% 
  pivot_longer(ETV5, names_to='genes', values_to = 'value') %>% 
  group_by(Diagnosis, Diagnosis_tag2, new_submain) %>%
  summarize(exp_avg = mean(value)) 

plot <- gene_data %>%
  ggplot(aes(x = Diagnosis, y = exp_avg), color = c("blue", "red", "purple")) + 
  geom_boxplot() + 
  geom_jitter(aes(color=Diagnosis)) + 
  facet_wrap(new_submain ~. ) + 
  theme_bw() + 
  labs(title = 'late_AT1_maturation ') + 
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") +theme_classic()
plot
