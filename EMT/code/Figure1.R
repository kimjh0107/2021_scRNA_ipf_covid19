library(here)
library(Seurat)
library(tidyverse)
library(ggpubr)


# Fig1B 
df <- readRDS(here("EMT/RDS/df.RDS"))
DefaultAssay(df) <- "RNA"

DimPlot(df, group.by = "cell_type_submain", label = T, repel = T) + NoLegend()
ggsave(here("EMT/figure/Figure1/Fig1b_DimPlot_celltype.pdf"), width = 7, height = 5)

DimPlot(df, group.by = "Diagnosis", label = F, repel = T, shuffle = T, cols = c("#00A0FA", "red", "purple")) + ggplot2::theme(legend.position = "bottom")
ggsave(here("EMT/figure/Figure1/Fig1b_DimPlot_diagnosis.pdf"), width = 7, height = 5)


# Fig1C
# Boxplot - lung fibrosis gene expression level 
genes <- read_csv(here('EMT/csv/lung_fibrosis_genesets.txt'))
genes
gene_list <- genes %>% pull()
my_comparisons <- list(c("Control", "COVID"), c("COVID", "IPF"), c("Control", "IPF"))
exp_data <- FetchData(object = df, vars = c(gene_list,"Diagnosis","Diagnosis_tag", "new_submain"))

gene_data <- exp_data %>% 
  pivot_longer(ATP11A:TNF, names_to='genes', values_to = 'value') %>% 
  group_by(Diagnosis, Diagnosis_tag, new_submain) %>%
  summarize(exp_avg = mean(value), .groups = "Diagnosis") %>%
  summarize(exp_avg = mean(value)) 
gene_data


consistentcolors = colors <- c("#006E82", "#AA0A3C", "#8214A0", "#FA5078", "#00A0FA", 
                               "#005AC8", "#CC79A7", "#FAE6BE", "#0072B2", "#A0FA82", "#F0F032", "#0AB45A", 
                               "#FA7850", "#14D2DC", "#FA78FA")

get_plot <- function(data){
  ggboxplot(data, x = "Diagnosis", y = "exp_avg", color = "Diagnosis", add = "jitter") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    facet_wrap(new_submain ~. ) + 
    theme_bw() + 
    labs(title = 'Fibrosis gene expressions') +
    scale_colour_manual(values = consistentcolors[5:2]) +
    stat_compare_means(comparisons = my_comparisons, label = "p.signif") + 
    theme(legend.position = "none") + 
    scale_x_discrete(limits = c("Control", "COVID", "IPF"))
}

plot <- get_plot(gene_data)
plot
ggsave(here("EMT/figure/Figure1/Fig1c_Boxplot_lung_fibrosis_expression.pdf"), height = 18, width = 10)

