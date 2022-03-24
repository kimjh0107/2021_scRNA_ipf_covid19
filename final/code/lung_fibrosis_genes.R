library(Seurat)
library(tidyverse)
library(here)
library(ggpubr)

genes <- read_csv(here('final/csv/lung_fibrosis_genesets.txt'))
gene_list <- genes %>% pull()

df <- readRDS(here("final/RDS/df.RDS"))
DefaultAssay(df) <- "RNA"

my_comparisons <- list(c("Control", "COVID"), c("COVID", "IPF"), c("Control", "IPF"))

exp_data <- FetchData(object = df, vars = c(gene_list,"Diagnosis","Diagnosis_tag", "new_submain"))

gene_data <- exp_data %>% 
    pivot_longer(ATP11A:TNF, names_to='genes', values_to = 'value') %>% 
    group_by(Diagnosis, Diagnosis_tag, new_submain) %>%
    summarize(exp_avg = mean(value)) 

plot <- gene_data %>%
        ggplot(aes(x = Diagnosis, y = exp_avg), color = "Diagnosis") + 
        geom_boxplot() + 
        geom_jitter() + 
        facet_wrap(new_submain ~. ) + 
        theme_bw() + 
        labs(title = 'Fibrosis gene expressions') + 
        stat_compare_means(comparisons = my_comparisons, label = "p.signif") 
      #  scale_x_discrete(limits = c("Control", "COVID", "IPF")
ggsave(here("final/figure/Boxplot_lung_fibrosis_expression.pdf"), height = 18, width = 10)






consistentcolors = colors <- c("#006E82", "#AA0A3C", "#8214A0", "#00A0FA", "#FA5078", 
                               "#005AC8", "#CC79A7", "#FAE6BE", "#0072B2", "#A0FA82", "#F0F032", "#0AB45A", 
                               "#FA7850", "#14D2DC", "#FA78FA")


get_plot <- function(data){
   ggboxplot(data, x = "Diagnosis", y = "exp_avg", color = "Diagnosis", add = "jitter") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    facet_wrap(new_submain ~. ) + 
    theme_bw() + 
    labs(title = 'Fibrosis gene expressions') +
    scale_colour_manual(values = consistentcolors[5:2]) +
    stat_compare_means(comparisons = my_comparisonsggsave(plot = plot_arrange, here("final/figure/Boxplot_total_celltype_type2.pdf"), height = 18, width = 10)
≈) + 
    theme(legend.position = "none") + 
    scale_x_discrete(limits = c("Control", "COVID", "IPF"))
}

plot <- get_plot(gene_data)
plot
ggsave(here("final/figure/Boxplot_lung_fibrosis_expression_ver2.pdf"), height = 18, width = 10)
