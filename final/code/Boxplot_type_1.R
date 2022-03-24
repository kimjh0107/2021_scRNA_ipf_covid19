library(Seurat)
library(tidyverse)
library(here)
library(ggpubr)

df <- readRDS(here("final/RDS/df.RDS"))

# Total 

consistentcolors = colors <- c("#006E82", "#AA0A3C", "#8214A0", "#00A0FA", "#FA5078", 
                               "#005AC8", "#CC79A7", "#FAE6BE", "#0072B2", "#A0FA82", "#F0F032", "#0AB45A", 
                               "#FA7850", "#14D2DC", "#FA78FA")


# calculate frequencies of cell_type_main classes in each sample
df_tobesummed = data.frame(orig.ident = df$Diagnosis_tag, 
                           group = df$Diagnosis, 
                           cell_type_main = df$cell_type_main, 
                           cell_type_submain = df$cell_type_newmain)

df_summed = df_tobesummed %>% group_by(orig.ident, cell_type_main, cell_type_submain, group) %>% tally()
df_summed = df_summed %>% group_by(orig.ident) %>% mutate(freq = n/sum(n))


# make boxplots of cell type frequencies in COVID-19 and Control, using cell type
ggboxplot(df_summed, x = "cell_type_submain", y = "freq", color = "group", add = "jitter") + ylim(0, 0.8) + 
  stat_compare_means(aes(group = group), label = "p.format", method = "wilcox.test", size = 2,  label.y = 0.8) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_colour_manual(values = consistentcolors[1:6])
ggsave(here("final/figure/Boxplot_total_type1.pdf"), width = 11, height = 7)


