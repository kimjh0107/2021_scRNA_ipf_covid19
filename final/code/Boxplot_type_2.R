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

unique(df_summed$cell_type_submain)


at1_summed <- df_summed %>% filter(cell_type_submain %in% c("Endothelia"))
at2_summed <- df_summed %>% filter(cell_type_submain %in% c("AT"))
macro_summed <- df_summed %>% filter(cell_type_submain %in% c("Epithelia"))
nk_summed <- df_summed %>% filter(cell_type_submain %in% c("APC"))
mono_summed <- df_summed %>% filter(cell_type_submain %in% c("B/Plasma"))
fibro_summed <- df_summed %>% filter(cell_type_submain %in% c("Mast"))
tcell_summed <- df_summed %>% filter(cell_type_submain %in% c("Neutrophil"))
nkt_summed <- df_summed %>% filter(cell_type_submain %in% c("T/NK"))
nkt_summed2 <- df_summed %>% filter(cell_type_submain %in% c("Fibroblasts"))
nkt_summed3 <- df_summed %>% filter(cell_type_submain %in% c("Smooth Muscle"))


my_comparisons <- list(c("Control", "COVID"), c("COVID", "IPF"), c("Control", "IPF"))

# function for plot 
#get_plot <- function(data){
#   ggboxplot(data, x = "group", y = "freq", color = "group", add = "jitter") + 
#    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
#    scale_colour_manual(values = consistentcolors[5:2]) +
#    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
#    theme(legend.position = "none") + 
#    scale_x_discrete(limits = c("Control", "COVID", "IPF"))
#}


get_plot <- function(data){
   ggboxplot(data, x = "group", y = "freq", color = "group", add = "jitter") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    scale_colour_manual(values = consistentcolors[5:2]) +
    stat_compare_means(comparisons = my_comparisons, label = "p.signif") + 
    stat_compare_means(label.y=50) +
    theme(legend.position = "none") + 
    scale_x_discrete(limits = c("Control", "COVID", "IPF"))
}


at1_plot <- get_plot(at1_summed) + labs(title="Endothelia", x =" ", y = " ") 
at2_plot <- get_plot(at2_summed) + labs(title="AT", x =" ", y = " ") 
macro_plot <- get_plot(macro_summed) + labs(title="Epithelia", x =" ", y = "Percentage") 
nk_plot <- get_plot(nk_summed) + labs(title="APC", x =" ", y = " ") 
mono_plot <- get_plot(mono_summed) + labs(title="B/Plasma", x =" ", y = " ") 
fibro_plot <- get_plot(fibro_summed) + labs(title="Mast", x =" ", y = " ") 
tcell_plot <- get_plot(tcell_summed) + labs(title="Neutrophil", x =" ", y = " ") 
nkt_plot <- get_plot(nkt_summed) + labs(title="T/NK", x =" ", y = "Percentage") 
nkt_plot2 <- get_plot(nkt_summed2) + labs(title="Fibroblasts", x =" ", y = "") 
nkt_plot3 <- get_plot(nkt_summed3) + labs(title="Smooth Muscle", x =" ", y = "") 




plot_arrange <- ggarrange(macro_plot, mono_plot, tcell_plot, nk_plot, nkt_plot2,
                          nkt_plot, at1_plot, at2_plot, fibro_plot,nkt_plot3,
                          ncol = 5, nrow = 5)
#at1_plot
ggsave(plot = plot_arrange, here("final/figure/Boxplot_total_type2_ver2.pdf"), height = 18, width = 10)






