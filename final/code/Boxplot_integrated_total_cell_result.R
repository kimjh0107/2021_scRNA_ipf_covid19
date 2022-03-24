library(Seurat)
library(tidyverse)
library(here)
library(ggpubr)


df <- readRDS(here("snRNA/RDS/df.RDS"))

consistentcolors = colors <- c("#006E82", "#AA0A3C", "#8214A0", "#00A0FA", "#FA5078", 
                               "#005AC8", "#CC79A7", "#FAE6BE", "#0072B2", "#A0FA82", "#F0F032", "#0AB45A", 
                               "#FA7850", "#14D2DC", "#FA78FA")

# calculate frequencies of cell_type_main classes in each sample
df_tobesummed = data.frame(orig.ident = df$Diagnosis_tag, 
                           group = df$Diagnosis, 
                           cell_type_main = df$cell_type_main, 
                           cell_type_submain = df$cell_type_submain)

df_summed = df_tobesummed %>% group_by(orig.ident, cell_type_main, cell_type_submain, group) %>% tally()
df_summed = df_summed %>% group_by(orig.ident) %>% mutate(freq = n/sum(n))



summed_1 <- df_summed %>% filter(cell_type_submain %in% c("Ciliated cells"))
summed_2 <- df_summed %>% filter(cell_type_submain %in% c("Macrophages"))
summed_3 <- df_summed %>% filter(cell_type_submain %in% c("AT2"))

summed_4 <- df_summed %>% filter(cell_type_submain %in% c("Smooth muscle cells"))
summed_5 <- df_summed %>% filter(cell_type_submain %in% c("AT1"))
summed_6 <- df_summed %>% filter(cell_type_submain %in% c("Fibroblasts"))

summed_7 <- df_summed %>% filter(cell_type_submain %in% c("Epithelial cells"))
summed_8 <- df_summed %>% filter(cell_type_submain %in% c("Dendritic cells"))
summed_9 <- df_summed %>% filter(cell_type_submain %in% c("Immune cells"))

summed_10 <- df_summed %>% filter(cell_type_submain %in% c("Endothelial cells"))
summed_11 <- df_summed %>% filter(cell_type_submain %in% c("Mast cells"))
summed_12 <- df_summed %>% filter(cell_type_submain %in% c("Cycling NKT cells"))

summed_13 <- df_summed %>% filter(cell_type_submain %in% c("Neuronal cells"))
summed_14 <- df_summed %>% filter(cell_type_submain %in% c("Plasma cells"))
summed_15 <- df_summed %>% filter(cell_type_submain %in% c("T cells"))

summed_16 <- df_summed %>% filter(cell_type_submain %in% c("NK cells"))
summed_17 <- df_summed %>% filter(cell_type_submain %in% c("B cells"))
summed_18 <- df_summed %>% filter(cell_type_submain %in% c("Monocytes"))
summed_19 <- df_summed %>% filter(cell_type_submain %in% c("NKT cells"))
summed_20 <- df_summed %>% filter(cell_type_submain %in% c("Neutrophils"))
summed_21 <- df_summed %>% filter(cell_type_submain %in% c("Club cells"))


my_comparisons <- list(c("Control", "COVID-19"), c("COVID-19", "IPF"), c("Control", "IPF"))
# function for plot 
get_plot <- function(data){
   ggboxplot(data, x = "group", y = "freq", color = "group", add = "jitter") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    scale_colour_manual(values = consistentcolors[5:2]) +
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
    theme(legend.position = "none") + 
    scale_x_discrete(limits = c("Control", "COVID-19", "IPF"))
}


plot <- get_plot(summed_1) + labs(title="Ciliated cells", x =" ", y = " ") 
plot2 <- get_plot(summed_2) + labs(title="Macrophages", x =" ", y = " ") 
plot3 <- get_plot(summed_3) + labs(title="AT2", x =" ", y = "Percentage") 
plot4 <- get_plot(summed_4) + labs(title="Smooth muscle cells", x =" ", y = " ") 
plot5 <- get_plot(summed_5) + labs(title="AT1", x =" ", y = " ") 
plot6 <- get_plot(summed_6) + labs(title="Fibroblasts", x =" ", y = " ") 
plot7 <- get_plot(summed_7) + labs(title="Epithelail cells", x =" ", y = " ") 
plot8 <- get_plot(summed_8) + labs(title="Dendritic cells", x =" ", y = "Percentage") 
plot9 <- get_plot(summed_9) + labs(title="Immune cells", x =" ", y = "") 
plot10 <- get_plot(summed_10) + labs(title="Endothelial cells", x =" ", y = "") 
plot11 <- get_plot(summed_11) + labs(title="Mast cells", x =" ", y = "") 
plot12 <- get_plot(summed_12) + labs(title="Cycling NKT cells", x =" ", y = "") 
plot13 <- get_plot(summed_13) + labs(title="Neuronal cells", x =" ", y = "") 
plot14 <- get_plot(summed_14) + labs(title="Plasma cells", x =" ", y = "") 
plot15 <- get_plot(summed_15) + labs(title="T cells", x =" ", y = "") 
plot16 <- get_plot(summed_16) + labs(title="NK cells", x =" ", y = "") 
plot17 <- get_plot(summed_17) + labs(title="B cells", x =" ", y = "") 
plot18 <- get_plot(summed_18) + labs(title="Monocytes", x =" ", y = "") 
plot19 <- get_plot(summed_19) + labs(title="NKT cells", x =" ", y = "") 
plot20 <- get_plot(summed_20) + labs(title="Neutrophils", x =" ", y = "") 



plot_arrange <- ggarrange(plot, plot2, plot3, plot4, plot5,
                          plot6, plot7, plot8, plot9,plot10,
                          plot11, plot12, plot13, plot14, plot15,
                          plot16, plot17, plot18, plot19, plot20
                          ncol = 5, nrow = 4)
#at1_plot
ggsave(plot = plot_arrange, here("final/figure/Boxplot_integrated_total_type2.pdf"), height = 18, width = 10)








