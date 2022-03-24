library(Seurat)
library(tidyverse)
library(here)
library(ggpubr)

df <- readRDS(here("final/RDS/df.RDS"))
DefaultAssay(df) <- "RNA"

consistentcolors = colors <- c("#006E82", "#AA0A3C", "#8214A0", "#00A0FA", "#FA5078", 
                               "#005AC8", "#CC79A7", "#FAE6BE", "#0072B2", "#A0FA82", "#F0F032", "#0AB45A", 
                               "#FA7850", "#14D2DC", "#FA78FA")

# calculate frequencies of cell_type_main classes in each sample
df_tobesummed = data.frame(orig.ident = df$Diagnosis_tag, 
                           group = df$Diagnosis, 
                           cell_type_main = df$cell_type_main, 
                           cell_type_submain = df$new_submain)

df_summed = df_tobesummed %>% group_by(orig.ident, cell_type_main, cell_type_submain, group) %>% tally()
df_summed = df_summed %>% group_by(orig.ident) %>% mutate(freq = n/sum(n))



summed_1 <- df_summed %>% filter(cell_type_submain %in% c("AT2"))
summed_2 <- df_summed %>% filter(cell_type_submain %in% c("Monocytes"))
summed_3 <- df_summed %>% filter(cell_type_submain %in% c("Endothelial cells"))
summed_4 <- df_summed %>% filter(cell_type_submain %in% c("Macrophages"))
summed_5 <- df_summed %>% filter(cell_type_submain %in% c("NK cells"))
summed_6 <- df_summed %>% filter(cell_type_submain %in% c("Ciliated cells"))
summed_7 <- df_summed %>% filter(cell_type_submain %in% c("Club cells"))
summed_8 <- df_summed %>% filter(cell_type_submain %in% c("mDCs"))
summed_9 <- df_summed %>% filter(cell_type_submain %in% c("B cells"))
summed_10 <- df_summed %>% filter(cell_type_submain %in% c("AT1"))
summed_11 <- df_summed %>% filter(cell_type_submain %in% c("Fibroblasts"))
summed_12 <- df_summed %>% filter(cell_type_submain %in% c("Basal cells"))
summed_13 <- df_summed %>% filter(cell_type_submain %in% c("Plasma cells"))
summed_14 <- df_summed %>% filter(cell_type_submain %in% c("Gobelt cells"))
summed_15 <- df_summed %>% filter(cell_type_submain %in% c("CD4 Tcells"))
summed_16 <- df_summed %>% filter(cell_type_submain %in% c("Mast cells"))
summed_17 <- df_summed %>% filter(cell_type_submain %in% c("Smooth Muscle cells"))
summed_18 <- df_summed %>% filter(cell_type_submain %in% c("CD8 Tcells"))
summed_19 <- df_summed %>% filter(cell_type_submain %in% c("pDCs"))
summed_20 <- df_summed %>% filter(cell_type_submain %in% c("Neutrophils"))


my_comparisons <- list(c("Control", "COVID"), c("COVID", "IPF"), c("Control", "IPF"))
# function for plot 
get_plot <- function(data){
   ggboxplot(data, x = "group", y = "freq", color = "group", add = "jitter") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    scale_colour_manual(values = consistentcolors[5:2]) +
    stat_compare_means(comparisons = my_comparisons, label = "p.signif") + 
    theme(legend.position = "none") + 
    scale_x_discrete(limits = c("Control", "COVID", "IPF"))
}



plot <- get_plot(summed_1) + labs(title="AT2", x =" ", y = " ") 
plot2 <- get_plot(summed_2) + labs(title="Monocytes", x =" ", y = " ") 
plot3 <- get_plot(summed_3) + labs(title="Endothelial cells", x =" ", y = "Percentage") 
plot4 <- get_plot(summed_4) + labs(title="Macrophages", x =" ", y = " ") 
plot5 <- get_plot(summed_5) + labs(title="NK cells", x =" ", y = " ") 
plot6 <- get_plot(summed_6) + labs(title="Ciliated cells", x =" ", y = " ") 
plot7 <- get_plot(summed_7) + labs(title="Club cells", x =" ", y = " ") 
plot8 <- get_plot(summed_8) + labs(title="mDCs", x =" ", y = "Percentage") 
plot9 <- get_plot(summed_9) + labs(title="B cells", x =" ", y = "") 
plot10 <- get_plot(summed_10) + labs(title="AT1", x =" ", y = "") 
plot11 <- get_plot(summed_11) + labs(title="Fibroblasts", x =" ", y = "") 
plot12 <- get_plot(summed_12) + labs(title="Basal cells", x =" ", y = "") 
plot13 <- get_plot(summed_13) + labs(title="Plasma cells", x =" ", y = "") 
plot14 <- get_plot(summed_14) + labs(title="Gobelt cells", x =" ", y = "") 
plot15 <- get_plot(summed_15) + labs(title="CD4 Tcells", x =" ", y = "") 
plot16 <- get_plot(summed_16) + labs(title="Mast cells", x =" ", y = "") 
plot17 <- get_plot(summed_17) + labs(title="Smooth Muscle cells", x =" ", y = "") 
plot18 <- get_plot(summed_18) + labs(title="CD8 Tcells", x =" ", y = "") 
plot19 <- get_plot(summed_19) + labs(title="pDCs", x =" ", y = "") 
plot20 <- get_plot(summed_20) + labs(title="Neutrophils", x =" ", y = "") 

plot_arrange <- ggarrange(plot10, plot, plot9, plot12, plot15,
                          plot18, plot6, plot7, plot3, plot11, 
                          plot14, plot4, plot16, plot8, plot2, 
                          plot20, plot5, plot19, plot13, plot17, ncol = 5, nrow = 4)
ggsave(plot = plot_arrange, here("final/figure/Boxplot_total_celltype_type2.pdf"), height = 18, width = 10)
    
    
   