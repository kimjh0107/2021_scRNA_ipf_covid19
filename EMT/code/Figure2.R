library(here)
library(Seurat)
library(tidyverse)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(reshape2)

# Dimplot - porportion 적으로 나누는 것이 가능할지 







# COVID-19 volcano 
markers_up <- read_csv(here("EMT/csv/fibro_upregulated_covid_min0.csv"))
markers_down <- read_csv(here("EMT/csv/fibro_downregulated_covid_min0.csv"))


markers_up$Diagnosis <- "COVID"
markers_down$Diagnosis <- "Control"
markers <- rbind(markers_up, markers_down)

markers$p_val_adj <- as.numeric(markers$p_val_adj)
markers$avg_log2FC <- as.numeric(markers$avg_log2FC)
markers

markers <- markers %>% mutate(avg_logFC_new = ifelse(Diagnosis == 'COVID',avg_log2FC,-avg_log2FC))
markers <- markers[c(order(-markers$avg_logFC_new)), ]

markers$color = '#418bbe'
markers$group = 'Other'
markers$label = markers$gene
rownames(markers) <- markers$gene

group_up = markers %>% filter(., Diagnosis == 'COVID' & p_val_adj < 0.05 )
group_down = markers %>% filter(.,Diagnosis == 'Control' & p_val_adj < 0.05 )


tcell_up <- read_csv(here("EMT/csv/fibro_upregulated_covid_min0.csv"))
tcell_up <- tcell_up %>% filter(p_val_adj < 0.05)



markers[intersect(group_up$gene,tcell_up$gene),'color'] = '#d53d49'




merge_up <- read_csv(here("EMT/csv/fibro_upregulated_merge_min0.csv"))
merge_up <- merge_up[c(order(-merge_up$avg_log2FC.x)),]
merge_up <- head(merge_up$gene, 20)

merge_down <- read_csv(here("EMT/csv/fibroell_downregulated_merge_min0.csv"))
merge_down <- merge_down[c(order(-merge_down$avg_log2FC.x)),]
merge_down <- head(merge_down$gene, 20)

merge_geneset <- c(merge_up, merge_down)
my_gene_list = merge_geneset



markers$p_val_adj = -log10(markers$p_val_adj)
markers <- markers[c(order(-markers$avg_logFC_new)),]

markers[!markers$gene %in% c(my_gene_list),'label'] = ''
markers

my_not_gene_list = c()
markers$Diagnosis = factor(markers$Diagnosis,labels = c('COVID','Control'),levels =c('COVID','Control'))
my_not_gene_list = c()
markers[!markers$gene %in% c(my_gene_list),'label'] = ''
markers = with(markers, markers[order(Diagnosis),])
markers



p1 <- ggplot(markers, aes(avg_logFC_new, p_val_adj, color = Diagnosis, order = Diagnosis)) + geom_point(size = 0.8) + 
  labs(x = "logFC", y = "-Log10(p.adjust)") + theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.text = element_text(size = 5),axis.title = element_text(size = 5),
        axis.title.x = element_text(margin = unit(c(4, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 4, 0, 0), "mm")),
        axis.line = element_line(colour = 'black',size = 0.4),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_line(colour = 'grey',size=0.4,linetype = 2),
        axis.ticks = element_line(colour = 'black',size = 0.4),
        legend.title = element_blank(),legend.text = element_text(size = 5,family = 'sans',face='plain'),
        legend.position = 'top',legend.key.size = unit(1.5, 'lines')) + theme_classic() + 
  scale_color_manual(values = setNames(markers$color, markers$Diagnosis)) #+ 
#geom_text_repel(aes(label = markers$label), verbose = T, seed = 123, max.time = 1, max.iter = Inf, size = 3) + NoLegend()
p1 <- p1 + geom_text_repel(aes(label = markers$label), max.overlaps = 200, , color = "black", size = 2, box.padding = 0.25, label.padding = 0.25,point.padding = 1e-06, segment.color = "grey80", verbose = T, max.iter = Inf, max.time = 1) + NoLegend() + labs(title = "COVID-19 DGE")
p1
ggsave(here("EMT/figure/Figure2/Volcano_fibro_covid.pdf"), width = 7, height = 5)








# IPF volcano 
markers_up <- read_csv(here("EMT/csv/fibro_upregulated_ipf_min0.csv"))
markers_down <- read_csv(here("EMT/csv/fibro_downregulated_ipf_min0.csv"))


markers_up$Diagnosis <- "IPF"
markers_down$Diagnosis <- "Control"
markers <- rbind(markers_up, markers_down)

markers$p_val_adj <- as.numeric(markers$p_val_adj)
markers$avg_log2FC <- as.numeric(markers$avg_log2FC)
markers

markers <- markers %>% mutate(avg_logFC_new = ifelse(Diagnosis == 'IPF',avg_log2FC,-avg_log2FC))
markers <- markers[c(order(-markers$avg_logFC_new)), ]

markers$color = '#418bbe'
markers$group = 'Other'
markers$label = markers$gene
rownames(markers) <- markers$gene

group_up = markers %>% filter(., Diagnosis == 'IPF' & p_val_adj < 0.05 )
group_down = markers %>% filter(.,Diagnosis == 'Control' & p_val_adj < 0.05 )


tcell_up <- read_csv(here("EMT/csv/fibro_upregulated_ipf_min0.csv"))
tcell_up <- tcell_up %>% filter(p_val_adj < 0.05)
tcell_down <- read_csv(here("EMT/csv/fibro_downregulated_ipf_min0.csv"))
tcell_down <- tcell_down %>% filter(p_val_adj < 0.05)

markers[intersect(group_up$gene,tcell_up$gene),'color'] = '#d53d49'

merge_up <- read_csv(here("EMT/csv/fibro_upregulated_merge_min0.csv"))
merge_up <- merge_up[c(order(-merge_up$avg_log2FC.x)),]
merge_up <- head(merge_up$gene, 20)

merge_down <- read_csv(here("EMT/csv/fibroell_downregulated_merge_min0.csv"))
merge_down <- merge_down[c(order(-merge_down$avg_log2FC.x)),]
merge_down <- head(merge_down$gene, 20)

merge_geneset <- c(merge_up, merge_down)
my_gene_list = merge_geneset



markers$p_val_adj = -log10(markers$p_val_adj)
markers <- markers[c(order(-markers$avg_logFC_new)),]

markers[!markers$gene %in% c(my_gene_list),'label'] = ''
markers

my_not_gene_list = c()
markers$Diagnosis = factor(markers$Diagnosis,labels = c('IPF','Control'),levels =c('IPF','Control'))
my_not_gene_list = c()
markers[!markers$gene %in% c(my_gene_list),'label'] = ''
markers = with(markers, markers[order(Diagnosis),])
markers



p2 <- ggplot(markers, aes(avg_logFC_new, p_val_adj, color = Diagnosis, order = Diagnosis)) + geom_point(size = 0.8) + 
  labs(x = "logFC", y = "-Log10(p.adjust)") + theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.text = element_text(size = 5),axis.title = element_text(size = 5),
        axis.title.x = element_text(margin = unit(c(4, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 4, 0, 0), "mm")),
        axis.line = element_line(colour = 'black',size = 0.4),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_line(colour = 'grey',size=0.4,linetype = 2),
        axis.ticks = element_line(colour = 'black',size = 0.4),
        legend.title = element_blank(),legend.text = element_text(size = 5,family = 'sans',face='plain'),
        legend.position = 'top',legend.key.size = unit(1.5, 'lines')) + theme_classic() + 
  scale_color_manual(values = setNames(markers$color, markers$Diagnosis)) #+ 
#geom_text_repel(aes(label = markers$label), verbose = T, seed = 123, max.time = 1, max.iter = Inf, size = 3) + NoLegend()
p2 <- p2 + geom_text_repel(aes(label = markers$label), max.overlaps = 200, , color = "black", size = 2, box.padding = 0.25, label.padding = 0.25,point.padding = 1e-06, segment.color = "grey80", verbose = T, max.iter = Inf, max.time = 1) + NoLegend() + labs(title = "IPF DGE")
p2
ggsave(here("EMT/figure/Figure2/Volcano_fibro_ipf.pdf"), width = 7, height = 5)






#### DGE Heatmap #### 
fibro <- readRDS(here("EMT/RDS/fibro.RDS"))
DefaultAssay(fibro) <- "RNA"
Idents(fibro) <- "Diagnosis"

fibro <- FindVariableFeatures(fibro, selection.method = "vst", nfeatures = 2000)
fibro <- ScaleData(fibro)

# load csv data - avg_log2FC > 0.5
fibro_dge_up <- read_csv(here("EMT/csv/fibro_upregulated_merge_min0.csv"))
fibro_dge_up = fibro_dge_up %>% filter(p_val_adj.x < 0.05 & p_val_adj.y < 0.05 & avg_log2FC.x > 0.5 & avg_log2FC.y > 0.5)

fibro_dge_down <- read_csv(here("EMT/csv/fibroell_downregulated_merge_min0.csv"))
fibro_dge_down = fibro_dge_down %>% filter(p_val_adj.x < 0.05 & p_val_adj.y < 0.05 & avg_log2FC.x > 0.5 & avg_log2FC.y > 0.5)

# fibro
DoHeatmap(fibro, features = fibro_dge_up$gene)
ggsave(here("EMT/figure/Figure2/Heatmap_fibro_up_filter_05.pdf"), width = 20, height = 20)
DoHeatmap(fibro, features = fibro_dge_down$gene)
ggsave(here("EMT/figure/Figure2/Heatmap_fibro_down_filter_05.pdf"), width = 20, height = 20)


# load csv data - avg_log2FC > 0.25
fibro_dge_up <- read_csv(here("EMT/csv/fibro_upregulated_merge_min0.csv"))
fibro_dge_up = fibro_dge_up %>% filter(p_val_adj.x < 0.05 & p_val_adj.y < 0.05 & avg_log2FC.x > 0.25 & avg_log2FC.y > 0.25)

fibro_dge_down <- read_csv(here("EMT/csv/fibroell_downregulated_merge_min0.csv"))
fibro_dge_down = fibro_dge_down %>% filter(p_val_adj.x < 0.05 & p_val_adj.y < 0.05 & avg_log2FC.x > 0.25 & avg_log2FC.y > 0.25)

# fibro
DoHeatmap(fibro, features = fibro_dge_up$gene)
ggsave(here("EMT/figure/Figure2/Heatmap_fibro_up_filter_025.pdf"), width = 20, height = 20)
DoHeatmap(fibro, features = fibro_dge_down$gene)
ggsave(here("EMT/figure/Figure2/Heatmap_fibro_down_filter_025.pdf"), width = 20, height = 20)