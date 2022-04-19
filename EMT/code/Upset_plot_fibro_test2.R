library(Seurat)
library(dplyr)
library(tidyr)
library(viridis)
library(matrixStats)
library(UpSetR)


fibro <- readRDS(here("EMT/RDS/fibro.RDS"))
DefaultAssay(fibro) <- "RNA"
Idents(fibro) <- "Diagnosis_tag2"

hc1_markers <- FindMarkers(fibro, ident.1 = "HC1")
hc2_markers <- FindMarkers(fibro, ident.1 = "HC2")
hc1_markers$gene <- rownames(hc1_markers)
hc2_markers$gene <- rownames(hc2_markers)


gene_list <- data.frame(Gene = unique(c(hc1_markers$gene, hc2_markers$gene)))

gene_list$hc1 <- gene_list$Gene %in% hc1_markers$gene
gene_list$hc2 <- gene_list$Gene %in% hc2_markers$gene
gene_list[,2:ncol(gene_list)] <- lapply(gene_list[,2:ncol(gene_list)], as.numeric)
gene_list
upset(gene_list, nsets=2, nintersects = 4, keep.order=T,
      sets=c("hc1", "hc2"),
      #intersections=intersection_list,
      order.by = "freq", point.size=2.5,
      mainbar.y.label = "Intersection Size",
      sets.x.label = "Differentially Expressed\nGene Count",
      text.scale=c(1.25, 1.25),
      mb.ratio=c(0.5,0.5))
ggsave(here("EMT/figure/upset_test.pdf"))


gene_list[which(rowSums(gene_list[,2:ncol(gene_list)]) == 2),] #Which genes are significant in all datasets

summary_table <- data.frame(Gene = gene_list$Gene,
                            SampleCount = rowSums(gene_list[,2:ncol(gene_list)]))
summary_table <- arrange(summary_table, desc(SampleCount))
summary_table$Index <- 1:nrow(summary_table)
summary_table$SampleCount <- factor(summary_table$SampleCount)

summary_plot <- ggplot(summary_table, aes(x=Index, y=SampleCount)) +
  geom_point() +
  xlab("Gene Index") + ylab("Significance Counts") +
  scale_x_continuous(expand=c(0.01,0)) +
  theme_classic() +
  theme(axis.text=element_text(size=12, color="black"),
        axis.title=element_text(size=14))
ggsave(here("EMT/figure/test.pdf"))


summary_plot
