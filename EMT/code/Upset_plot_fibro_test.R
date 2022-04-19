library(Seurat)
library(here)
library(tidyverse)
library(UpSetR)
library(sctransform)

fibro <- readRDS(here("EMT/RDS/fibro.RDS"))
DefaultAssay(fibro) <- "RNA"

fibro <- FindVariableFeatures(fibro, selection.method = "vst", nfeatures = 2000)
fibro <- ScaleData(fibro)
fibro <- RunPCA(object = fibro)
fibro <- FindNeighbors(fibro, dims = 1:30)
fibro <- FindClusters(fibro, resolution = 0.6)
fibro <- RunUMAP(object = fibro, dims = 1:30)


# Subset into each Diagnosis tag 
HC1 <- subset(fibro, subset = Diagnosis_tag2 == "HC1")
HC2 <- subset(fibro, subset = Diagnosis_tag2 == "HC2")
HC3 <- subset(fibro, subset = Diagnosis_tag2 == "HC3")
HC4 <- subset(fibro, subset = Diagnosis_tag2 == "HC4")
HC5 <- subset(fibro, subset = Diagnosis_tag2 == "HC5")
HC6 <- subset(fibro, subset = Diagnosis_tag2 == "HC6")
HC7 <- subset(fibro, subset = Diagnosis_tag2 == "HC7")
HC8 <- subset(fibro, subset = Diagnosis_tag2 == "HC8")
HC9 <- subset(fibro, subset = Diagnosis_tag2 == "HC9")
HC10 <- subset(fibro, subset = Diagnosis_tag2 == "HC10")

C1_1 <- subset(fibro, subset = Diagnosis_tag2 %in% c("C1_1"))
C1_2 <- subset(fibro, subset = Diagnosis_tag2 %in% c("C1_2"))
C1_3 <- subset(fibro, subset = Diagnosis_tag2 %in% c("C1_3"))
C1_4 <- subset(fibro, subset = Diagnosis_tag2 %in% c("C1_4"))

I1 <- subset(fibro, subset = Diagnosis_tag2 == "I1")
I2 <- subset(fibro, subset = Diagnosis_tag2 == "I2")
I3 <- subset(fibro, subset = Diagnosis_tag2 == "I3")
I4 <- subset(fibro, subset = Diagnosis_tag2 == "I4")
I5 <- subset(fibro, subset = Diagnosis_tag2 == "I5")
I6 <- subset(fibro, subset = Diagnosis_tag2 == "I6")
I7 <- subset(fibro, subset = Diagnosis_tag2 == "I7")
I8 <- subset(fibro, subset = Diagnosis_tag2 == "I8")
I9 <- subset(fibro, subset = Diagnosis_tag2 == "I9")
I10 <- subset(fibro, subset = Diagnosis_tag2 == "I10")
I11 <- subset(fibro, subset = Diagnosis_tag2 == "I11")
I12 <- subset(fibro, subset = Diagnosis_tag2 == "I12")




DefaultAssay(HC1) <- "RNA"
DefaultAssay(HC2) <- "RNA"
DefaultAssay(HC3) <- "RNA"
DefaultAssay(HC4) <- "RNA"
DefaultAssay(HC5) <- "RNA"
DefaultAssay(HC6) <- "RNA"
DefaultAssay(HC7) <- "RNA"
DefaultAssay(HC8) <- "RNA"
DefaultAssay(HC9) <- "RNA"
DefaultAssay(HC10) <- "RNA"


DefaultAssay(C1_1) <- "RNA"
DefaultAssay(C1_2) <- "RNA"
DefaultAssay(C1_3) <- "RNA"
DefaultAssay(C1_4) <- "RNA"


DefaultAssay(I1) <- "RNA"
DefaultAssay(I2) <- "RNA"
DefaultAssay(I3) <- "RNA"
DefaultAssay(I4) <- "RNA"
DefaultAssay(I5) <- "RNA"
DefaultAssay(I6) <- "RNA"
DefaultAssay(I7) <- "RNA"
DefaultAssay(I8) <- "RNA"
DefaultAssay(I9) <- "RNA"
DefaultAssay(I10) <- "RNA"
DefaultAssay(I11) <- "RNA"
DefaultAssay(I12) <- "RNA"


getHVG_data <- function(seurat_object, top_n){
  gene_list <- seurat_object[["RNA"]]@meta.features
  gene_list$Gene <- rownames(gene_list)
  gene_list <- arrange(gene_list, desc(vst.variance.standardized))[1:top_n,]
}



gene_count <- 1000

hvg_control_1 <- getHVG_data(HC1, gene_count)
hvg_control_2 <- getHVG_data(HC2, gene_count)
hvg_control_3 <- getHVG_data(HC3, gene_count)
hvg_control_4 <- getHVG_data(HC4, gene_count)
hvg_control_5 <- getHVG_data(HC5, gene_count)
hvg_control_6 <- getHVG_data(HC6, gene_count)
hvg_control_7 <- getHVG_data(HC7, gene_count)
hvg_control_8 <- getHVG_data(HC8, gene_count)
hvg_control_9 <- getHVG_data(HC9, gene_count)
hvg_control_10 <- getHVG_data(HC10, gene_count)


hvg_covid_1 <- getHVG_data(C1_1, gene_count)
hvg_covid_2 <- getHVG_data(C1_2, gene_count)
hvg_covid_3 <- getHVG_data(C1_3, gene_count)
hvg_covid_4 <- getHVG_data(C1_4, gene_count)


hvg_ipf_1 <- getHVG_data(I1, gene_count)
hvg_ipf_2 <- getHVG_data(I2, gene_count)
hvg_ipf_3 <- getHVG_data(I3, gene_count)
hvg_ipf_4 <- getHVG_data(I4, gene_count)
hvg_ipf_5 <- getHVG_data(I5, gene_count)
hvg_ipf_6 <- getHVG_data(I6, gene_count)
hvg_ipf_7 <- getHVG_data(I7, gene_count)
hvg_ipf_8 <- getHVG_data(I8, gene_count)
hvg_ipf_9 <- getHVG_data(I9, gene_count)
hvg_ipf_10 <- getHVG_data(I10, gene_count)
hvg_ipf_11 <- getHVG_data(I11, gene_count)
hvg_ipf_12 <- getHVG_data(I12, gene_count)


hvg_list <- data.frame(Gene = unique(c(hvg_control_1$Gene, 
                                       hvg_control_2$Gene,
                                       hvg_control_3$Gene,
                                       hvg_control_4$Gene,
                                       hvg_control_5$Gene,
                                       hvg_control_6$Gene,
                                       hvg_control_7$Gene,
                                       hvg_control_8$Gene,
                                       hvg_control_9$Gene,
                                       hvg_control_10$Gene,
               
                                       hvg_covid_1$Gene,
                                       hvg_covid_2$Gene,
                                       hvg_covid_3$Gene,
                                       hvg_covid_4$Gene,
                                       
                                       hvg_ipf_1$Gene,
                                       hvg_ipf_2$Gene,
                                       hvg_ipf_3$Gene,
                                       hvg_ipf_4$Gene,
                                       hvg_ipf_5$Gene,
                                       hvg_ipf_6$Gene,
                                       hvg_ipf_7$Gene,
                                       hvg_ipf_8$Gene,
                                       hvg_ipf_9$Gene,
                                       hvg_ipf_10$Gene,
                                       hvg_ipf_11$Gene,
                                       hvg_ipf_12$Gene
                                       )))



hvg_list$control_1 <- hvg_list$Gene %in% hvg_control_1$Gene
hvg_list$control_2 <- hvg_list$Gene %in% hvg_control_2$Gene
hvg_list$control_3 <- hvg_list$Gene %in% hvg_control_3$Gene
hvg_list$control_4 <- hvg_list$Gene %in% hvg_control_4$Gene
hvg_list$control_5 <- hvg_list$Gene %in% hvg_control_5$Gene
hvg_list$control_6 <- hvg_list$Gene %in% hvg_control_6$Gene
hvg_list$control_7 <- hvg_list$Gene %in% hvg_control_7$Gene
hvg_list$control_8 <- hvg_list$Gene %in% hvg_control_8$Gene
hvg_list$control_9 <- hvg_list$Gene %in% hvg_control_9$Gene
hvg_list$control_10 <- hvg_list$Gene %in% hvg_control_10$Gene

hvg_list$covid_1 <- hvg_list$Gene %in% hvg_covid_1$Gene
hvg_list$covid_2 <- hvg_list$Gene %in% hvg_covid_2$Gene
hvg_list$covid_3 <- hvg_list$Gene %in% hvg_covid_3$Gene
hvg_list$covid_4 <- hvg_list$Gene %in% hvg_covid_4$Gene

hvg_list$ipf_1 <- hvg_list$Gene %in% hvg_ipf_1$Gene
hvg_list$ipf_2 <- hvg_list$Gene %in% hvg_ipf_2$Gene
hvg_list$ipf_3 <- hvg_list$Gene %in% hvg_ipf_3$Gene
hvg_list$ipf_4 <- hvg_list$Gene %in% hvg_ipf_4$Gene
hvg_list$ipf_5 <- hvg_list$Gene %in% hvg_ipf_5$Gene
hvg_list$ipf_6 <- hvg_list$Gene %in% hvg_ipf_6$Gene
hvg_list$ipf_7 <- hvg_list$Gene %in% hvg_ipf_7$Gene
hvg_list$ipf_8 <- hvg_list$Gene %in% hvg_ipf_8$Gene
hvg_list$ipf_9 <- hvg_list$Gene %in% hvg_ipf_9$Gene
hvg_list$ipf_10 <- hvg_list$Gene %in% hvg_ipf_10$Gene
hvg_list$ipf_11 <- hvg_list$Gene %in% hvg_ipf_11$Gene
hvg_list$ipf_12 <- hvg_list$Gene %in% hvg_ipf_12$Gene


#Convert TRUE/FALSE to 1/0
hvg_list[,2:ncol(hvg_list)] <- lapply(hvg_list[,2:ncol(hvg_list)], as.numeric)

#Remove cell cycle genes because they are variable across all conditions
hvg_list <- hvg_list[-which(hvg_list$Gene %in% c(cc.genes[[1]], cc.genes[[2]])),]

upset(hvg_list, nsets=26, nintersects=15, keep.order=T,
      sets=c("control_1", "control_2","control_3","control_4","control_5","control_6","control_7","control_8","control_9","control_10",
             
             "covid_1", "covid_2","covid_3","covid_4",
             "ipf_1", "ipf_2", "ipf_3", "ipf_4", "ipf_5","ipf_6","ipf_7","ipf_8","ipf_9","ipf_10",
             "ipf_11", "ipf_12"),
      order.by = "freq", point.size=2.5,
      mainbar.y.label = "Variable Gene\nIntersection Size",
      sets.x.label = "Variable Gene Count",
      text.scale=c(1.25, 1.25, 1.25,1.25, 1.5, 1.25),
      mb.ratio=c(0.5,0.5))




gene_list[which(rowSums(gene_list[,2:ncol(gene_list)]) == 26),] #Which genes are significant in all datasets


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


ggsave(here("EMT/figure/Upset_plot_fibro.pdf"), width = 6, height = 5)




