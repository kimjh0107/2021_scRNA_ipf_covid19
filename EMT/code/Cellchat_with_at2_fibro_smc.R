library(here)
library(Seurat)
library(tidyverse)
library(CellChat)

df <- readRDS(here("EMT/RDS/at2_fibro_smc.RDS"))
DefaultAssay(df) <- "RNA"

                
df@meta.data$tag <- ifelse(df@meta.data$seurat_clusters %in% c(9), "AT2_1", "NA")
df@meta.data$tag <- ifelse(df@meta.data$seurat_clusters %in% c(6), "AT2_2", df@meta.data$tag)
df@meta.data$tag <- ifelse(df@meta.data$seurat_clusters %in% c(4), "AT2_3", df@meta.data$tag)
df@meta.data$tag <- ifelse(df@meta.data$seurat_clusters %in% c(3), "AT2_4", df@meta.data$tag)
df@meta.data$tag <- ifelse(df@meta.data$seurat_clusters %in% c(1), "AT2_5", df@meta.data$tag)
df@meta.data$tag <- ifelse(df@meta.data$seurat_clusters %in% c(0), "AT2_6", df@meta.data$tag)
df@meta.data$tag <- ifelse(df@meta.data$seurat_clusters %in% c(16), "AT2_7", df@meta.data$tag)
df@meta.data$tag <- ifelse(df@meta.data$seurat_clusters %in% c(7), "AT2_8", df@meta.data$tag)
df@meta.data$tag <- ifelse(df@meta.data$seurat_clusters %in% c(12), "AT2_9", df@meta.data$tag)
df@meta.data$tag <- ifelse(df@meta.data$seurat_clusters %in% c(14), "AT2_10", df@meta.data$tag)

df@meta.data$tag <- ifelse(df@meta.data$seurat_clusters %in% c(8), "SMC_1", df@meta.data$tag)
df@meta.data$tag <- ifelse(df@meta.data$seurat_clusters %in% c(11), "SMC_2", df@meta.data$tag)

df@meta.data$tag <- ifelse(df@meta.data$seurat_clusters %in% c(15), "F_1", df@meta.data$tag)
df@meta.data$tag <- ifelse(df@meta.data$seurat_clusters %in% c(13), "F_2", df@meta.data$tag)
df@meta.data$tag <- ifelse(df@meta.data$seurat_clusters %in% c(5), "F_3", df@meta.data$tag)
df@meta.data$tag <- ifelse(df@meta.data$seurat_clusters %in% c(10), "F_4", df@meta.data$tag)
df@meta.data$tag <- ifelse(df@meta.data$seurat_clusters %in% c(2), "F_5", df@meta.data$tag)


control <- subset(df, subset = Diagnosis == "Control")
covid <- subset(df, subset = Diagnosis == "COVID")
ipf <- subset(df, subset = Diagnosis == "IPF")



#### subset by each diagnosis
# control
data_input <- GetAssayData(control, assay = "RNA", slot = "data")
Idents(control) <- "tag"
labels <- Idents(control)
meta <- data.frame(labels = labels, row.names = names(labels))
cellchat_control <- createCellChat(object = data_input, meta = meta, group.by = "labels")

# covid
data_input <- GetAssayData(covid, assay = "RNA", slot = "data")
Idents(covid) <- "tag"
labels <- Idents(covid)
meta <- data.frame(labels = labels, row.names = names(labels))
cellchat_covid <- createCellChat(object = data_input, meta = meta, group.by = "labels")

# ipf
data_input <- GetAssayData(ipf, assay = "RNA", slot = "data")
Idents(ipf) <- "tag"
labels <- Idents(ipf)
meta <- data.frame(labels = labels, row.names = names(labels))
cellchat_ipf <- createCellChat(object = data_input, meta = meta, group.by = "labels")


#### Add meta 
# control
data_input <- GetAssayData(control, assay = "RNA", slot = "data")
Idents(control) <- "tag"
labels <- Idents(control)
meta <- data.frame(labels = labels, row.names = names(labels))
cellchat_control <- createCellChat(object = data_input, meta = meta, group.by = "labels")

# covid
data_input <- GetAssayData(covid, assay = "RNA", slot = "data")
Idents(covid) <- "tag"
labels <- Idents(covid)
meta <- data.frame(labels = labels, row.names = names(labels))
cellchat_covid <- createCellChat(object = data_input, meta = meta, group.by = "labels")

# ipf
data_input <- GetAssayData(ipf, assay = "RNA", slot = "data")
Idents(ipf) <- "tag"
labels <- Idents(ipf)
meta <- data.frame(labels = labels, row.names = names(labels))
cellchat_ipf <- createCellChat(object = data_input, meta = meta, group.by = "labels")



# Calculate the aggregated cell-cell communication network
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling

cellchat_control@DB <- CellChatDB.use
cellchat_covid@DB <- CellChatDB.use
cellchat_ipf@DB <- CellChatDB.use

cellchat_control <- subsetData(cellchat_control) 
cellchat_covid <- subsetData(cellchat_covid) 
cellchat_ipf <- subsetData(cellchat_ipf) 

cellchat_control <- identifyOverExpressedGenes(cellchat_control)
cellchat_covid <- identifyOverExpressedGenes(cellchat_covid)
cellchat_ipf <- identifyOverExpressedGenes(cellchat_ipf)

cellchat_control <- identifyOverExpressedInteractions(cellchat_control)
cellchat_covid <- identifyOverExpressedInteractions(cellchat_covid)
cellchat_ipf <- identifyOverExpressedInteractions(cellchat_ipf)

cellchat_control <- projectData(cellchat_control, PPI.human)
cellchat_covid <- projectData(cellchat_covid, PPI.human)
cellchat_ipf <- projectData(cellchat_ipf, PPI.human)

cellchat_control <- computeCommunProb(cellchat_control)
cellchat_covid <- computeCommunProb(cellchat_covid)
cellchat_ipf <- computeCommunProb(cellchat_ipf)

cellchat_control <- filterCommunication(cellchat_control, min.cells = 100)
cellchat_covid <- filterCommunication(cellchat_covid, min.cells = 100)
cellchat_ipf <- filterCommunication(cellchat_ipf, min.cells = 100)

cellchat_control <- computeCommunProbPathway(cellchat_control)
cellchat_covid <- computeCommunProbPathway(cellchat_covid)
cellchat_ipf <- computeCommunProbPathway(cellchat_ipf)

cellchat_control <- aggregateNet(cellchat_control)
saveRDS(cellchat_control, here("EMT/RDS/cellchat_control_at2_fibro_smc.RDS"))

cellchat_covid <- aggregateNet(cellchat_covid)
saveRDS(cellchat_covid, here("EMT/RDS/cellchat_covid_at2_fibro_smc.RDS"))

cellchat_ipf <- aggregateNet(cellchat_ipf)
saveRDS(cellchat_ipf, here("EMT/RDS/cellchat_ipf_at2_fibro_smc.RDS"))
