library(here)
library(Seurat)
library(tidyverse)
library(CellChat)

df <- readRDS(here("EMT/RDS/emt_snRNA_new_annotation.RDS"))

DefaultAssay(df) <- "RNA"
control <- subset(df, subset = Diagnosis == "Control")
covid <- subset(df, subset = Diagnosis == "COVID-19")

#### subset by each diagnosis
# control
data_input <- GetAssayData(control, assay = "RNA", slot = "data")
Idents(control) <- "new_submain3"
labels <- Idents(control)
meta <- data.frame(labels = labels, row.names = names(labels))
cellchat_control <- createCellChat(object = data_input, meta = meta, group.by = "labels")

# covid
data_input <- GetAssayData(covid, assay = "RNA", slot = "data")
Idents(covid) <- "new_submain3"
labels <- Idents(covid)
meta <- data.frame(labels = labels, row.names = names(labels))
cellchat_covid <- createCellChat(object = data_input, meta = meta, group.by = "labels")




#### Add meta 
# control
data_input <- GetAssayData(control, assay = "RNA", slot = "data")
Idents(control) <- "new_submain3"
labels <- Idents(control)
meta <- data.frame(labels = labels, row.names = names(labels))
cellchat_control <- createCellChat(object = data_input, meta = meta, group.by = "labels")

# covid
data_input <- GetAssayData(covid, assay = "RNA", slot = "data")
Idents(covid) <- "new_submain3"
labels <- Idents(covid)
meta <- data.frame(labels = labels, row.names = names(labels))
cellchat_covid <- createCellChat(object = data_input, meta = meta, group.by = "labels")


# Calculate the aggregated cell-cell communication network
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling

cellchat_control@DB <- CellChatDB.use
cellchat_covid@DB <- CellChatDB.use

cellchat_control <- subsetData(cellchat_control) 
cellchat_covid <- subsetData(cellchat_covid) 

cellchat_control <- identifyOverExpressedGenes(cellchat_control)
cellchat_covid <- identifyOverExpressedGenes(cellchat_covid)

cellchat_control <- identifyOverExpressedInteractions(cellchat_control)
cellchat_covid <- identifyOverExpressedInteractions(cellchat_covid)

cellchat_control <- projectData(cellchat_control, PPI.human)
cellchat_covid <- projectData(cellchat_covid, PPI.human)

cellchat_control <- computeCommunProb(cellchat_control)
cellchat_covid <- computeCommunProb(cellchat_covid)

cellchat_control <- filterCommunication(cellchat_control, min.cells = 100)
cellchat_covid <- filterCommunication(cellchat_covid, min.cells = 100)

cellchat_control <- computeCommunProbPathway(cellchat_control)
cellchat_covid <- computeCommunProbPathway(cellchat_covid)

cellchat_control <- aggregateNet(cellchat_control)
saveRDS(cellchat_control, here("EMT/RDS/snRNA_cellchat_new_emt_control.RDS"))
cellchat_covid <- aggregateNet(cellchat_covid)
saveRDS(cellchat_covid, here("EMT/RDS/snRNA_cellchat_new_emt_covid.RDS"))












library(here)
library(Seurat)
library(tidyverse)
library(CellChat)

df <- readRDS(here("EMT/RDS/emt_new_annotation.RDS"))

DefaultAssay(df) <- "RNA"
control <- subset(df, subset = Diagnosis == "Control")
covid <- subset(df, subset = Diagnosis == "COVID")
ipf <- subset(df, subset = Diagnosis == "IPF")

#### subset by each diagnosis
# control
data_input <- GetAssayData(control, assay = "RNA", slot = "data")
Idents(control) <- "new_submain2"
labels <- Idents(control)
meta <- data.frame(labels = labels, row.names = names(labels))
cellchat_control <- createCellChat(object = data_input, meta = meta, group.by = "labels")

# covid
data_input <- GetAssayData(covid, assay = "RNA", slot = "data")
Idents(covid) <- "new_submain2"
labels <- Idents(covid)
meta <- data.frame(labels = labels, row.names = names(labels))
cellchat_covid <- createCellChat(object = data_input, meta = meta, group.by = "labels")

# covid
data_input <- GetAssayData(ipf, assay = "RNA", slot = "data")
Idents(ipf) <- "new_submain2"
labels <- Idents(ipf)
meta <- data.frame(labels = labels, row.names = names(labels))
cellchat_ipf <- createCellChat(object = data_input, meta = meta, group.by = "labels")



#### Add meta 
# control
data_input <- GetAssayData(control, assay = "RNA", slot = "data")
Idents(control) <- "new_submain2"
labels <- Idents(control)
meta <- data.frame(labels = labels, row.names = names(labels))
cellchat_control <- createCellChat(object = data_input, meta = meta, group.by = "labels")

# covid
data_input <- GetAssayData(covid, assay = "RNA", slot = "data")
Idents(covid) <- "new_submain2"
labels <- Idents(covid)
meta <- data.frame(labels = labels, row.names = names(labels))
cellchat_covid <- createCellChat(object = data_input, meta = meta, group.by = "labels")


data_input <- GetAssayData(ipf, assay = "RNA", slot = "data")
Idents(ipf) <- "new_submain2"
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
saveRDS(cellchat_control, here("EMT/RDS/ellchat_new_emt_control.RDS"))
cellchat_covid <- aggregateNet(cellchat_covid)
saveRDS(cellchat_covid, here("EMT/RDS/cellchat_new_emt_covid.RDS"))
cellchat_ipf <- aggregateNet(cellchat_ipf)
saveRDS(cellchat_ipf, here("EMT/RDS/cellchat_new_emt_ipf.RDS"))