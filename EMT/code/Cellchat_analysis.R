library(here)
library(Seurat)  
library(tidyverse)
library(CellChat)


control <- readRDS(here("EMT/RDS/cellchat_EMT_control.RDS"))
covid <- readRDS(here("EMT/RDS/cellchat_EMT_covid.RDS"))
ipf <- readRDS(here("EMT/RDS/cellchat_EMT_ipf.RDS"))

object_list <- list(Control = control, COVID = covid)
control_covid <- mergeCellChat(object_list, add.names = names(object_list))

object_list <- list(Control = control, IPF = ipf)
control_ipf <- mergeCellChat(object_list, add.names = names(object_list))

object_list <- list(COVID = covid, IPF = ipf)
covid_ipf <- mergeCellChat(object_list, add.names = names(object_list))

object_list <- list(IPF = ipf, COVID = covid)
ipf_covid <- mergeCellChat(object_list, add.names = names(object_list))


par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(control_covid, weight.scale = T, measure = "count.merged", label.edge = T)
