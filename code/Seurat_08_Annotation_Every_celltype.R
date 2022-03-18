library(Seurat)
library(here)
library(tidyverse)

refquery <- readRDS(here('explore/RDS/05_New_Integrate_df.RDS'))

Immune <- read.csv(here('explore/RDS/new_celltype/New_Immune_celltype'))
Epithelial <- read.csv(here('explore/RDS/new_celltype/New_Epithelial_celltype.csv'))
Endothelial <- read.csv(here('explore/RDS/new_celltype/New_Endothelial_celltype'))
Mesenchymal <- read.csv(here('explore/RDS/new_celltype/New_Mesenchymal_celltype'))


# Add celltypes 
celltypes <- rbind.data.frame(Immune, Epithelial, Endothelial, Mesenchymal)

# left join
refquery@meta.data <- left_join(refquery@meta.data, celltypes, by = 'cellbarcodes')
rownames(refquery@meta.data) <- refquery@meta.data$cellbarcodes

saveRDS(refquery, file = "07_df.RDS" )
