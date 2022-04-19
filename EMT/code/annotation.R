library(here)
library(Seurat)
library(tidyverse)

df <- readRDS(here("EMT/RDS/df.RDS"))

df@meta.data

unique(df$Diagnosis_tag)


# add donor_id 
df$row <- substr(rownames(df@meta.data), 1,6)
df$row2 <- substr(rownames(df@meta.data), 1, 4)


df@meta.data$Diagnosis_tag2 <- ifelse(df@meta.data$row %in% c("F01157"), "HC1", "NA")
df@meta.data$Diagnosis_tag2 <- ifelse(df@meta.data$row %in% c("F01174"), "HC2", df@meta.data$Diagnosis_tag2)
df@meta.data$Diagnosis_tag2 <- ifelse(df@meta.data$row %in% c("F01365", "F01366", "F01367"), "HC3", df@meta.data$Diagnosis_tag2)
df@meta.data$Diagnosis_tag2 <- ifelse(df@meta.data$row %in% c("HD65_A", "HD65_C", "HD65_G", "HD65_T"), "HC4", df@meta.data$Diagnosis_tag2)
df@meta.data$Diagnosis_tag2 <- ifelse(df@meta.data$row %in% c("HD66_A", "HD66_C", "HD66_G", "HD66_T"), "HC5", df@meta.data$Diagnosis_tag2)
df@meta.data$Diagnosis_tag2 <- ifelse(df@meta.data$row %in% c("HD67_A", "HD67_C", "HD67_G", "HD67_T"), "HC6", df@meta.data$Diagnosis_tag2)
df@meta.data$Diagnosis_tag2 <- ifelse(df@meta.data$row %in% c("HD68_A", "HD68_C", "HD68_G", "HD68_T"), "HC7", df@meta.data$Diagnosis_tag2)
df@meta.data$Diagnosis_tag2 <- ifelse(df@meta.data$row %in% c("F00409"), "HC8", df@meta.data$Diagnosis_tag2)
df@meta.data$Diagnosis_tag2 <- ifelse(df@meta.data$row %in% c("HD70_A", "HD70_C", "HD70_G", "HD70_T"), "HC9", df@meta.data$Diagnosis_tag2)
df@meta.data$Diagnosis_tag2 <- ifelse(df@meta.data$row %in% c("F01394"), "HC10", df@meta.data$Diagnosis_tag2)


df@meta.data$Diagnosis_tag2 <- ifelse(df@meta.data$row %in% c("F00431"), "I1", df@meta.data$Diagnosis_tag2)
df@meta.data$Diagnosis_tag2 <- ifelse(df@meta.data$row %in% c("F01302", "F01303"), "I2", df@meta.data$Diagnosis_tag2)
df@meta.data$Diagnosis_tag2 <- ifelse(df@meta.data$row %in% c("F01214"), "I3", df@meta.data$Diagnosis_tag2)
df@meta.data$Diagnosis_tag2 <- ifelse(df@meta.data$row %in% c("F01172", "F01173"), "I4", df@meta.data$Diagnosis_tag2)
df@meta.data$Diagnosis_tag2 <- ifelse(df@meta.data$row %in% c("F01379", "F01380"), "I5", df@meta.data$Diagnosis_tag2)
df@meta.data$Diagnosis_tag2 <- ifelse(df@meta.data$row %in% c("ILD53_"), "I6", df@meta.data$Diagnosis_tag2)
df@meta.data$Diagnosis_tag2 <- ifelse(df@meta.data$row %in% c("ILD59-"), "I7", df@meta.data$Diagnosis_tag2)
df@meta.data$Diagnosis_tag2 <- ifelse(df@meta.data$row %in% c("ILD60-"), "I8", df@meta.data$Diagnosis_tag2)
df@meta.data$Diagnosis_tag2 <- ifelse(df@meta.data$row %in% c("ILD61-"), "I9", df@meta.data$Diagnosis_tag2)
df@meta.data$Diagnosis_tag2 <- ifelse(df@meta.data$row %in% c("ILD63_"), "I10", df@meta.data$Diagnosis_tag2)
df@meta.data$Diagnosis_tag2 <- ifelse(df@meta.data$row %in% c("F01390"), "I11", df@meta.data$Diagnosis_tag2)
df@meta.data$Diagnosis_tag2 <- ifelse(df@meta.data$row %in% c("F01391", "F01392"), "I12", df@meta.data$Diagnosis_tag2)




df@meta.data$Diagnosis_tag2 <- ifelse(df@meta.data$row %in% c("obj1_A", "obj1_C", "obj1_G" ,"obj1_T"), "C1_1", df@meta.data$Diagnosis_tag2)
df@meta.data$Diagnosis_tag2 <- ifelse(df@meta.data$row %in% c("obj2_A", "obj2_C", "obj2_G", "obj2_T"), "C1_2", df@meta.data$Diagnosis_tag2)
df@meta.data$Diagnosis_tag2 <- ifelse(df@meta.data$row %in% c("obj3_A" ,"obj3_C", "obj3_G" ,"obj3_T"), "C1_3", df@meta.data$Diagnosis_tag2)
df@meta.data$Diagnosis_tag2 <- ifelse(df@meta.data$row %in% c("obj4_A", "obj4_C" ,"obj4_G", "obj4_T"), "C1_4", df@meta.data$Diagnosis_tag2)

df <- subset(df, subset = seurat_clusters != 29)
saveRDS(df, file = here("EMT/RDS/df.RDS"))








