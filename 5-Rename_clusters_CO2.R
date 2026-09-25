#Author: Debora Desantis (DD)
#Email: desantis@connect.hku.hk (DD)
#Affilitions: The Hong Kong University (HKU) - The Swire Institute for Marine Sciences (SWIMS)
#Copyright: Copyright 2026, Debora Desantis

###Rename clusters based on FindAllMarkers results and custom gene list###

setwd("C:/Users/debor/Desktop/PhD/snRNAseq_project/snRNAseq_seurat_analysis/Results_MT/")
dir <- "C:/Users/debor/Desktop/PhD/snRNAseq_project/snRNAseq_seurat_analysis/Results_MT/"

#Load required libraries
library(readxl)
library(stringr)
library(tibble)
library(tidyr)
library(dplyr)
library(Seurat)
library(ggplot2)
library(patchwork)

#Set seed for reproducibility
set.seed(42)

#Load data
Integrated <- readRDS("Int-CO2/Integrated_dataset_CO2_FindAllMarkers_ready_1.3.rds")

###---Rename clusters----
Integrated <- RenameIdents(object = Integrated, `0` = "Glut-0")
Integrated <- RenameIdents(object = Integrated, `1` = "GABA-1")
Integrated <- RenameIdents(object = Integrated, `2` = "Glut-2")
Integrated <- RenameIdents(object = Integrated, `3` = "Glut-3")
Integrated <- RenameIdents(object = Integrated, `4` = "Pallium-4")
Integrated <- RenameIdents(object = Integrated, `5` = "Glut-5")
Integrated <- RenameIdents(object = Integrated, `6` = "Glut-6")
Integrated <- RenameIdents(object = Integrated, `7` = "GABA-7")
Integrated <- RenameIdents(object = Integrated, `8` = "RGc-8")
Integrated <- RenameIdents(object = Integrated, `9` = "GABA-9")
Integrated <- RenameIdents(object = Integrated, `10` = "RGc-10")
Integrated <- RenameIdents(object = Integrated, `11` = "Glut-11")
Integrated <- RenameIdents(object = Integrated, `12` = "Glut-12")
Integrated <- RenameIdents(object = Integrated, `13` = "Glut-13")
Integrated <- RenameIdents(object = Integrated, `14` = "GABA-14")
Integrated <- RenameIdents(object = Integrated, `15` = "Glut-15")
Integrated <- RenameIdents(object = Integrated, `16` = "Glut-16")
Integrated <- RenameIdents(object = Integrated, `17` = "Glut-17")
Integrated <- RenameIdents(object = Integrated, `18` = "Glut-18")
Integrated <- RenameIdents(object = Integrated, `19` = "Glut-19")
Integrated <- RenameIdents(object = Integrated, `20` = "Glut-20")
Integrated <- RenameIdents(object = Integrated, `21` = "Glut-21")
Integrated <- RenameIdents(object = Integrated, `22` = "Glut-22")
Integrated <- RenameIdents(object = Integrated, `23` = "Glut-23")
Integrated <- RenameIdents(object = Integrated, `24` = "GABA-24")
Integrated <- RenameIdents(object = Integrated, `25` = "Oligo")
Integrated <- RenameIdents(object = Integrated, `26` = "GABA-26")
Integrated <- RenameIdents(object = Integrated, `27` = "Glut-27")
Integrated <- RenameIdents(object = Integrated, `28` = "Glut-28")
Integrated <- RenameIdents(object = Integrated, `29` = "Glut-29")
Integrated <- RenameIdents(object = Integrated, `30` = "GABA-30")
Integrated <- RenameIdents(object = Integrated, `31` = "OPCs")
Integrated <- RenameIdents(object = Integrated, `32` = "Perithelial")

#Add metadata column for CellType
Integrated$CellType <- Idents(Integrated)

saveRDS(Integrated, "Int-CO2/Integrated_dataset_CO2-FindAllMarkers_ready_1.3_names.rds")