#Author: Debora Desantis (DD)
#Email: desantis@connect.hku.hk (DD)
#Affilitions: The Hong Kong University (HKU) - The Swire Institute for Marine Sciences (SWIMS)
#Copyright: Copyright 2026, Debora Desantis

###scDblFinder###
#Samples: Int = Ambient CO2; CO2 = Elevated CO2

setwd("C:/Users/debor/Desktop/PhD/snRNAseq_project/snRNAseq_seurat_analysis/Results_MT/")
dir <- "C:/Users/debor/Desktop/PhD/snRNAseq_project/snRNAseq_seurat_analysis/Results_MT/"

library(Seurat)
library(ggplot2)
library(dplyr)
library(sctransform)
library(glmGamPoi)
library(scDblFinder)
library(BiocParallel)

parallel::mclapply

#Set seed for reproducibility
set.seed(42)

#Load data previously pre-processed with SoupX and preQC filters
LD24Int <- readRDS("LD24-Int/LD24Int_scDblFinder_co2.rds")
LD26Int <- readRDS("LD26-Int/LD26Int_scDblFinder_co2.rds")
LD28Int <- readRDS("LD28-Int/LD28Int_scDblFinder_co2.rds")
LD1co2 <- readRDS("LD1-CO2/LD1CO2_scDblFinder.rds")
LD3co2 <- readRDS("LD3-CO2/LD3CO2_scDblFinder.rds")
LD6co2 <- readRDS("LD6-CO2/LD6CO2_scDblFinder.rds")

#Assign metadata
LD1co2$Treatment <- "CO2"
LD3co2$Treatment <- "CO2"
LD6co2$Treatment <- "CO2"
LD24Int$Treatment <- "Int"
LD26Int$Treatment <- "Int"
LD28Int$Treatment <- "Int"

#Assign metadata
LD1co2$sample <- "LD1CO2"
LD3co2$sample <- "LD3CO2"
LD6co2$sample <- "LD6CO2"
LD24Int$sample <- "LD24Int"
LD26Int$sample <- "LD26Int"
LD28Int$sample <- "LD28Int"

#Merge function in Seurat v5 keeps them separate as layers inside the RNA assay (already preserves per-sample identity via layers)
dt.list <- merge(LD1co2, c(LD3co2,LD6co2,LD24Int,LD26Int,LD28Int),
                 add.cell.ids = c("LD1co2", "LD3co2", "LD6co2",
                                  "LD24Int", "LD26Int", "LD28Int"))
#clean memory
rm(LD1co2, LD3co2, LD6co2, LD24Int, LD26Int, LD28Int)

#Violin plot to check data distribution
v1 <- VlnPlot(dt.list, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
              group.by = "orig.ident", pt.size = 0.01, ncol = 3) 

###---Run scDblFinder----
#Join layers and create data input
dt.list[["RNA"]] <- JoinLayers(dt.list[["RNA"]])
Layers(dt.list[["RNA"]])

sce <- as.SingleCellExperiment(dt.list)

###Run scDblFinder
sce <- scDblFinder(sce, samples = "sample", BPPARAM = SerialParam(RNGseed = 42))

#Visualize results
table(sce$scDblFinder.class)

#Distribution of doublet score (you want most values close to 0, a few to 1 and very little in between)
h <- hist(sce$scDblFinder.score)

#Explore results and add doublet info to Seurat object
meta_scdblfinder <- sce@colData@listData %>% as.data.frame() %>% 
  dplyr::select(starts_with("scDblFinder"))
head(meta_scdblfinder)
rownames(meta_scdblfinder) <- sce@colData@rownames
head(meta_scdblfinder)
dt.list <- AddMetaData(object = dt.list, metadata = meta_scdblfinder %>% 
                     dplyr::select("scDblFinder.class", "scDblFinder.score"))

#Re-check the number of doublets in the Seurat obj
table(dt.list$scDblFinder.class)

#Create a table with % of doublets per sample
doublets_summary <- dt.list@meta.data %>% 
  group_by(orig.ident, scDblFinder.class) %>% 
  summarise(total_count = n(),.groups = 'drop') %>% as.data.frame() %>% ungroup() %>%
  group_by(orig.ident) %>%
  mutate(countT = sum(total_count)) %>%
  group_by(scDblFinder.class, .add = TRUE) %>%
  mutate(percent = paste0(round(100 * total_count/countT, 2),'%')) %>%
  dplyr::select(-countT)
doublets_summary

#Remove doublets in the Seurat object
dt.list <- subset(dt.list, subset = scDblFinder.class == 'singlet')

#Check that the n. in your new Seurat object correspond to the n. of Singlets from DoubletFinder
check <- nrow(dt.list@meta.data)
check

#Re-split layers
dt.list[["RNA"]] <- split(dt.list[["RNA"]], f = dt.list$orig.ident)
Layers(dt.list[["RNA"]])

###----Save the output----
saveRDS(dt.list, file = "Int-CO2/co2_Int_noscDblets.rds")