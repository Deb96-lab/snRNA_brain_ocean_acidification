#Author: Debora Desantis (DD)
#Email: desantis@connect.hku.hk (DD)
#Affilitions: The Hong Kong University (HKU) - The Swire Institute for Marine Sciences (SWIMS)
#Copyright: Copyright 2026, Debora Desantis

###SoupX and data Preprocessing###
#Samples: Int = Ambient CO2; CO2 = Elevated CO2

setwd("C:/Users/debor/Desktop/PhD/snRNAseq_project/snRNAseq_seurat_analysis/Results_MT/")
dir <- "C:/Users/debor/Desktop/PhD/snRNAseq_project/snRNAseq_seurat_analysis/Results_MT/"

library(Seurat)
library(ggplot2)
library(dplyr)
library(SoupX)
library(sctransform)
library(glmGamPoi)
library(cowplot)

#Set seed for reproducibility
set.seed(42)

###----SoupX to remove ambient RNA from snRNA data. Done per each sample----
###----SoupX Auto
sc = load10X("LD1-CO2/outs")
sc = autoEstCont(sc)
out = adjustCounts(sc)

###----Pre-process QC per individual sample----
LD1co2_seurat <- CreateSeuratObject(counts = out, project = "LD1co2", 
                                    min.cells = 3)
LD1co2_seurat

#Add the %MT to the metadata
LD1co2_seurat[["percent.mt"]] <- PercentageFeatureSet(LD1co2_seurat, 
                                                      features = c("COXI", "COXII", "COXIII", "ND1", "ND2", 
                                                                   "ND3", "ND4", "ND5", "ND6", "ND4L",
                                                                   "ATPase8", "ATPase6", "Cytb"))

###----Generate QC plots to inspect data---- 
plot1 <- VlnPlot(LD1co2_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                 group.by = "orig.ident", ncol = 3)
plot1 

#plot the number of counts against the number of features 
plot2 <- FeatureScatter(LD1co2_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot2

#Plot the number of counts against the number of MT genes
plotMT <- FeatureScatter(LD1co2_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plotMT

###----Filter dataset----
LD1co2_seurat <- subset(LD1co2_seurat, subset = nFeature_RNA > 250 & nFeature_RNA < 5000 &
                          percent.mt < 3)
LD1co2_seurat

#Plot again to check the data distribution
plot3 <- VlnPlot(LD1co2_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                 ncol = 3, pt.size = 0.01, group.by = "orig.ident")
plot3

###----Save the output----
saveRDS(LD1co2_seurat, file = "LD1-CO2/LD1co2_scDblFinder.rds")

#Check cluster distribution within the sample
LD1co2_seurat <- SCTransform(LD1co2_seurat, method = "glmGamPoi", 
                             vars.to.regress = c("percent.mt"), verbose = TRUE)

LD1co2_seurat <- RunPCA(LD1co2_seurat, verbose = TRUE)
LD1co2_seurat <- RunUMAP(LD1co2_seurat, reduction = "pca", dims = 1:30, verbose = TRUE)
LD1co2_seurat <- FindNeighbors(LD1co2_seurat, reduction = "pca", dims = 1:30, verbose = TRUE)
LD1co2_seurat <- FindClusters(LD1co2_seurat, verbose = TRUE)

#Plot UMAP
un <- DimPlot(LD1co2_seurat, reduction = "umap", label = TRUE, repel = TRUE)

###Repeat for each individual sample and move to scDblFinder step###