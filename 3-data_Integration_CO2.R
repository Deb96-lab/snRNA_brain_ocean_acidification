#Author: Debora Desantis (DD)
#Email: desantis@connect.hku.hk (DD)
#Affilitions: The Hong Kong University (HKU) - The Swire Institute for Marine Sciences (SWIMS)
#Copyright: Copyright 2026, Debora Desantis

###Integrate samples with RPCA pipeline###
#Samples: Int = Ambient CO2; CO2 = Elevated CO2

setwd("C:/Users/debor/Desktop/PhD/snRNAseq_project/snRNAseq_seurat_analysis/Results_MT/")
dir <- "C:/Users/debor/Desktop/PhD/snRNAseq_project/snRNAseq_seurat_analysis/Results_MT/"

library(Seurat)
library(ggplot2)
library(dplyr)
library(sctransform)
library(glmGamPoi)

parallel::mclapply

#Set seed for reproducibility
set.seed(42)

###----Normalize and integrate libraries----
#Load data
dt.list <- readRDS("Int-CO2/co2_Int_noscDblets.rds")

###---Normalize with SCTransform and identify variable features per each dataset independently----
dt.list <- SCTransform(dt.list, method = "glmGamPoi", 
                       vars.to.regress = c("percent.mt"), verbose = TRUE) 

#Violin plot to check data distribution after SCTransform
v1 <- VlnPlot(dt.list, features = c("nFeature_SCT", "nCount_SCT", "percent.mt"), 
              group.by = "orig.ident", pt.size = 0.01, ncol = 3) 

#Run PCA
dt.list <- RunPCA(dt.list, verbose = TRUE)

##Check how many PCs you need to use with an Elbowplot
ElbowPlot(dt.list, ndims= 50)

#Run UMAP before integration to have a look at the non-integrated data
dt.list <- FindNeighbors(dt.list, reduction = "pca",
                         dims = 1:30, verbose = TRUE)

dt.list <- FindClusters(dt.list, cluster.name = "non_Integrated_cluster", 
                        verbose = TRUE)

dt.list <- RunUMAP(dt.list, reduction = "pca", reduction.name = "umap.non_Integrated",
                   dims = 1:30, verbose = TRUE)

un <- DimPlot(dt.list, reduction = "umap.non_Integrated", split.by = "orig.ident",  
              label = TRUE, repel = TRUE)

###----Integrated Analysis----
##---RPCA Integration----
dt.list <- IntegrateLayers(object = dt.list, method = RPCAIntegration, 
                               normalization.method = "SCT", orig.reduction = "pca", 
                               new.reduction = "integrated.rpca")

#Join layers 
dt.list[["RNA"]] <- JoinLayers(dt.list[["RNA"]])
Layers(dt.list[["RNA"]])

#Save RPCA integrated dataset
saveRDS(dt.list, "Int-CO2/Integrated_dataset_CO2_rpca.rds")