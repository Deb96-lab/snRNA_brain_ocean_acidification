#Author: Debora Desantis (DD)
#Email: desantis@connect.hku.hk (DD)
#Affilitions: The Hong Kong University (HKU) - The Swire Institute for Marine Sciences (SWIMS)
#Copyright: Copyright 2026, Debora Desantis

###Define cluster resolution and run FindAllMarkers###
#Samples: Int = Ambient CO2; CO2 = Elevated CO2

setwd("C:/Users/debor/Desktop/PhD/snRNAseq_project/snRNAseq_seurat_analysis/Results_MT/")
dir <- "C:/Users/debor/Desktop/PhD/snRNAseq_project/snRNAseq_seurat_analysis/Results_MT/"

library(Seurat)
library(ggplot2)
library(dplyr)
library(sctransform)
library(glmGamPoi)
library(presto)
library(clustree)

parallel::mclapply

#Set seed for reproducibility
set.seed(42)

###----Cluster analysis and visualization of integrated data----
##Load RPCA Integrated data----
Integrated <- readRDS("Int-CO2/Integrated_dataset_CO2_rpca.rds")

#Apply a standard Seurat workflow
#Join layers 
Integrated[["RNA"]] <- JoinLayers(Integrated[["RNA"]])
Layers(Integrated[["RNA"]])

#Cluster cells RPCA Integration
Integrated <- FindNeighbors(Integrated, reduction = "integrated.rpca", assay = "SCT",
                            dims = 1:30, verbose = TRUE)

###---Run FindCluster at different resolutions----
Integrated <- FindClusters(Integrated, resolution = 0.5, cluster.name = "rpca_clusters_0.5", 
                           verbose = TRUE)
Integrated <- FindClusters(Integrated, resolution = 0.6, cluster.name = "rpca_clusters_0.6", 
                           verbose = TRUE)
Integrated <- FindClusters(Integrated, resolution = 0.7, cluster.name = "rpca_clusters_0.7", 
                           verbose = TRUE)
Integrated <- FindClusters(Integrated, resolution = 0.8, cluster.name = "rpca_clusters_0.8", 
                           verbose = TRUE)
Integrated <- FindClusters(Integrated, resolution = 1, cluster.name = "rpca_clusters_1", 
                           verbose = TRUE)
Integrated <- FindClusters(Integrated, resolution = 1.2, cluster.name = "rpca_clusters_1.2", 
                           verbose = TRUE)
Integrated <- FindClusters(Integrated, resolution = 1.3, cluster.name = "rpca_clusters_1.3", 
                           verbose = TRUE)
Integrated <- FindClusters(Integrated, resolution = 1.4, cluster.name = "rpca_clusters_1.4", 
                           verbose = TRUE)
Integrated <- FindClusters(Integrated, resolution = 1.5, cluster.name = "rpca_clusters_1.5", 
                           verbose = TRUE)
Integrated <- FindClusters(Integrated, resolution = 1.6, cluster.name = "rpca_clusters_1.6", 
                           verbose = TRUE)

#Plot the relationship between clustering at different resolution 
c <- clustree(Integrated, prefix = "rpca_clusters_", dims = 1:30)

#Select res = 1.3
Integrated <- FindClusters(Integrated, resolution = 1.3, cluster.name = "rpca_clusters_1.3", 
                           verbose = TRUE)  
Integrated <- RunUMAP(Integrated, reduction = "integrated.rpca", reduction.name = "umap.rpca",
                      dims = 1:30, verbose = TRUE)

#Plot UMAP 
p <- DimPlot(Integrated, reduction = "umap.rpca", label = TRUE, repel =TRUE)

#Save integrated dataset with clusters info
saveRDS(Integrated, "Int-CO2/Integrated_dataset_CO2_clusterInfo_rpca_1.3.rds")

###---FindAllMarkers----
##Prep dataset for FindAllMarkers function
Integrated <- PrepSCTFindMarkers(Integrated, assay = "SCT")

##Run FindAllMarkers function
Int_markers <- FindAllMarkers(Integrated, only.pos = TRUE, 
                              min.pct = 0.25, logfc.threshold = 0.25, test.use = "wilcox") 

#Save dataset ready for FindAllMArkers
saveRDS(Integrated, "Int-CO2/Integrated_dataset_CO2_FindAllMarkers_ready_1.3.rds")

#Select the top 10 markers and save the output
Top10 <- Int_markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)

write.csv(Top10, file = "Int-CO2/outs/Integrated_CO2_FindAllMarkers_rpca_top10.csv")

#Visualize to Top 10 markers
d1 <- DotPlot(Integrated, features = Top10, cols = c("lightgrey", "blue"),
              col.min = 0, dot.scale = 6)+ 
  RotatedAxis()