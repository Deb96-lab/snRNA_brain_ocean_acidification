#Author: Debora Desantis (DD)
#Email: desantis@connect.hku.hk (DD)
#Affilitions: The Hong Kong University (HKU) - The Swire Institute for Marine Sciences (SWIMS)
#Copyright: Copyright 2025, Debora Desantis

###hdWGCNA analysis from "https://smorabit.github.io/hdWGCNA/articles/basic_tutorial.html"###
###EnrichGO analysis based on the hdWGCNA results###
#Samples: Int = Ambient CO2; CO2 = Elevated CO2

setwd("C:/Users/debor/Desktop/PhD/snRNAseq_project/snRNAseq_seurat_analysis/Results_MT/Int-CO2/")
dir <- "C:/Users/debor/Desktop/PhD/snRNAseq_project/snRNAseq_seurat_analysis/Results_MT/Int-CO2/"

library(org.Ldimidiatus.eg.db)
library(clusterProfiler)
library(ggplot2)
library(Seurat)
library(dplyr)
library(tidyverse)
library(tidytext) 
library(dplyr)
library(forcats)
library(rcartocolor)
library(WGCNA)
library(hdWGCNA)

#Font for Science
plot_font <- "sans"   

#Set seed for reproducibility
set.seed(42)

#Load data
Integrated <- readRDS("Integrated_dataset_CO2-FindAllMarkers_ready_1.3_names.rds")

###---Setup Seurat object for WGCNA----
Integrated <- SetupForWGCNA(Integrated,
                            gene_select = "fraction",
                            fraction = 0.05,
                            assay = "SCT",
                            wgcna_name = "Coexpression_Analysis")

#Construct Metacells
Integrated <- MetacellsByGroups(
  seurat_obj = Integrated,
  group.by = c("CellType", "sample"), 
  reduction = "integrated.rpca", 
  k = 25,          
  max_shared = 10, 
  min_cells = 100,
  ident.group = 'CellType')

#Normalize Metacells
Integrated <- NormalizeMetacells(Integrated)

#Save file for future use
saveRDS(Integrated, file='outs/hdWGCNA_CO2/hdWGCNA_object_CO2_MetacellsByGroups.rds')

#Set the expression data for these specific clusters
Integrated <- SetDatExpr(Integrated,
                         group_name = "GABA-14", 
                         group.by = "CellType",
                         assay = 'SCT', 
                         layer = 'data')
#Soft powers
Integrated <- TestSoftPowers(Integrated, networkType = 'signed')

#Plot the results to choose a power (look for the "elbow" where Scale-free Topology hits 0.8+)
plot_list <- PlotSoftPowers(Integrated)
wrap_plots(plot_list, ncol=2)

#Set the soft power (cluster specific)
soft_power <- 6

#Construct the network
Integrated <- ConstructNetwork(Integrated,
                               soft_power = soft_power,
                               overwrite_tom = TRUE,
                               tom_name = "GABA-14")
# Visualize the dendrogram
PlotDendrogram(Integrated, main='GABA-14 clusters hdWGCNA Dendrogram')

###---Relate the modules to my Treatment----
#Compute Module Eigengenes (MEs)
Integrated <- ScaleData(Integrated, features = VariableFeatures(Integrated))
Integrated <- ModuleEigengenes(Integrated) 

#Extract MEs and merge with metadata 
#ME with row as cells and columns as modules
MEs <- GetMEs(Integrated, harmonized=FALSE)
metadata <- GetMetacellObject(Integrated)@meta.data
Integrated@meta.data <- cbind(Integrated@meta.data, MEs)

#Compute eigengene-based connectivity (kME)
#Computes pairwise correlations between genes and module eigengenes
Integrated <- ModuleConnectivity(Integrated,
                                 group.by = 'CellType', 
                                 group_name = "GABA-14")

#Rename the modules
Integrated <- ResetModuleNames(Integrated, new_name = "GABA-14-M")

#Plot genes ranked by kME for each module
PlotKMEs(Integrated, ncol=4)

#Get the module assignment table
modules <- GetModules(Integrated) %>% subset(module != 'grey')
write.csv(modules, "outs/hdWGCNA_CO2/GABA-14/kME_ModuleConnectivity_GABA-14.csv")

#Get hub genes
hub_df <- GetHubGenes(Integrated, n_hubs = 20)
write.csv(hub_df, "outs/hdWGCNA_CO2/GABA-14/kME_ModuleConnectivity_HubGenes_GABA-14.csv")

###---DME analysis----
#Identify cell barcodes for each experimental group
group1 <- colnames(Integrated)[Integrated@meta.data$Treatment == "CO2"]
group2 <- colnames(Integrated)[Integrated@meta.data$Treatment == "Int"]

#Run the Differential Module Eigengene (DME) test
#This performs a Wilcoxon rank-sum test between the two barcode groups
DMEs_df <- FindDMEs(Integrated,
                    barcodes1 = group1,
                    barcodes2 = group2,
                    test.use = "wilcox", 
                    min.pct = 0.25,
                    wgcna_name = "Coexpression_Analysis")

#Save df
write.csv(DMEs_df, "outs/hdWGCNA_CO2/GABA-14/DMEs_CO2_GABA-14_hdWGCNA_logfc.csv")

#Generate a volcano plot for the modules
PlotDMEsVolcano(Integrated,
                DMEs_df,
                wgcna_name = "Coexpression_Analysis",
                plot_labels = TRUE)

###---Run EnrichGO funtion on the hdWGCNA results----
Integrated <- readRDS("Integrated_dataset_CO2-FindAllMarkers_ready_1.3_names.rds")
genes <- read.csv(file.path(dir, "outs/hdWGCNA_CO2/GABA-14/kME_ModuleConnectivity_GABA-14.csv"))

#Extract universe from Seurat
universe <- rownames(Integrated)

colnames(genes)[1] <- "GID"

genes_v <-genes$GID

### Gene Ontology (GO) over-representation analysis (ORA) from FindAllMarkers
ora_compare <- compareCluster(GID ~ module,
                              data = genes, 
                              fun = "enrichGO", 
                              OrgDb = org.Ldimidiatus.eg.db,
                              keyType = "GID",
                              ont = "ALL",
                              universe = universe,
                              pAdjustMethod = "BH",
                              pvalueCutoff = 0.05,
                              qvalueCutoff = 0.05,
                              minGSSize = 10,
                              maxGSSize = 500)

### Report results as data-frames 
results_all <- as.data.frame(ora_compare)

results_BP <- results_all[results_all$ONTOLOGY == "BP", ]

### Save results
write.csv(as.data.frame(results_all), 
          file = "outs/hdWGCNA_CO2/GABA-14/ORA_all_CO2_hdWGCNA_GABA-14.csv", row.names = F)

###---Repeat the analysis steps for all clusters of interest (RGc-8, RGc-10, GABA-9)----