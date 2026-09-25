#Author: Debora Desantis (DD)
#Email: desantis@connect.hku.hk (DD)
#Affilitions: The Hong Kong University (HKU) - The Swire Institute for Marine Sciences (SWIMS)
#Copyright: Copyright 2025, Debora Desantis

####CellChat analysis from "https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/Comparison_analysis_of_multiple_datasets.html"###
#Samples: Int = Ambient CO2; CO2 = Elevated CO2

setwd("C:/Users/debor/Desktop/PhD/snRNAseq_project/snRNAseq_seurat_analysis/Results_MT/")
dir <- "C:/Users/debor/Desktop/PhD/snRNAseq_project/snRNAseq_seurat_analysis/Results_MT/"

#Load required libraries
library(Seurat)
library(tidyverse)
library(CellChat)
library(patchwork)
library(ggplot2)

options(stringsAsFactors = FALSE)

#Set seed for reproducibility
set.seed(42)

#Load dataset
Integrated <- readRDS("Int-CO2/Integrated_dataset_CO2-FindAllMarkers_ready_1.3_names.rds")

###---First integrate Seurat obj with mouse annotation---- 
#Add mouse annotation
ortho_map <- read.csv("orthologs_with_genes_mouse.csv")

#Extract raw counts from your Seurat object
exp_matrix <- GetAssayData(Integrated, assay = "RNA", layer = "counts")

#Filter and Rename Genes based on Orthologs  to keep only genes present in mapping file
genes_to_keep <- intersect(rownames(exp_matrix), ortho_map$Ldim_ID)
exp_matrix <- exp_matrix[genes_to_keep, ]

#Create a named vector for remapping
mapping_vector <- ortho_map$Gene_Name
names(mapping_vector) <- ortho_map$Ldim_ID

#Update row names to Mouse Symbols
new_rownames <- mapping_vector[rownames(exp_matrix)]
exp_matrix <- exp_matrix[!is.na(new_rownames), ] # remove unmapped 
rownames(exp_matrix) <- new_rownames[!is.na(new_rownames)]

#Handle duplicates before creating Seurat obj
if (any(duplicated(rownames(exp_matrix)))) {
  message("Duplicate gene symbols detected. Aggregating counts...")
  
  #Aggregate using rowsum 
  exp_matrix <- rowsum(
    as.matrix(exp_matrix),
    group = rownames(exp_matrix))
}

#Create a clean Seurat object for CellChat
seurat_object <- CreateSeuratObject(counts = exp_matrix)

#Transfer metadata 
seurat_object@meta.data <- Integrated@meta.data[colnames(seurat_object), ]

#Split layers
seurat_object[["RNA"]] <- split(seurat_object[["RNA"]], f = seurat_object$orig.ident)
Layers(seurat_object[["RNA"]])
gc()

###---Run SCTransform pipeline----
seurat_object <- SCTransform(seurat_object, method = "glmGamPoi", 
                             vars.to.regress = c("percent.mt"), verbose = TRUE)

#Run PCA
seurat_object <- RunPCA(seurat_object, verbose = TRUE)

#Join layers 
seurat_object[["RNA"]] <- JoinLayers(seurat_object[["RNA"]])
Layers(seurat_object[["RNA"]])

###----Start CellChat Workflow----
#Subset and run CellChat for each treatment separately
split_list <- SplitObject(Integrated, split.by = "Treatment")

#Extract the individual objects from the list
seurat_co2 <- split_list[["CO2"]]
seurat_int  <- split_list[["Int"]]

Idents(seurat_co2) <- "CellType"
Idents(seurat_int)  <- "CellType"

#Exclude the unwanted clusters----
clusters_to_remove <- c("Perithelial", "Pallium-4")

seurat_int  <- subset(seurat_int,  idents = clusters_to_remove, invert = TRUE)
seurat_co2  <- subset(seurat_co2,  idents = clusters_to_remove, invert = TRUE)

#Drop the empty factor levels
Idents(seurat_int)  <- droplevels(Idents(seurat_int))
Idents(seurat_co2)  <- droplevels(Idents(seurat_co2))

seurat_int$CellType  <- droplevels(seurat_int$CellType)
seurat_co2$CellType  <- droplevels(seurat_co2$CellType)

###---Ambient CO2 treatment----
data.input_int <- seurat_int[["SCT"]]@data # normalized data matrix
labels_int <- Idents(seurat_int)
meta_int <- data.frame(labels = labels_int, row.names = names(labels_int))

seurat_int$samples <- seurat_int$orig.ident
cellChat_int <- createCellChat(object = seurat_int, group.by = "ident", assay = "SCT")

CellChatDB <- CellChatDB.mouse 

#showDatabaseCategory(CellChatDB)
CellChatDB.use <- CellChatDB 

#Assign the database to the cellChat object
cellChat_int@DB <- CellChatDB.use
cellChat_int <- subsetData(cellChat_int)
cellChat_int <- identifyOverExpressedGenes(cellChat_int)
cellChat_int <- identifyOverExpressedInteractions(cellChat_int)

ptm = Sys.time()
cellChat_int <- computeCommunProb(cellChat_int, type = "triMean", population.size = FALSE)

cellChat_int <- filterCommunication(cellChat_int, min.cells = 10)

cellChat_int <- computeCommunProbPathway(cellChat_int)

cellChat_int <- aggregateNet(cellChat_int)

#Extracting the inferred communications at the level of ligands/receptors
df.net_int <- subsetCommunication(cellChat_int)

#Extracting the inferred communications at the level of signaling pathways
df.netP_int <- subsetCommunication(cellChat_int, slot.name = "netP")

write.csv(df.netP_int, "Int-co2/outs/CellChat/CellChat_netP_Int.csv")
write.csv(df.net_int, "Int-co2/outs/CellChat/CellChat_net_Int.csv")
saveRDS(cellChat_int, file = "Int-co2/outs/CellChat/cellchat_Int.rds")

###---Elevated CO2 treatment----
data.input_co2 <- seurat_co2[["SCT"]]@data # normalized data matrix
labels_co2 <- Idents(seurat_co2)
meta_co2 <- data.frame(labels = labels_co2, row.names = names(labels_co2))

seurat_co2$samples <- seurat_co2$orig.ident
cellChat_co2 <- createCellChat(object = seurat_co2, group.by = "ident", assay = "SCT")

CellChatDB <- CellChatDB.mouse 
#showDatabaseCategory(CellChatDB)
CellChatDB.use <- CellChatDB 

#Assign the database to the cellChat object
cellChat_co2@DB <- CellChatDB.use
cellChat_co2 <- subsetData(cellChat_co2)
cellChat_co2 <- identifyOverExpressedGenes(cellChat_co2)
cellChat_co2 <- identifyOverExpressedInteractions(cellChat_co2)

ptm = Sys.time()
cellChat_co2 <- computeCommunProb(cellChat_co2, type = "triMean", population.size = FALSE)

cellChat_co2 <- filterCommunication(cellChat_co2, min.cells = 10)

cellChat_co2 <- computeCommunProbPathway(cellChat_co2)

cellChat_co2 <- aggregateNet(cellChat_co2)

#Extracting the inferred communications at the level of ligands/receptors
df.net_co2 <- subsetCommunication(cellChat_co2)

#Extracting the inferred communications at the level of signaling pathways
df.netP_co2 <- subsetCommunication(cellChat_co2, slot.name = "netP")

write.csv(df.netP_co2, "Int-co2/outs/CellChat/CellChat_netP_co2.csv")
write.csv(df.net_co2, "Int-co2/outs/CellChat/CellChat_net_co2.csv")
saveRDS(cellChat_co2, file = "Int-co2/outs/CellChat/cellchat_co2.rds")

###---Comparison between two datasets----
cellChat_int <- readRDS("Int-CO2/outs/CellChat/cellchat_Int.rds")
cellChat_co2 <- readRDS("Int-CO2/outs/CellChat/cellchat_co2.rds")

#compute the network centrality scores
cellChat_int <- netAnalysis_computeCentrality(cellChat_int, slot.name = "netP")
cellChat_co2 <- netAnalysis_computeCentrality(cellChat_co2, slot.name = "netP")

##merge objects
object.list <- list(Int = cellChat_int, CO2 = cellChat_co2)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
cellchat

##1. Compare the total number of interactions and interaction strength
ptm = Sys.time()
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
p1 <- gg1 + gg2

##2. Identify cell populations with significant changes in sending or receiving signals
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count) - diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size across datasets

gg_custom <- list()

for (i in 1:length(object.list)) {
  gg_custom[[i]] <- netAnalysis_signalingRole_scatter(
    object.list[[i]], 
    title = names(object.list)[i],
    weight.MinMax = weight.MinMax) +
    labs(
      x = "Outgoing interaction strength",
      y = "Incoming interaction strength") +
    theme_classic(base_size = 8, base_family = "sans") +
    theme(
      panel.grid.minor = element_blank(),
      axis.title.x = element_text(face = "bold", size = 7, color = "black"),
      axis.title.y = element_text(face = "bold", size = 7, color = "black"),
      axis.text = element_text(color = "black", size = 7),
      axis.text.x = element_text(angle = 0, hjust = 0.5), 
      axis.line = element_line(color = "black", linewidth = 0.4),
      axis.ticks = element_line(color = "black", linewidth = 0.4),
      legend.position = "right") +
    coord_cartesian(xlim = x.lim, ylim = y.lim)
}

#Combining Plots with patchwork
final_plot <- patchwork::wrap_plots(plots = gg_custom, ncol = 2) 
final_plot

##3. Compare the overall information flow of each signaling pathway or ligand-receptor pair
gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, 
               targets.use = NULL, stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, 
               targets.use = NULL, stacked = F, do.stat = TRUE)

gg1 + gg2

##4. Looking at incoming signals from target selected cell to various cluster cells
v <- netVisual_bubble(
  object = cellchat, 
  sources.use = "OPCs", 
  targets.use = c(1:31), 
  comparison = c(1, 2), 
  angle.x = 45,
  remove.isolate = TRUE) # Removes L-R pairs with no significant interactions

plot_font <- "sans" 

v <- v +
  theme_bw(base_size = 8, base_family = plot_font) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6, color = "cyan3"),
    axis.text.y = element_text(angle = 0, hjust = 1, size = 6, color = "black"),
    axis.title.x = element_text(size = 8, face = "bold", color = "black", margin = margin(t = 8)),
    legend.title = element_text(size = 8, angle = 0),
    legend.text  = element_text(size = 6, color = "black"),
    panel.grid   = element_blank(),
    axis.ticks   = element_blank()) +
  guides(fill = guide_colourbar(
    title.position = "top",
    barheight = unit(1.5, "cm"),
    barwidth  = unit(0.3, "cm")))
v

##5. Chord diagram for selected pathways across all cell clusters
pathways.show <- c("PTPR") 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", 
                      signaling.name = paste(pathways.show, names(object.list)[i]))
}
