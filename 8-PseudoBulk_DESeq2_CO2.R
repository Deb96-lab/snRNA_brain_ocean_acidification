#Author: Debora Desantis (DD)
#Email: desantis@connect.hku.hk (DD)
#Affilitions: The Hong Kong University (HKU) - The Swire Institute for Marine Sciences (SWIMS)
#Copyright: Copyright 2026, Debora Desantis

###DESeq2 pseudobulk analysis###
#Samples: Int = Ambient CO2; CO2 = Elevated CO2

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
library(DESeq2)

#Set seed for reproducibility
set.seed(42)

#Load dataset
Integrated <- readRDS("Int-CO2/Integrated_dataset_CO2-FindAllMarkers_ready_1.3_names.rds")

#Set the default assay for DE analysis
DefaultAssay(Integrated) <- "RNA"

#Aggregate counts by sample and CellType
pseudo_counts <- AggregateExpression(Integrated,
                                     group.by = c("CellType", "sample"),
                                     assays = "RNA",
                                     slot = "counts",
                                     return.seurat = FALSE)$RNA

#Prepare metadata for the Pseudobulk object 
metadata_pb <- Integrated@meta.data %>%
  select(sample, Treatment, CellType) %>%
  distinct() %>%
  mutate(pseudo_ID = paste(CellType, sample, sep = "_"))

#Align counts matrix columns with metadata rows
pseudo_counts <- pseudo_counts[, metadata_pb$pseudo_ID]

#Calculate cell counts per group for filtering
cell_counts <- as.data.frame(table(Integrated$CellType, Integrated$sample))
colnames(cell_counts) <- c("CellType", "sample", "n_cells")

#Initialize an empty list to hold the full results coming from the loop
all_results_list <- list()

###---Start the loop----
#Loop through each cell type assigned in the CellType column
for (ct in unique(Integrated$CellType)) {
  
  #1. Identify which columns in pseudo_counts belong to this CellType (>= 10 cells)
  #Deletes a whole sample from the analysis if it has <10 cells for that cluster.
  valid_samples <- cell_counts %>% 
    filter(CellType == ct, n_cells >= 10) %>%  
    pull(sample)
  
  #Split by treatment to have a more robust verification 
  ctrl_samples <- valid_samples[grepl("Int", valid_samples)]
  int_samples <- valid_samples[grepl("CO2", valid_samples)]
  
  #Ensure both groups have at least 2 replicates
  if (length(ctrl_samples) < 2 | length(int_samples) < 2) {
    message(paste("Skipping", ct, ": insufficient replicates in one or both groups."))
    next
  }
  
  #2. Extract specific columns using the exact names
  target_cols <- paste0(ct, "_", valid_samples)
  ct_counts <- pseudo_counts[, target_cols]
  
  #3. Align Metadata 
  ct_meta <- data.frame(
    pseudo_ID = target_cols,
    sample = valid_samples,
    stringsAsFactors = FALSE)
  
  #4. Pull the Treatment info from the original Integrated object metadata
  treatment_lookup <- Integrated@meta.data %>%
    select(sample, Treatment) %>%
    distinct()
  
  #Join them together
  ct_meta <- ct_meta %>%
    left_join(treatment_lookup, by = "sample")
  
  #Set the rownames to match the counts matrix columns
  rownames(ct_meta) <- ct_meta$pseudo_ID
  
  #5. Keep genes expressed in at least 5% of cells in the current cell type
  cells_in_ct <- which(Integrated$CellType == ct)
  pct_expr <- rowMeans(GetAssayData(Integrated, assay = "RNA", layer = "counts")[, cells_in_ct] > 0)
  keep_genes <- names(pct_expr[pct_expr >= 0.05]) 
  
  #Subset counts
  ct_counts <- ct_counts[keep_genes, ]
  ct_counts <- round(ct_counts)
  
  #6. Clean up Treatment column (remove any accidental NAs and ensures Treatment is a factor)
  ct_meta <- ct_meta[!is.na(ct_meta$Treatment), ]
  ct_counts <- ct_counts[, rownames(ct_meta)] 
  ct_meta$Treatment <- factor(ct_meta$Treatment, levels = c("Int", "CO2"))
  
  #6. Build DESeq2 object and run DESeq2
  dds <- DESeqDataSetFromMatrix(countData = ct_counts,
                                colData = ct_meta,
                                design = ~ Treatment)
  
  dds$Treatment <- relevel(dds$Treatment, ref = "Int")
  dds <- DESeq(dds, test = "Wald", quiet = TRUE)
  
  #7. Obtain results and apply filters
  res <- results(dds, 
                 contrast = c("Treatment", "CO2", "Int"), 
                 alpha = 0.05) 
  
  #Create a DF
  df_res <- as.data.frame(res)
  df_res$gene <- rownames(df_res)
  df_res$cell_type <- ct 
  
  #Order by padj
  df_res <- df_res[order(df_res$padj), ]
  
  #8. Store the full result in the empty list created before
  all_results_list[[ct]] <- df_res
  
  #Summary while running
  message(paste("Summary for Cell Type:", ct))
  summary(res)
}
##Loop ended

###---Create a DF with all results----
final_df <- do.call(rbind, all_results_list)

#Subset for significance across all cell types (DF with significant res only)
final_sig <- subset(final_df, padj < 0.05)

#Save files
write.csv(final_df, "Int-CO2/outs/DESeq2/Integrated_CO2_pseudobulk_full_res_1.3_padj_0.05.csv", row.names = FALSE)

###---Run EnrichGO analysis on the GABA-14 cluster DESeq2 results to find putative enriched pathways----
library(org.Ldimidiatus.eg.db)
library(clusterProfiler)

genes <- read.csv(file.path(dir, "Int-CO2/outs/DESeq2/gaba-14.csv"))

#Extract universe from Seurat
universe <- rownames(Integrated)
colnames(genes)[1] <- "GID"
genes_v <-genes$GID

###Gene Ontology (GO) over-representation analysis (ORA) 
ora <- enrichGO(gene = genes_v,
                        OrgDb = org.Ldimidiatus.eg.db,
                        keyType = "GID",
                        minGSSize = 10,
                        maxGSSize = 500,
                        universe = universe,
                        ont = "ALL",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05,
                        readable = TRUE)

### Report results as data-frames 
results_all <- as.data.frame(ora)

### Save results
write.csv(as.data.frame(results_all), 
          file = "Int-CO2/outs/DESeq2/GABA-14_ORA_res.csv", row.names = F)
