#Author: Debora Desantis (DD)
#Email: desantis@connect.hku.hk (DD)
#Affilitions: The Hong Kong University (HKU) - The Swire Institute for Marine Sciences (SWIMS)
#Copyright: Copyright 2026, Debora Desantis

###scGSEA analysis from DESeq2 results###
#Samples: Int = Ambient CO2; CO2 = Elevated CO2
#Genome annotation file located at "https://www.schunterlab.com/resources"

setwd("C:/Users/debor/Desktop/PhD/snRNAseq_project/snRNAseq_seurat_analysis/Results_MT/")
dir <- "C:/Users/debor/Desktop/PhD/snRNAseq_project/snRNAseq_seurat_analysis/Results_MT/"

library(fgsea)
library(data.table)
library(ggplot2)
library(presto)
library(Seurat)
library(dplyr)
library(stringr)
library(tidyr)
library(readxl)
library(stringr)
library(tibble)
library(GO.db)
library(tibble)
library(purrr)
library(tidyverse)   
library(tidytext)

#Set seed for reproducibility
set.seed(42)

#Load data
Integrated <- readRDS("Integrated_dataset_ctrl-Int_FindAllMarkers_ready_0.8_NoCluster13_names.rds")

###---Obtain custom GO list (custom_GO_set) for L. dimidiatus from the original annotation file present at "https://www.schunterlab.com/resources"----
###---DESeq2 GSEA----
res_all <- read.csv("Int-CO2/outs/DESeq2/Integrated_CO2_pseudobulk_full_res_1.3_padj_0.05.csv",
                    stringsAsFactors = FALSE)

colnames(res_all)[1] <- "GID"
gene_id_col <- "GID"  
stat_col <- "stat"

make_ranked_stats <- function(df) {
  df <- df %>%
    filter(
      !is.na(.data[[gene_id_col]]), .data[[gene_id_col]] != "",
      !is.na(.data[[stat_col]]),
      is.finite(.data[[stat_col]]))   # fgsea() cannot take NA/Inf/-Inf
  df <- df %>% filter(baseMean >= 10) #remove lowly expressed genes

  ##collapse duplicate gene IDs - keep whichever row has the larger |stat|
  df <- df %>%
    group_by(.data[[gene_id_col]]) %>%
    slice_max(order_by = abs(.data[[stat_col]]), n = 1, with_ties = FALSE) %>%
    ungroup()
  ranks <- setNames(df[[stat_col]], df[[gene_id_col]])
  sort(ranks, decreasing = TRUE)            # positive (up) at the top
}

ranked_stats_list <- res_all %>%
  split(.$cell_type) %>%
  map(make_ranked_stats)

##sanity check: no NA/Inf, no duplicate names, in any cluster
walk(ranked_stats_list, function(r) {
  stopifnot(all(is.finite(r)), !anyDuplicated(names(r)))
})

#Run fgsea analysis
fgsea_results_list <- imap(ranked_stats_list, function(ranks, ct) {
    
    res <- fgsea(
      pathways = custom_GO_set, #obtained from "https://www.schunterlab.com/resources"
      stats    = ranks,
      minSize  = 25,
      maxSize  = 500,
      eps      = 0, 
      nproc = 1)
    
    res %>%
      arrange(padj) %>%
      mutate(
        cell_type   = ct,
        leadingEdge = map_chr(leadingEdge, ~ paste(.x, collapse = ","))  # flatten list-col for export
      ) %>%
      relocate(cell_type, .before = 1)
  })
  
fgsea_all <- bind_rows(fgsea_results_list)

#Save results  
write.csv(fgsea_all, "Int-CO2/outs/DESeq2/fgsea_results/DEseq2_fgsea_significant_only_descriptions.csv", 
          row.names = FALSE)  #Table S11 within Supplementary Datasets

###---Plot scGSEA heatmap----
input_file    <- "Int-CO2/outs/DESeq2/fgsea_results/DEseq2_fgsea_significant_only_descriptions.csv"
padj_cutoff   <- 0.05
min_clusters  <- 5      # a pathway must be significant in >=
top_n         <- 40
plot_font <- "sans"    

df <- read_csv(input_file, show_col_types = FALSE) %>%
  mutate(
    padj = as.numeric(padj),   
    NES  = as.numeric(NES)) %>%
  filter(!is.na(padj), !is.na(NES))

#Identify pathways significant in >= min_clusters cell types
sig_df <- df %>% filter(padj < padj_cutoff)

pathway_hits <- sig_df %>%
  distinct(GO_Description, cell_type) %>%
  count(GO_Description, name = "n_clusters") %>%
  filter(n_clusters >= min_clusters)

#Rank pathways by the smallest padj 
pathway_rank <- sig_df %>%
  filter(GO_Description %in% pathway_hits$GO_Description) %>%
  group_by(GO_Description) %>%
  summarise(best_padj = min(padj), .groups = "drop") %>%
  left_join(pathway_hits, by = "GO_Description") %>%
  arrange(best_padj, desc(n_clusters))

selected_pathways <- pathway_rank %>%
  slice_head(n = top_n) %>%
  pull(GO_Description)

#Build the pathway x cell_type NES matrix
mat_df <- df %>%
  filter(GO_Description %in% selected_pathways) %>%
  mutate(NES_masked = ifelse(padj < padj_cutoff, NES, NA_real_)) %>%
  select(GO_Description, cell_type, NES_masked) %>%
  distinct(GO_Description, cell_type, .keep_all = TRUE) %>%
  pivot_wider(names_from = cell_type, values_from = NES_masked)

mat <- as.matrix(mat_df[, -1])
rownames(mat) <- mat_df$GO_Description
mat <- mat[selected_pathways, , drop = FALSE]
mat_for_clustering <- mat
mat_for_clustering[is.na(mat_for_clustering)] <- 0

row_dist <- dist(mat_for_clustering, method = "euclidean")
col_dist <- dist(t(mat_for_clustering), method = "euclidean")

row_hc <- hclust(row_dist, method = "complete")
col_hc <- hclust(col_dist, method = "complete")

row_order <- rownames(mat)[row_hc$order]
col_order <- colnames(mat)[col_hc$order]

#Reshape to long format for ggplot 
plot_df <- mat %>%
  as.data.frame() %>%
  rownames_to_column("GO_Description") %>%
  pivot_longer(-GO_Description, names_to = "cell_type", values_to = "NES") %>%
  mutate(GO_Description = factor(GO_Description, levels = rev(row_order)),
         cell_type = factor(cell_type, levels = col_order))

#Colour scale 
max_abs <- max(abs(plot_df$NES), na.rm = TRUE)
n_row <- nlevels(plot_df$GO_Description)

#Plot
p <- ggplot(plot_df, aes(x = cell_type, y = GO_Description, fill = NES)) +
  geom_tile(colour = "white", linewidth = 0.28) +
  scale_fill_gradient2(
    low      = "purple4",      
    mid      = "lightyellow1",   
    high     = "darkgreen",       
    midpoint = 0,
    limits   = c(-max_abs, max_abs),
    na.value = "white",     
    name     = "NES",
    guide    = guide_colourbar(barwidth = unit(0.3, "cm"), barheight = unit(2, "cm"))) +
  labs(x = "Cell type", y = "Pathways") +
  theme_minimal(base_family = plot_font, base_size = 8)