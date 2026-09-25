#Author: Debora Desantis (DD)
#Email: desantis@connect.hku.hk (DD)
#Affilitions: The Hong Kong University (HKU) - The Swire Institute for Marine Sciences (SWIMS)
#Copyright: Copyright 2026, Debora Desantis

###GLMM with Wald Test for cell cluster proportion###
#Samples: Int = Ambient CO2; CO2 = Elevated CO2

setwd("C:/Users/debor/Desktop/PhD/snRNAseq_project/snRNAseq_seurat_analysis/Results_MT/")
dir <- "C:/Users/debor/Desktop/PhD/snRNAseq_project/snRNAseq_seurat_analysis/Results_MT/"

library(lme4)
library(dplyr)
library(purrr)
library(broom.mixed) 
library(Matrix)
library(tibble)
library(ggplot2)
library(patchwork)

#Set seed for reproducibility
set.seed(42)

#Load data
Integrated <- readRDS("Int-CO2/Integrated_dataset_CO2-FindAllMarkers_ready_1.3_names.rds")

#Get number of cells per cluster and per sample of origin
cell_counts <- table(Integrated@meta.data$CellType, Integrated@meta.data$orig.ident)
cellNum <- as.matrix(cell_counts)

write.csv(cellNum, file = "Int-CO2/outs/cell_cluster_numbers.csv")

#Extract metadata and set ambient CO2 (Int) the reference baseline
meta_all <- Integrated@meta.data %>%
  mutate(Treatment = factor(Treatment, levels = c("Int", "CO2")),
         sample = as.factor(sample))
clusters <- unique(meta_all$CellType)

#Iterate the analysis through each cluster
stat_results <- map_dfr(clusters, function(cl) {
  
  df_cl <- meta_all %>%
    mutate(is_target = ifelse(CellType == cl, 1, 0))
  
  #Use tryCatch because GLMMs often struggle to converge with small sample sizes
  result <- tryCatch({
    
    #Fit Model
    model_full <- glmer(is_target ~ Treatment + (1|sample), 
                        data = df_cl, 
                        family = binomial,
                        control = glmerControl(optimizer = "bobyqa"))
    
    #Extract the coefficients table (Wald Test results)
    coef_summary <- summary(model_full)$coefficients
    
    #Extract stats for the "TreatmentCO2" coefficient
    log_OR <- coef_summary["TreatmentCO2", "Estimate"]
    p_val_wald <- coef_summary["TreatmentCO2", "Pr(>|z|)"]
    
    tibble(
      CellType = cl,
      log_odds_ratio = log_OR,
      p_value = p_val_wald)
    
  }, error = function(e) {
    return(tibble(CellType = cl, p_value = NA, log_odds_ratio = NA))
  })
  
  return(result)
})
#Loop ended

#Multiple Testing Correction
stat_results <- stat_results %>%
  filter(!is.na(p_value)) %>%
  mutate(
    q_value = p.adjust(p_value, method = "BH"),
    significant = q_value < 0.05,
    direction = case_when(
      log_odds_ratio > 0 & significant ~ "Increased in CO2",
      log_odds_ratio < 0 & significant ~ "Decreased in CO2",
      TRUE ~ "No significant change")) %>%
  arrange(q_value)

##Plot data
#Generate proportions from metadata
prop_df <- meta_all %>%
  group_by(sample, Treatment, CellType) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(sample) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()

#Add significance labels and metadata to stat_results
stat_results <- stat_results %>%
  mutate(
    sig_label = case_when(
      q_value < 0.001 ~ "***",
      q_value < 0.01  ~ "**",
      q_value < 0.05  ~ "*",
      TRUE            ~ "ns"),
    test_used = "GLMM-Wald") 

#Save 
write.csv(stat_results, "Int-CO2/outs/Integrated_CO2_GLMM_Wald_cell_stats.csv", row.names = FALSE)

#Font for Science
plot_font <- "sans" 

#Plot
clusters_to_plot <- stat_results$CellType 

plot_list <- map(clusters_to_plot, function(cl) {
  
  df_cl <- prop_df %>% filter(CellType == cl)
  sig   <- stat_results %>% filter(CellType == cl)
  
  y_pos <- max(df_cl$proportion) * 1.15
  
    ggplot(df_cl, aes(x = Treatment, y = proportion, fill = Treatment)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) +
      geom_point(aes(color = Treatment), position = position_dodge(0.8), size = 1.5) +
      scale_fill_manual(values = c("Int" = "cyan3", "CO2" = "darkorange3")) +
      scale_color_manual(values = c("Int" = "cyan4", "CO2" = "darkorange4")) +
    scale_y_continuous(
      labels = scales::percent_format(accuracy = 0.1),
      expand = expansion(mult = c(0.1, 0.2))) + 
    labs(
      title = cl,
      x = NULL, 
      y = "Proportion") +
    annotate("text", x = 1.5, y = y_pos,
             label = sig$sig_label, size = 6, fontface = "bold") +
      theme_classic(base_size = 8, base_family = plot_font)
})

#Arrange into grid and save plot
p_grid <- wrap_plots(plot_list, ncol = 4)