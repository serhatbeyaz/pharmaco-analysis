library(ggplot2)
source("notebooks/tools.R")
library("dplyr")


create_scatter_plot_pubs <- function(drug_name) {
  
  # Load the data
  

  
  gdsc_auc = read.csv("Auc_values_GDSC2_published.csv", header = TRUE, sep = ",") %>% mutate_at(vars(drug), tolower)
  
  colnames(gdsc_auc) = tolower(colnames(gdsc_auc))
  
  gdsc_auc = gdsc_auc[gdsc_auc$drug == drug_name, ] %>% as.vector()
  
  ccle_auc = read.csv("Auc_values_CCLE_published.csv", header = TRUE, sep = ",") %>% mutate_at(vars(drug), tolower)
  
  colnames(ccle_auc) = tolower(colnames(ccle_auc))
  ccle_auc = ccle_auc[ccle_auc$drug == drug_name, ] %>% as.vector()
  
  #our_auc = compute_auc_across_cells(drug_name)
  
  intersects <- intersect(names(gdsc_auc), names(ccle_auc))
  
  # Filter both datasets to include only overlapping cell lines
  gdsc_auc <- gdsc_auc[intersects]  
  ccle_auc <- ccle_auc[intersects]
  
  # Remove NA values before creating the data frame
  na_filtered_indices <- !(is.na(unlist(gdsc_auc)) | is.na(unlist(ccle_auc)))
  gdsc_auc <- gdsc_auc[na_filtered_indices]
  ccle_auc <- ccle_auc[na_filtered_indices]
  
  # Prepare the data frame for plotting
  data <- data.frame(Published_CCLE = unlist(ccle_auc, use.names = FALSE),
                     Published_GDSC = unlist(gdsc_auc, use.names = FALSE),
                     #RecomputedAUCs = unlist(our_auc),
                     CellLines = names(ccle_auc))
  
  data = data[2:nrow(data), ]
  
  
  # Add labels based on the specified conditions
  data$label <- with(data, case_when(
    (Published_GDSC <= 0.2) & (Published_CCLE <= 0.2) ~ "Both resistant",
    (Published_GDSC > 0.2) & (Published_CCLE > 0.2) ~ "Both sensitive",
    (Published_GDSC > 0.2) & (Published_CCLE <= 0.2) ~ "GDSC sensitive, CCLE resistant",
    (Published_GDSC <= 0.2) & (Published_CCLE > 0.2) ~ "GDSC resistant, CCLE sensitive"
  ))
  
  label_colors <- c("Both resistant" = "red",
                    "Both sensitive" = "blue",
                    "GDSC sensitive, CCLE resistant" = "orange",
                    "GDSC resistant, CCLE sensitive" = "purple")
  
  # Fix the colors
  data <- data %>%
    mutate(Published_CCLE = as.numeric(Published_CCLE),
           Published_GDSC = as.numeric(Published_GDSC),
           label = factor(label))
  
  
  pearson_correlation <- cor(data$Published_GDSC, data$Published_CCLE, method = "pearson")
  spearman_correlation <- cor(data$Published_GDSC, data$Published_CCLE, method = "spearman")
  
  
  # Create the ggplot2 scatter plot
  plot <- ggplot(data, aes(x = Published_GDSC, y = Published_CCLE)) +
    geom_point(aes(color = label), size = 2.5, alpha = 0.6) +
    geom_hline(yintercept = 0.2, linetype = "dashed", color = "black") +
    geom_vline(xintercept = 0.2, linetype = "dashed", color = "black") +
    geom_abline(slope=1, intercept = 0)+
    labs(title = paste0("Correlation between Published AUCs and Recomputed AUCs, ", drug_name),
         x = "GDSC Published AUCs",
         y = "CCLE Published AUCs") +
    scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.1))+
    scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.1))+
    scale_color_manual(values = label_colors) +
    theme_bw() +
    theme(legend.position = "right")+
    annotate("text", x = 0.95, y = 0.15, label = paste("Pearson:", round(pearson_correlation, 2)), size = 4, hjust = 1) +
    annotate("text", x = 0.95, y = 0.05, label = paste("Spearman:", round(spearman_correlation, 2)), size = 4, hjust = 1)
  
  return(plot)
}


create_scatter_plot_recomp <- function(drug_name, common_dose_ranges, p_filter) {
  
  # Load the data
  ccle_auc = compute_auc_across_cells(drug_name, source = "ccle", common_dose_ranges, p_filter)
  gdsc_auc = compute_auc_across_cells(drug_name, source = "gdsc2", common_dose_ranges, p_filter)
  
  intersects <- intersect(names(gdsc_auc), names(ccle_auc))
  
  # Filter both datasets to include only overlapping cell lines
  gdsc_auc <- gdsc_auc[intersects]  
  ccle_auc <- ccle_auc[intersects]
  
  # Remove NA values before creating the data frame
  na_filtered_indices <- !(is.na(unlist(gdsc_auc)) | is.na(unlist(ccle_auc)))
  gdsc_auc <- gdsc_auc[na_filtered_indices]
  ccle_auc <- ccle_auc[na_filtered_indices]
  
  # Prepare the data frame for plotting
  data <- data.frame(Recomputed_CCLE = unlist(ccle_auc, use.names = FALSE),
                     Recomputed_GDSC = unlist(gdsc_auc, use.names = FALSE),
                     #RecomputedAUCs = unlist(our_auc),
                     CellLines = names(ccle_auc))
  
  
  # Add labels based on the specified conditions
  data$label <- with(data, case_when(
    (Recomputed_GDSC <= 0.2) & (Recomputed_CCLE <= 0.2) ~ "Both resistant",
    (Recomputed_GDSC > 0.2) & (Recomputed_CCLE > 0.2) ~ "Both sensitive",
    (Recomputed_GDSC > 0.2) & (Recomputed_CCLE <= 0.2) ~ "GDSC sensitive, CCLE resistant",
    (Recomputed_GDSC <= 0.2) & (Recomputed_CCLE > 0.2) ~ "GDSC resistant, CCLE sensitive"
  ))
  
  label_colors <- c("Both resistant" = "red",
                    "Both sensitive" = "blue",
                    "GDSC sensitive, CCLE resistant" = "orange",
                    "GDSC resistant, CCLE sensitive" = "purple")
  
  # Fix the colors
  data <- data %>%
    mutate(Recomputed_CCLE = as.numeric(Recomputed_CCLE),
           Recomputed_GDSC = as.numeric(Recomputed_GDSC),
           label = factor(label))

  pearson_correlation <- cor(data$Recomputed_GDSC, data$Recomputed_CCLE, method = "pearson")
  spearman_correlation <- cor(data$Recomputed_GDSC, data$Recomputed_CCLE, method = "spearman")
  
  cat("nrow", nrow(data))
  
  p_value = round(10**-p_filter, digits = 2)
  
  # Create the ggplot2 scatter plot
  plot <- ggplot(data, aes(x = Recomputed_GDSC, y = Recomputed_CCLE)) +
    geom_point(aes(color = label), size = 2.5, alpha = 0.6) +
    geom_hline(yintercept = 0.2, linetype = "dashed", color = "black") +
    geom_vline(xintercept = 0.2, linetype = "dashed", color = "black") +
    geom_abline(slope=1, intercept = 0)+
    labs(title = paste0("Correlation between Recomputed AUCs, ", drug_name),
         x = "GDSC Recomputed AUCs",
         y = "CCLE Recomputed AUCs") +
    scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.1))+
    scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.1))+
    scale_color_manual(values = label_colors) +
    theme_bw() +
    theme(legend.position = "right")+
    annotate("text", x = 0.95, y = 0.15, label = paste("Pearson:", round(pearson_correlation, 2)), size = 4, hjust = 1) +
    annotate("text", x = 0.95, y = 0.05, label = paste("Spearman:", round(spearman_correlation, 2)), size = 4, hjust = 1) + 
    annotate("text", x = 0.15, y = 0.95, label = paste0("P-filtering: < ", p_value), size = 3, hjust = 1)
  
  
  return(plot)
}


# Create a boxplot for MADs
create_mad_boxplot <- function(mads, source) {
  # Convert the list to a data frame and rename columns
  mad_df <- stack(mads) %>%
    rename(AUC = values, Drug = ind) %>%
    mutate(class = case_when(AUC > 0.13 ~ "cytotoxic",
                             AUC <= 0.13 ~ "targated"))
  
  # Create the boxplot
  plot <- ggplot(mad_df, aes(x = class, y = AUC, fill = class)) +
    geom_boxplot() +
    labs(title = paste0("Comparison of MAD of AUC values, ", source),
         x = "Drug class",
         y = "AUC") +
    geom_hline(yintercept = 0.13, linetype = "dashed", color = "red") +
    scale_y_continuous(breaks = c(seq(0, 1, 0.1), 0.13)) +
    theme_minimal()
  
  return(plot)
}


drugs = unique(common_conc$drug)

for (drug in drugs){
  
  recomp = create_scatter_plot_recomp(drug, common_dose_ranges = common_conc, p_filter = 1.3)
  
  pub = create_scatter_plot_pubs(drug_name = drug)
  
  grid = cowplot::plot_grid(pub, recomp, ncol = 1)
  
  pdf(paste0("correlation_plot/", drug, "_pub_vs_recomp.pdf"))
  print(grid)
  dev.off()
  
}


