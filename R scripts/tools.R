#### FUNCTIONS THAT'S USED DURING THE ANALYSIS ####
## Written by Serhat Beyaz / serhattbeyaz@gmail.com 


## Load required packages

library(pracma)
library(dplyr)
library(tidyr)



# Find overlapping cell lines using a mapping table
#
# This function finds the overlapping cell lines in two datasets after converting cell names based on a mapping table.
#
# Args:
#   dataset1: A data.frame containing the first dataset with a 'cell_line' column.
#   dataset2: A data.frame containing the second dataset with a 'cell_line' column.
#   mapping_table: A data.frame containing the mapping of cell line names with 'ID' and 'AC' columns.
#   cell_line_column: A string representing the name of the column to merge on in the datasets (default is "cell_line").
#   id_column: A string representing the name of the ID column in the mapping_table (default is "ID").
#   ac_column: A string representing the name of the AC column in the mapping_table (default is "AC").
#   relevant_columns: A character vector of relevant column names to keep in the final data frame.
#
# Returns:
#   A data.frame containing the overlapping cell lines in both datasets after mapping.

find_overlapping_cell_lines <- function(dataset1, dataset2, mapping_table, cell_line_column = "cell_line", id_column = "ID", ac_column = "AC", relevant_columns) {
  
  # Convert cell line names to lowercase in all datasets
  dataset1[, cell_line_column] <- tolower(dataset1[, cell_line_column])
  dataset2[, cell_line_column] <- tolower(dataset2[, cell_line_column])
  mapping_table[, id_column] <- tolower(mapping_table[, id_column])
  
  # Remove duplicates from the datasets and the mapping table
  dataset1 <- dataset1[!duplicated(dataset1[, cell_line_column]), ]
  dataset2 <- dataset2[!duplicated(dataset2[, cell_line_column]), ]
  mapping_table <- mapping_table[!duplicated(mapping_table[, id_column]), ]
  
  # Rename the columns in the mapping table to match the cell_line column name
  mapping_table_renamed <- mapping_table %>%
    dplyr::rename(!!cell_line_column := !!as.name(id_column),
                  new_cell_line := !!as.name(ac_column))
  
  # Merge dataset1 with mapping_table to convert cell_line names
  dataset1_mapped <- merge(dataset1, mapping_table_renamed, by = cell_line_column, all.x = TRUE) %>%
    dplyr::select(-cell_line, cell_line = new_cell_line)
  
  # Merge dataset2 with mapping_table to convert cell_line names
  dataset2_mapped <- merge(dataset2, mapping_table_renamed, by = cell_line_column, all.x = TRUE) %>%
    dplyr::select(-cell_line, cell_line = new_cell_line)
  
  # Find overlapping cell lines in the mapped datasets
  overlapping_cell_lines <- merge(dataset1_mapped, dataset2_mapped, by = cell_line_column)
  
  
  return(overlapping_cell_lines)
}





#' @name find_common_concentration_ranges
#'
#' @description This function reads two input files, combines them, and calculates the common
#' concentration ranges for each cell line and drug combination.
#'
#' @params:
#'   [file1]: A string representing the path to the first input file (e.g., "ccle_normalized.tsv").
#'   [file2]: A string representing the path to the second input file (e.g., "gdsc2_normalized.tsv").
#'
#' @return:
#'   A data.frame containing the common concentration ranges for each cell line and drug combination.

find_common_concentration_ranges <- function(file1, file2) {

  # Read input files
  data1 <- read.delim(file1) %>% drop_na(cell_line)
  data2 <- read.delim(file2) %>% drop_na(cell_line)
  
  # Add dataset identifier columns
  data1$dataset <- "data1"
  data2$dataset <- "data2"
  
  # Combine datasets and convert 'cell_line' and 'drug' columns to lower case
  combined_data <- rbind(data1, data2) %>% mutate_at(vars(cell_line, drug), tolower)
  
  # Calculate common dose ranges
  common_dose_ranges <- combined_data %>%
    # Filter out rows with non-positive dose values
    filter(dose > 0) %>%
    # Group data by 'cell_line', 'drug', and 'dataset'
    group_by(cell_line, drug, dataset) %>%
    # Calculate the minimum and maximum dose for each group
    summarise(min_dose = min(dose), max_dose = max(dose)) %>%
    # Pivot the dataset to have columns for min_dose and max_dose from both datasets
    pivot_wider(names_from = dataset, values_from = c(min_dose, max_dose)) %>%
    # Calculate the common minimum and maximum dose
    mutate(common_min_dose = pmax(min_dose_data1, min_dose_data2),
           common_max_dose = pmin(max_dose_data1, max_dose_data2)) %>%
    # Filter out rows with missing common_min_dose or common_max_dose values
    filter(!is.na(common_min_dose) & !is.na(common_max_dose)) %>%
    # Select the relevant columns
    select(cell_line, drug, common_min_dose, common_max_dose)
  
  return(common_dose_ranges)
}

common_conc = find_common_concentration_ranges("gdsc2_normalized_overlapping.tsv", "ccle_normalized_overlapping.tsv")


#####

# Define 4-parameter logistic regression model 
# Taken from the model.py of CurveCurator
logistic_4PL <- function(x, log_ec50, slope, front, back) {
  
  model = (front - back) / (1 + 10 ** (slope * (x - log_ec50))) + back
  
  return(model)
}


#' Compute AUC using fitted parameters from CurveCurator
#'
#' This function computes the Area Under the Curve (AUC) for a 4-parameter logistic 
#' (4PL) model using the fitted parameters generated by CurveCurator. It also 
#' provides an option to convert the concentration unit from micromolar (uM) to molar (M).
#'
#' @param min_conc A numeric value representing the minimum concentration.
#' @param max_conc A numeric value representing the maximum concentration.
#' @param ec50 A numeric value representing the half maximal effective concentration (EC50).
#' @param slope A numeric value representing the slope of the curve.
#' @param front A numeric value representing the response value at minimum concentration.
#' @param back A numeric value representing the response value at maximum concentration.
#' @param conc_unit A string representing the unit of concentration (default: "uM").
#'                  Supported units are "uM" (micromolar) and "M" (molar).
#'
#' @return A numeric value representing the normalized AUC of the 4PL model in the range of 0 to 1.
#'         Returns NA if the minimum concentration is missing (NA).
#'
#' @examples
#' compute_AUC(0.001, 100, 5, 1, 0, 100)
#' compute_AUC(1e-9, 1e-6, 5e-8, 1, 0, 100, conc_unit = "M")
compute_AUC = function(min_conc, max_conc, ec50, slope, front, back, conc_unit = "uM"){

  # Convert the unit
  if (conc_unit == "uM"){
    min_conc = min_conc * 1e-06
    max_conc = max_conc * 1e-06
  } 
  
  if(!is.na(min_conc)){ 
   # Concentration values
  min_concentration <- log10(min_conc) # Replace with your min tested concentration
  max_concentration <- log10(max_conc) # Replace with your max tested concentration
  num_points <- 1000 # Number of points to evaluate the function at
  
 
    # Generate a sequence of concentration points
    concentration_points <- seq(min_concentration, max_concentration, length.out = num_points)
    
    # Evaluate the 4PL model at the concentration points
    response_values <- logistic_4PL(concentration_points, ec50, slope, front, back)
    
    # Calculate AUC using the trapezoidal rule
    auc <- trapz(concentration_points, response_values) 
    
    model_front = response_values[1]
    
    # Normalize AUC to the range of 0 to 1
    max_auc <- model_front * (max(concentration_points) - min(concentration_points))
    normalized_auc <-(auc / max_auc)
    
    return(1 - normalized_auc)
    
  }else{
    return(NA)
  }
  
  
}




#' Compute AUC across cell lines for a specific drug
#'
#' This function computes the Area Under the Curve (AUC) for a specific drug across
#' multiple cell lines using the fitted curve parameters from a given source.
#' It takes into account the common dose ranges for each cell line and drug combination.
#'
#' @param drug A string representing the name of the drug.
#' @param source A string representing the the source, which is used as path to the directory containing the curve files.
#' @param common_dose_ranges A data.frame containing the common dose ranges for each cell line and drug combination. Output of the find_common_concentration_ranges()
#' @param p_filter An integer to be used as a threshold for filtering the curves based on their -Log p-values
#'
#' @return A list containing the computed AUC values for each cell line.
#'         Returns NA if no curve with a log.P > p_filter is found.
#'
#' @examples
#' # Assuming the common_dose_ranges data.frame is already available
#' auc_results = compute_auc_across_cells("lapatinib", "ccle", common_dose_ranges, 2)
compute_auc_across_cells = function(drug, source, common_dose_ranges, p_filter){
  

  curves = read.delim(paste0(source, "_overlaps", "/curves/", drug, "_curves.txt"))  
  
  curves = curves %>% filter(Curve.Log.P_Value.Corrected > p_filter) %>% mutate_at(vars(Name), tolower)
  
  if (nrow(curves) == 0){
    cat("No curve log.P > 2 is found!")
    return(NA)
  }else{
    
    auc_list = list()
    
    for (i in 1:nrow(curves)){
      
      if(drug %in% common_dose_ranges[common_dose_ranges$cell_line == curves$Name[i], ]$drug){
        
        conc_min =  common_dose_ranges[common_dose_ranges$cell_line == curves$Name[i] & common_dose_ranges$drug == drug, ]$common_min_dose
        
        conc_max =  common_dose_ranges[common_dose_ranges$cell_line == curves$Name[i] & common_dose_ranges$drug == drug, ]$common_max_dose
        
        params = curves[i, c("Log.EC50", "Curve.Slope", "Curve.Front", "Curve.Back")]
        
        
        auc = compute_AUC(min_conc = conc_min, max_conc = conc_max, ec50 = params$Log.EC50,
                          slope = params$Curve.Slope, front = params$Curve.Front, back = params$Curve.Back)
        
        
        auc_list[[curves$Name[i]]] = auc
        
      }else{
        next
      }
    }
    return(auc_list)
    
  }
}


#' Compute Median Absolute Deviation (MAD) 
#' 
#' This function computes MAD for a given drug.
#' 
#' @param auc_list List of AUC metrics across cell lines for a given drug
#' 
#' @return MAD value
compute_mad <- function(auc_list) {
  
  if (is.null(auc_list)) {
    return(NA)
  }
  
  # Filter out NULL or non-numeric values
  numeric_values <- sapply(auc_list, is.numeric)
  
  #filtered_values <- auc_list[which(numeric_values)]
  
  # Calculate MAD
  mads = mad(unlist(auc_list))

  return(mads)
}


#' Compute AUC and MAD for multiple drugs
#'
#' This function computes the Area Under the Curve (AUC) and Median Absolute Deviation (MAD) for multiple drugs.
#'
#' 
#' @param drug_list: A vector of strings containing the unique drug names.
#' @param source: A string representing the path to the source directory containing the curve files.
#'
#' @returns:
#'   A list containing two elements:
#'     - 'aucs': A named list containing the computed AUC values for each drug across cell lines.
#'     - 'mads': A named list containing the computed MAD values for each drug.

compute_auc_and_mad_for_drugs <- function(drug_list, source, common_dose_ranges, p_filter) {
  # Calculate the AUC for each drug
  aucs_test <- lapply(drug_list, function(x) {
    cat(x, " is started!\n")
    tryCatch({
      compute_auc_across_cells(x, source, common_dose_ranges, p_filter)
    },
    error = function(e) {
      message("Error occurred for ", x, e$message)
    })
  })
  
  # Name the list
  names(aucs_test) <- drug_list
  
  # Calculate MADs of AUCs for all drugs
  mads <- lapply(aucs_test, compute_mad)
  names(mads) <- names(aucs_test)
  
  return(list(aucs = aucs_test, mads = mads))
}

compute_auc_across_cells("erlotinib", source = "ccle", common_dose_ranges = common_conc, p_filter = 1.3)

