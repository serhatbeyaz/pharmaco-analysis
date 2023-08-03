################################################################################
# Copyright (c) 2015, 2016, 2017, 2018 Genome Research Ltd. 
# Copyright (c) 2015, 2016, 2017, 2018 The Netherlands Cancer Institute (NKI)
#  
# Author: Howard Lightfoot <cancerrxgene@sanger.ac.uk> 
# Author: Dieudonne van der Meer
# Author: Daniel J Vis
# 
# This file is part of gdscIC50. 
# 
# gdscIC50 is free software: you can redistribute it and/or modify it under 
# the terms of the GNU General Public License as published by the Free Software 
# Foundation; either version 3 of the License, or (at your option) any later 
# version. 
#  
# This program is distributed in the hope that it will be useful, but WITHOUT 
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more 
# details. 
#  
# You should have received a copy of the GNU General Public License along with 
# this program. If not, see <http://www.gnu.org/licenses/>. 
################################################################################

# I (Serhat BEYAZ), copied the code from https://github.com/CancerRxGene/gdscIC50/blob/master/R/nlme_fit_prep.R
# and adapted it to normalize flurosence intensity and convert it into cell viability measures.



#' Removes well positions where the drug is now considered unsuitable for
#'  screening from GDSC raw data, i.e., TAG = 'FAIL'
#' 
#' \code{removeFailedDrugs} removes rows from GDSC raw data where the
#'  \code{TAG == 'UN-USED'}.
#' 
#' @param myDat a GDSC raw data data frame.
#' 
#' @seealso \code{\link{removeMissingDrugs}}, \code{\link{normalizeData}},
#'  \code{\link{setConcsForNlme}}, \code{\link{prepNlmeData}}
#'
#' @examples
#' data("gdsc_example")
#' gdsc_example <- removeFailedDrugs(gdsc_example)
#' gdsc_example <- removeMissingDrugs(gdsc_example)
#' gdsc_example <- normalizeData(gdsc_example)
#' gdsc_example <- setConcsForNlme(gdsc_example)
#' nlme_data <- prepNlmeData(gdsc_example, "COSMIC_ID")
#'  
#' @export
removeFailedDrugs <- function(myDat){
  # Remove the entire position because the failed drug might be part of a combination.
  failed_positions <- myDat %>% 
    filter_(~TAG == 'UN-USED') %>%
    select_(~SCAN_ID, ~POSITION)
  myDat <-  anti_join(myDat, failed_positions, by = c("SCAN_ID", "POSITION"))
  return(myDat)
}

#' Removes library drugs listed with a drug id of NA from GDSC raw data. 
#' 
#' In contrast to drugs with 'FAIL' tags, the drugset contained NA for 
#'  the DRUG_ID, i.e., no drug id was assigned for that tag in the experimental
#'  design.
#' 
#' \code{removeMissingDrugs} removes rows from GDSC raw data where the
#'  \code{DRUG_ID} is NA.
#' 
#' @param myDat a GDSC raw data data frame.
#' 
#' @seealso  \code{\link{removeFailedDrugs}},  \code{\link{normalizeData}},
#'   \code{\link{setConcsForNlme}},  \code{\link{prepNlmeData}}
#'  
#' @examples
#' data("gdsc_example")
#' gdsc_example <- removeFailedDrugs(gdsc_example)
#' gdsc_example <- removeMissingDrugs(gdsc_example)
#' gdsc_example <- normalizeData(gdsc_example)
#' gdsc_example <- setConcsForNlme(gdsc_example)
#' nlme_data <- prepNlmeData(gdsc_example, "COSMIC_ID")
#' 
#' @export
removeMissingDrugs <- function(myDat){
  na_libs <- myDat %>%
    filter_(~grepl("^(L|R|A)\\d+", TAG)) %>% 
    filter_(~is.na(DRUG_ID))
  myDat <- anti_join(myDat, na_libs, by = c("SCAN_ID", "POSITION"))
  return(myDat)
}


calcTagMean <- function(myDat, tag_name, mean_col_name = "tag_mean") {
  suppressMessages(check_for_tag <- left_join(myDat %>% 
                                                select_(~SCAN_ID) %>% 
                                                distinct(),
                                              myDat %>% 
                                                group_by_(~SCAN_ID) %>% 
                                                filter_(~TAG == tag_name) %>% 
                                                count()
  )
  )
  
  e1 <- simpleError(paste("calcTagMean:", 
                          tag_name, 
                          "is not present for some or all of the SCAN_IDs in your data.", 
                          sep = " "))
  if (any(is.na(check_for_tag$n))){
    stop(e1)
  }
  
  tag_means <- myDat %>% 
    group_by_(~SCAN_ID) %>% 
    filter_(~TAG == tag_name) %>% 
    summarise_(tag_mean = ~mean(INTENSITY, na.rm = T))
  
  e2 <- simpleError(paste("calcTagMean:", 
                          tag_name, 
                          "has a mean of NaN for some or all of the SCAN_IDs in your data.", 
                          sep = " "))
  if (any(is.nan(tag_means$tag_mean))){
    stop(e2)
  }
  
  tag_means <- tag_means %>% rename_(.dots = stats::setNames("tag_mean", mean_col_name))
  
  return(tag_means)
}

#' Normalizes GDSC raw data intensities with respect to controls.
#' 
#' \code{normalizeData} returns normalized intensities for the drug treated
#'   wells - column \code{normalized_intensity}. Replaces \code{TAG} column
#'   with \code{lib_drug} and \code{dose} columns.
#'
#' @param myDat a GDSC raw data data frame.
#' @param trim logical indicating whether to trim normalized values to the range
#' 0 to 1. default \code{(trim = T)}
#' @param neg_control The tag used to recognise a negative control well - the 
#'   upper end of the dynamic range.
#' @param pos_control The tag used to recognise a positive control well - the 
#'   lower end of the dynamic range.
#' 
#' @seealso  \code{\link{removeFailedDrugs}},  \code{\link{removeMissingDrugs}},
#'   \code{\link{setConcsForNlme}},  \code{\link{prepNlmeData}}
#' 
#' @examples
#' data("gdsc_example")
#' gdsc_example <- removeFailedDrugs(gdsc_example)
#' gdsc_example <- removeMissingDrugs(gdsc_example)
#' gdsc_example <- normalizeData(gdsc_example)
#' gdsc_example <- setConcsForNlme(gdsc_example)
#' nlme_data <- prepNlmeData(gdsc_example, "COSMIC_ID")
#'
#' @export
normalizeData <- function(myDat, trim = T, neg_control = 'NC-0',
                          pos_control = 'B'){
  
  nc1 <- calcTagMean(myDat, tag_name = neg_control, mean_col_name = "NC")
  pc1 <- calcTagMean(myDat, tag_name = pos_control, mean_col_name = "PC")
  
  # Take account of historic data with no MASTER_CELL_ID column - use ends_with
  normalized_data <- myDat %>% 
    filter_(~ grepl("(A|L|R)\\d+(-D\\d+)?-(S|C)", TAG)) %>%
    select_(~SCAN_ID, ~BARCODE,  ~RESEARCH_PROJECT, ~DATE_CREATED, ~DRUGSET_ID,
            ~CELL_LINE_NAME, ~ends_with("CELL_ID"),  ~COSMIC_ID, ~POSITION,
            ~TAG, ~DRUG_ID, ~CONC, ~INTENSITY) %>%
    mutate_(lib_drug = ~sub("((L|R)\\d+)(-D\\d+)?-(S|C)", "\\1", TAG),
            lib_drug = ~ifelse(grepl("^A.+", lib_drug), yes = NA, no = lib_drug),
            anchor = ~sub("(A\\d+)-(S|C)", "\\1", TAG),
            anchor = ~ifelse(grepl("^(L|R).+", anchor), yes = NA, no = anchor),
            dose = ~sub("(A|L|R)\\d+-?(D\\d+)?-(S|C)", "\\2", TAG),
            treatment = ~sub("((A|L|R)\\d+)(-D\\d+)?-(S|C)", "\\4", TAG)
    ) %>%
    select_(~-TAG) 
  
  libraries <- normalized_data %>% 
    filter_(~!is.na(lib_drug)) %>% 
    select_(~-anchor, DRUG_ID_lib = ~DRUG_ID, CONC = ~CONC)
  
  anchors <- normalized_data %>% 
    filter_(~!is.na(anchor)) %>% 
    select_(~-lib_drug, ~-dose, DRUG_ID_anch = ~DRUG_ID, CONC_anch = ~CONC)
  if (nrow(anchors) > 0){
    suppressMessages(normalized_data <- full_join(libraries, anchors))
  }
  else {
    normalized_data <- libraries
  }
  
  normalized_data <- left_join(normalized_data, nc1, by=("SCAN_ID"))
  normalized_data <- left_join(normalized_data, pc1, by=("SCAN_ID"))
  normalized_data <- normalized_data %>%
    mutate_(normalized_intensity = ~((INTENSITY - PC) / (NC - PC)))
  
  if(trim){
    normalized_data <- normalized_data %>%
      mutate_(normalized_intensity = 
                ~(ifelse(normalized_intensity > 1, 1, normalized_intensity))) %>%
      mutate_(normalized_intensity =
                ~(ifelse(normalized_intensity < 0, 0, normalized_intensity)))
  }
  
  normalized_data <- normalized_data %>% 
    mutate_(norm_neg_pos = ~paste(neg_control, pos_control, sep = "+"))
  
  normalized_data <- normalized_data %>% 
    mutate_(time_stamp = ~ Sys.time())
  
  return(normalized_data)
}


main_normalizer = function(gdsc_dataset){
  
  gdsc_dataset = removeFailedDrugs(gdsc_dataset) %>% removeMissingDrugs() %>% normalizeData()
  
  return(gdsc_dataset)
}

add_dmso = function(dataset){
  
  # Identify the unique cell_line-drug pairs
  unique_pairs <- dataset %>%
    select(cell_line, drug) %>%
    distinct()
  
  # Add the new rows with 0 dose and 1.0 response for each unique cell_line-drug pair
  new_rows <- unique_pairs %>%
    mutate(dose = 0, response = 1.0)
  
  # Combine the original dataset with the new rows
  final_dataset <- bind_rows(dataset, new_rows)
  
}

drug_name_mapper = function(dataset, drug_df){
  
  drug_df <- drug_df %>%
    select(DRUG_ID, DRUG_NAME) %>%
    rename(drug = DRUG_ID)
  
  merged_dataframe <- dataset %>%
    left_join(drug_df, by = "drug") %>%
    mutate(DRUG_NAME = ifelse(is.na(DRUG_NAME), as.character(drug), DRUG_NAME)) %>%
    select(-drug) %>%
    rename(drug = DRUG_NAME)
  
  return(merged_dataframe)
}

transform_for_curvecurator = function(df){
  
  # Remove references
  df = df[grep("^L", df$lib_drug), ]
  
  
  # Create new cell_name column, if there is more than 1 CELL_ID for that cell_line
  # it merges cell_name with cell_id
  
  # df <- df %>%
  #   group_by(CELL_LINE_NAME) %>%
  #   mutate(
  #     count_cell_id = n_distinct(CELL_ID)
  #   ) %>%
  #   ungroup() %>%
  #   mutate(
  #     CELL_NAME = ifelse(
  #       count_cell_id > 1,
  #       paste(CELL_LINE_NAME, CELL_ID, sep = "_"),
  #       as.character(CELL_LINE_NAME)
  #     )
  #   ) %>%
  #   select(-count_cell_id)
  
  

  # Extract related colunms and rename
  df = df %>% select(CELL_LINE_NAME, DRUG_ID_lib, CONC, normalized_intensity) %>%
              rename(cell_line = `CELL_LINE_NAME`,
                     drug = `DRUG_ID_lib`,
                     dose = `CONC`,
                     response = `normalized_intensity`)
  
  # Add DMSO for each cell-drug pair
  df = add_dmso(df)
  
  # Add map drug IDs to drug_names
  df = drug_name_mapper(df, drug_names)
  
  # Order the colunms
  col_order = c("cell_line", "drug", "dose", "response")
  
  df = df[, col_order] %>% arrange(cell_line, drug, dose) 

  
  return(df)
  
}

## GDSC1 

gdsc1_df = read.csv("Raw_data/GDSC1_public_raw_data_24Jul22.csv/GDSC1_public_raw_data_24Jul22.csv")

drug_names = read.csv("meta/screened_compounds_rel_8.4.csv")

gdsc1_norm = main_normalizer(gdsc1_df)

formatted_gdsc1 = transform_for_curvecurator(gdsc1_norm)

## GDSC2 


gdsc2_df = read.csv("Raw_data/GDSC2_public_raw_data_24Jul22.csv/GDSC2_public_raw_data_24Jul22.csv")

drug_names = read.csv("meta/screened_compounds_rel_8.4.csv")

gdsc2_norm = main_normalizer(gdsc2_df)

formatted_gdsc2 = transform_for_curvecurator(gdsc2_norm)

######

mapped = map_cell_line_names(formatted_gdsc1)

mapped_gdsc1 = mapped[[1]]

mapped_gdsc1$drug = sanitize_name(mapped_gdsc1$drug)

write.table(mapped_gdsc1, file = "processed_data/gdsc1_normalized_mapped.tsv", sep = "\t", row.names = FALSE)


