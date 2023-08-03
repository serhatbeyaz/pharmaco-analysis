
library(dplyr)

to_char_vector <- function(value) {
  value <- gsub("\\[|\\]|'", "", value)
  return(strsplit(value, ", ")[[1]])
}


# Helper function to sanitize cell line names
sanitize_name <- function(name) {
  # Convert to lower case
  name <- tolower(name)
  
  # Remove or replace special characters
  name <- gsub("[ _\\(\\)\\-]", "", name)  # removes spaces, hyphens, and underscores
    

  
  
  return(name)
}


# Map cell line names based on the mapping_table
map_cell_line_names <- function(dataset) {
  

  mapping_table = read.csv("cell_line_mapping2.csv")
  
  mapping_table$ID <- sapply(mapping_table$ID, sanitize_name)
  mapping_table$SY <- sapply(mapping_table$SY, sanitize_name)
  mapping_table$DepMap <- sapply(mapping_table$DepMap, sanitize_name)
  mapping_table$SY <- sapply(mapping_table$SY, to_char_vector)
  mapping_table$DepMap <- sapply(mapping_table$DepMap, to_char_vector)
  
  
  dataset <- as_tibble(dataset) # Convert dataset to tibble
  
  dataset$cell_line <- sapply(dataset$cell_line, sanitize_name, USE.NAMES = FALSE)
  
  initial_df = dataset
  
  dataset = dataset[!duplicated(dataset$cell_line), ]

  # Create a table that only includes ID and AC for the join operation
  mapping_id_ac <- mapping_table %>% select(AC, ID) %>% mutate(cell_line = tolower(ID))
  
  match_id <- dataset %>%
    left_join(mapping_id_ac, by = c(cell_line = "cell_line")) %>% mutate(mapped_cell_line = AC)
  
  no_match_id <- subset(match_id, is.na(mapped_cell_line))
  
  mapping_table_sy <- mapping_table %>% mutate(SY = strsplit(gsub("'", "", gsub("\"|\\[|\\]|\\s", "", SY)), ",")) %>% as_tibble()
  
  match_sy <- no_match_id %>%
    mutate(mapped_cell_line = purrr::map_chr(cell_line, function(x) {
      match_index <- which(purrr::map_lgl(mapping_table_sy$SY, function(y) x %in% y))
      if (length(match_index) > 0) {
        mapping_table_sy$AC[match_index][1]
      } else {
        NA_character_  # default value for character output
      }
    })) 
  
  final_df <- rbind(subset(match_id, !is.na(mapped_cell_line)), subset(match_sy, !is.na(mapped_cell_line)))
  
  final_df <- final_df[!duplicated(final_df$cell_line), ]
  
  # create a named vector where the names are the cell_line and the values are the mapped_cell_line
  conversion_vector <- setNames(final_df$mapped_cell_line, final_df$cell_line)
  
  # use this vector to replace the cell_line in df2
  initial_df$cell_line <- conversion_vector[initial_df$cell_line]
  
  before_ncell = length(unique(initial_df$cell_line))
  
  initial_df = initial_df[!is.na(initial_df$cell_line), ]
  
  after_ncell = length(unique(initial_df$cell_line))
  
  cat(after_ncell, "out of ", before_ncell, "cells are mapped!")
  
  
  return(list(initial_df, conversion_vector))

  


}


ccle = read.delim("ccle_normalized.tsv")
gdsc = read.delim("gdsc2_normalized.tsv")


### Drug map #####
library(dplyr)

drug_maps = read.csv("overlapping_compounds.csv", allowEscapes = FALSE)
drug_maps$id = drug_maps$id %>%tolower()
drug_maps$Synonyms = drug_maps$Synonyms %>% tolower()

ccle$drug = ccle$drug %>% tolower()
gdsc$drug = gdsc$drug %>% tolower()

gdsc_filtered = gdsc %>% filter(drug %in% drug_maps$id)


ccle_filtered = ccle %>% 
  left_join(drug_maps, by = c("drug" = "Synonyms")) %>% 
  mutate(drug = if_else(is.na(id), drug, id)) %>% 
  filter(drug %in% drug_maps$id) %>%
  select(-id)

### Map cells ### 

ccle_mapped = map_cell_line_names(ccle_filtered)

gdsc_mapped = map_cell_line_names(gdsc_filtered)

#### FIND OVERLAPPING CELL LINES
ccle_overlaps <- ccle_mapped %>%
  semi_join(gdsc_mapped, by = c("cell_line", "drug"))

gdsc2_overlaps <- gdsc_mapped %>%
  semi_join(ccle_mapped, by = c("cell_line", "drug"))


### Write table


write.table(ccle_overlaps, file = "ccle_normalized_overlapping.tsv", sep = "\t", row.names = FALSE)

write.table(gdsc2_overlaps, file = "gdsc2_normalized_overlapping.tsv", sep = "\t", row.names = FALSE)


