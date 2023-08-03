library(tidyverse)

cell_meta = read.delim("meta/CCLE_sample_info_file_2012-10-18.txt")

new_df = data.frame("cell_line" = cell_meta$Cell.line.primary.name)

## Map cells to CVCL

mapped_cells = map_cell_line_names(new_df)

mapping_vector = mapped_cells[[2]]

new_df$sanitazed =  sapply(new_df$cell_line, sanitize_name, USE.NAMES = FALSE) 


new_df$mapped = mapping_vector[new_df$sanitazed]

new_df = new_df[,-2]

colnames(new_df) = c("Cell.line.primary.name", "CVCL_ID")

mapped_meta = cell_meta %>% left_join(new_df, by = "Cell.line.primary.name")

filt_meta = mapped_meta %>% 
  select(CVCL_ID, Site.Primary, Histology, Hist.Subtype1) %>%
  rename("cell_line" = CVCL_ID) 

filt_meta$cell_line = tolower(filt_meta$cell_line)

## 

cpds = read.csv("meta/comp_target.csv")
head(cpds)

cpds = cpds %>%
  select(Name, Action, Targeted.process.pathway) %>%
  rename("drug" = Name)
cpds$drug = tolower(cpds$drug)

abc_seed = read.csv("meta/CVCL_GDSC_seeds.csv")

##
unique(cpds$drug)
unique(sanitize_name(cpds$drug))
unique(abc_seed$drug) %in% unique(cpds$drug)
##
merged = abc_seed %>%
  left_join(filt_meta, by = "cell_line") 

merged = merged %>%
  left_join(cpds, by = 'drug')


write.csv(merged, "cell_drug_ABC_with_tissue.csv", row.names = FALSEmk)


