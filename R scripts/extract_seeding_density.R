library(dplyr)

gdsc_raw = read.csv("Raw_data/GDSC2_public_raw_data_24Jul22.csv/GDSC2_public_raw_data_24Jul22.csv")

drug_names = read.csv("meta/screened_compounds_rel_8.4.csv")

gdsc1_norm = main_normalizer(gdsc_raw)

# Extract related colunms and rename
df = gdsc1_norm  %>% select(CELL_LINE_NAME, DRUG_ID_lib, CONC, normalized_intensity, SEEDING_DENSITY) %>%
  rename(cell_line = `CELL_LINE_NAME`,
         drug = `DRUG_ID_lib`,
         dose = `CONC`,
         response = `normalized_intensity`)




drug_map = drug_name_mapper(df, drug_names)


cell_mapped = map_cell_line_names(drug_map)

gdsc2_df = cell_mapped[[1]]

gdsc2_df$drug = sanitize_name(gdsc2_df$drug)


drug_maps = read.csv("overlapping_compounds.csv", allowEscapes = FALSE)
drug_maps$id = drug_maps$id %>%tolower()
drug_maps$Synonyms = drug_maps$Synonyms %>% tolower()

gdsc2_df$drug = gdsc2_df$drug %>% tolower()

gdsc2_filtered = gdsc2_df %>% filter(drug %in% drug_maps$id)

gdsc2_pairs = gdsc2_df %>% 
  select(drug, cell_line, SEEDING_DENSITY) %>%
  group_by(cell_line, drug) %>%
  distinct() %>%
  arrange(desc(SEEDING_DENSITY))

gdsc2_pairs


df_grouped <- gdsc2_pairs %>%
  group_by(cell_line, drug) %>%
  summarize(SEEDING_DENSITY = mean(SEEDING_DENSITY, na.rm = TRUE)) %>%
  arrange(desc(SEEDING_DENSITY))


write.csv(gdsc2_pairs, 'meta/gdsc1_seeding_density.csv', row.names = FALSE)

gdsc2_dens = read.csv("gdsc2_seeding_density_all.csv")

head(gdsc2_dens)

gdsc2_dens$drug = sanitize_name(gdsc2_dens$drug)

write.csv(gdsc2_dens, "meta/gdsc2_seeding_density.csv")
