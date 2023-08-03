

source("R scripts/mapper.R")
library("tidyr")
library("dplyr")

ccle_pubs = read.csv("Auc_values_published/ccle_pub_raw.csv", header = TRUE)


length(colnames(ccle_pubs))

df_long <- ccle_pubs %>% 
  gather(key = "cell_line", value = "AUC", -X) %>% 
  rename(drug = X)

head(df_long)

long_ccle = df_long[df_long$drug == "Crizotinib",] %>% filter(!is.na(AUC))

mapped_ccle = map_cell_line_names(long_ccle)


gdsc_pub = read.csv("Auc_values_published/AUC_values_gdsc2_published.csv")
head(gdsc_pub)

gdsc_long = gdsc_pub %>%
  gather(key = "cell_line", value = "AUC", -drug)

gdsc_crizo = gdsc_long %>% filter(drug == "crizotinib",
                                  !is.na(AUC))


inner_join(mapped_ccle, gdsc_crizo, by = 'cell_line')

mapped_df = map_cell_line_names(df_long)



drug_maps = read.csv("overlapping_compounds.csv", allowEscapes = FALSE)
drug_maps$id = drug_maps$id %>%tolower()
drug_maps$Synonyms = drug_maps$Synonyms %>% tolower()

mapped_df$drug = mapped_df$drug %>% tolower()

ccle_filtered = mapped_df %>% 
  left_join(drug_maps, by = c("drug" = "Synonyms")) %>% 
  mutate(drug = if_else(is.na(id), drug, id)) %>% 
  filter(drug %in% drug_maps$id) %>%
  select(-id)


ccle_filtered

ccle_wide = ccle_filtered %>%
  pivot_wider(names_from = 'cell_line',
              values_from = 'AUC')


write.csv(ccle_wide, "AUC_values_ccle_published.csv")
