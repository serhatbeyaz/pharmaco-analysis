library(dplyr)

df = read.delim("Raw_data/CTRPv2.0_2015_ctd2_ExpandedDataset/v20.data.curves_post_qc.txt")
meta_exp = read.delim("Raw_data/CTRPv2.0_2015_ctd2_ExpandedDataset/v20.meta.per_experiment.txt")
meta_cells = read.delim("Raw_data/CTRPv2.0_2015_ctd2_ExpandedDataset/v20.meta.per_cell_line.txt") %>%
  select(master_ccl_id, ccl_name, ccle_primary_site, ccle_primary_hist)

meta_cpd = read.delim("Raw_data/CTRPv2.0_2015_ctd2_ExpandedDataset/v20.meta.per_compound.txt") %>%
  select(master_cpd_id, cpd_name)

meta_exp = meta_exp[!duplicated(meta_exp$experiment_id), ]

## Scale auc values between 0-1

df = df %>% filter(area_under_curve < 20)

scale01 <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}


df %>% group_by(drug) %>%

scaled_auc = scale01(df$area_under_curve)

rob_scale = robust_scalar(df$area_under_curve)

df$area_under_curve = 1 - scaled_auc




with_cell = df %>% left_join(meta_exp, by = "experiment_id")

with_cell = with_cell %>% left_join(meta_cells, by = "master_ccl_id")

with_cpd = with_cell %>% left_join(meta_cpd, by = "master_cpd_id")



head(df)

head(with_cpd)

filtered = with_cpd %>%
  select(ccl_name, cpd_name, area_under_curve) %>%
  rename(cell_line = ccl_name,
         drug = cpd_name,
         AUC = area_under_curve)


dups = filtered %>%
  group_by(drug, cell_line) %>%
  summarise(count = n())

dupss = dups %>% filter(count > 1)


scaled_df = filtered %>%
  group_by(drug) %>%
  mutate(AUC_scaled = 1 - scale01(AUC)) %>%
  ungroup()


head(scaled_df)
hist(scaled_df$AUC_scaled)

###





### Map cell lines and drugs


mapped_df = map_cell_line_names(scaled_df)

mapped_df = mapped_df[[1]]


### Map drugs

conversion = readxl::read_xlsx("GDSC-CCLE-CTRP_conversion.xlsx", sheet = "Drugs") %>% 
  rename(drug = `CTRP name`)

conversion = mutate_all(conversion[,3:7], function(x) na_if(x, "N/A"))

conversion


test = mapped_df %>%
  filter(drug %in% conversion$drug) %>%
  left_join(conversion, by = 'drug')

new_drug = test %>%
  mutate(
    drug = if_else(!is.na(`GDSC name`), `GDSC name`, `CCLE name`))


drug_mapped = new_drug %>% 
  select(cell_line, drug, AUC_scaled)

drug_mapped$drug = tolower(drug_mapped$drug)


unique(dups$drug)

### 

library(tidyr)

pivot_auc = drug_mapped %>%
  pivot_wider(names_from = "cell_line",
              values_from = "AUC_scaled",
              values_fn = list(AUC_scaled = mean))


head(pivot_auc)

hist(drug_mapped$AUC_scaled)

write.csv(x = pivot_auc, file = "Auc_values_published/AUC_values_ctrp2_published.csv")
###
