library(dplyr)

### test

informer_ctrp1 = read.delim("Raw_data/CTRPv1.0_2013_pub_Cell_154_1151/v10.M1.informer_set.txt")

informer_ctrp2 = readxl::read_xlsx("Raw_data/CTRPv2.0_2015_ctd2_ExpandedDataset/CTRPv2.0._INFORMER_SET.xlsx")

head(informer_ctrp1)
head(informer_ctrp2)


intersect(informer_ctrp1$master_cpd_id, informer_ctrp2$master_cpd_id)


raw = read.delim("Raw_data/CTRPv1.0_2013_pub_Cell_154_1151/v10.D1.raw_viability_data.txt")

"Vehicle" %in% unique(raw$cpd_name)

## No vehicle is found in the raw data, so normalization is not possible.



## I will gather the percent viability directly from processed data


avg_per = read.delim("Raw_data/CTRPv1.0_2013_pub_Cell_154_1151/v10.D2.avg_pct_viability_data.txt")

avg_per = avg_per %>%
  rename("cell_line" = ccl_name,
         "drug" = cpd_name,
         "dose" = cpd_conc_umol,
         "response" = cpd_avg_pv) %>%
  select(cell_line, drug, dose, response)


merge_cpd_id = avg_per %>% left_join(informer_ctrp1, by = c("drug" = "cpd_name")) %>%
  select(cell_line, drug, dose, response, master_cpd_id, gene_symbol_of_protein_target, target_or_activity_of_compound)


merge_ctrp_id = ctrp2 %>% left_join(informer_ctrp2, by = c("drug" = "cpd_name")) 





per_dmso = add_dmso(avg_per) %>% arrange(cell_line, drug, dose)


mapped = map_cell_line_names(per_dmso)

mapped_ctr = mapped[[1]]

drug_conversion = readxl::read_xlsx("GDSC-CCLE-CTRP_conversion.xlsx", sheet = "Drugs")

head(drug_conversion)

merged = mapped_ctr %>% left_join(drug_conversion, by = c("drug" = "CTRP name"))

filtered_ctrp1 = merged %>%
  mutate(drug = if_else(!is.na(`GDSC name`), `GDSC name`, drug)) %>%
  select(cell_line, drug, dose, response) %>%
  mutate(drug = sanitize_name(drug))

write.table(filtered_ctrp1, "processed_data/ctrp1_normalized_mapped.tsv", sep = "\t", row.names = FALSE)

### 

