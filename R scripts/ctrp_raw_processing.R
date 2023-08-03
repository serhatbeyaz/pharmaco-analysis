
library(dplyr)


raw_df = read.delim("Raw_data/CTRPv2.0_2015_ctd2_ExpandedDataset/v20.data.per_cpd_well.txt")

meta_exp = read.delim("Raw_data/CTRPv2.0_2015_ctd2_ExpandedDataset/v20.meta.per_experiment.txt")
meta_cells = read.delim("Raw_data/CTRPv2.0_2015_ctd2_ExpandedDataset/v20.meta.per_cell_line.txt") %>%
  select(master_ccl_id, ccl_name, ccle_primary_site, ccle_primary_hist)

meta_cpd = read.delim("Raw_data/CTRPv2.0_2015_ctd2_ExpandedDataset/v20.meta.per_compound.txt") %>%
  select(master_cpd_id, cpd_name)


test_df = meta_exp %>% left_join(meta_cells, by = "master_ccl_id") %>% 
  select(experiment_id, baseline_signal, cells_per_well, growth_mode, master_ccl_id, ccl_name, ccle_primary_site, ccle_primary_hist)

nd_meta = test_df[!duplicated(test_df$experiment_id), ]

assay_df = read.delim("Raw_data/CTRPv2.0_2015_ctd2_ExpandedDataset/v20.meta.per_assay_plate.txt") %>%
  select(experiment_id, assay_plate_barcode, dmso_plate_avg_log2)

head(assay_df)

head(test_df)

viability_transform = function(raw, nc){
  
  raw = 2^raw
  nc = 2^nc
  
  viable_cells = (nc - raw)/nc
  
  return(1 - viable_cells)
  
  
}



raw_df = raw_df %>% left_join(meta_cpd, by = "master_cpd_id")


whole_df = raw_df %>% left_join(nd_meta, by = "experiment_id")

plate_df = whole_df %>% left_join(assay_df, by = c("assay_plate_barcode", "experiment_id"))


plate_df$viability = viability_transform(plate_df$raw_value_log2, plate_df$dmso_plate_avg_log2)


avg_response = plate_df %>% 
  group_by(cpd_name, ccl_name, cpd_conc_umol) %>%
  summarise(cell_line = ccl_name, 
            compound = cpd_name,
            dose = cpd_conc_umol,
            response = viability)

head(avg_response)

str(avg_response)

testt = avg_response %>% reframe(response = response)



final_df <- testt %>% 
  rename(cell_line = ccl_name, 
         drug = cpd_name, 
         dose = cpd_conc_umol)


head(final_df)
str(final_df)
### Map cell lines and drugs


mapped_df = map_cell_line_names(final_df)

mapped_df = mapped_df[[1]]


### Map drugs

conversion = readxl::read_xlsx("GDSC-CCLE-CTRP_conversion.xlsx", sheet = "Drugs") %>% 
  rename(drug = `CTRP name`)

conversion = mutate_all(conversion[,3:7], function(x) na_if(x, "N/A"))

conversion

filtered_ctrp2 = mapped_df %>%
  left_join(conversion, by = "drug") %>%
  mutate(drug = if_else(!is.na(`GDSC name`), `GDSC name`, drug)) %>%
  select(cell_line, drug, dose, response) %>%
  mutate(drug = sanitize_name(drug))


### 

cc_ready = add_dmso(filtered_ctrp2) %>% arrange(cell_line, drug, dose)

write.table(filtered_ctrp2, "processed_data/ctrp2_normalized_mapped.tsv", sep = "\t", row.names = FALSE)


####

post_qc = read.delim("Raw_data/CTRPv2.0_2015_ctd2_ExpandedDataset/v20.data.curves_post_qc.txt")

post_qc_join = post_qc %>% left_join(nd_meta, by = "experiment_id")
head(post_qc_join)

drug_cell = plate_df %>% select(master_ccl_id, master_cpd_id, cpd_conc_umol)
head(drug_cell)

max_min_conc = plate_df %>% group_by(master_ccl_id, master_cpd_id) %>%
  summarise(min_conc = min(cpd_conc_umol),
            max_conc = max(cpd_conc_umol))


aucs = post_qc_join %>% select(master_ccl_id, master_cpd_id, area_under_curve)

merged = aucs %>% left_join(max_min_conc, by = c("master_ccl_id", "master_cpd_id"))

test = merged %>%
  mutate(norm_auc = area_under_curve/(log10(max_conc*1e-6) - log10(min_conc*1e-6)))

hist(test$norm_auc)
summary(test$norm_auc)
hist(test$area_under_curve)

filtered_ctrp2 %>% filter(drug == "bi2536",
                          cell_line == "CVCL_1331") %>% as.data.frame()
