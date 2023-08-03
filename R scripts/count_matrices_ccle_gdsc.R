

library(PharmacoGx)
source("R scripts/mapper.R")


# Download the object for CCLE
CCLE <- downloadPSet('CCLE_2015')

GDSC2 <- downloadPSet('GDSC_2020(v2-8.2)')

saveRDS(CCLE, file = "CCLE_pSet_haibe.rds")
saveRDS(GDSC2, file = "GDSC2_pSet_haibe.rds")

# Extract the expression data to a matrix
CCLE.expression <- summarizeMolecularProfiles(CCLE, mDataType="rna")

CCLE.expression@assays@data$exprs[1:10, 1:10]

GDSC.expression <- summarizeMolecularProfiles(GDSC2, mDataType = "rna")

GDSC.expression@assays@data$exprs[1:10, 1:10]

#####

CCLE.cnv = summarizeMolecularProfiles(CCLE, mDataType = "cnv")

GDSC.cnv = summarizeMolecularProfiles(GDSC2, mDataType = "cnv")

CCLE.cnv@assays@data$exprs[1:10, 1:10]

GDSC.cnv@assays@data$exprs[1:10, 1:10]



# Map CCLE cell line names
count_ccle = CCLE.expression@assays@data$exprs

colnames(count_ccle) = sapply(colnames(count_ccle), sanitize_name, USE.NAMES = FALSE)

ccle_df = data.frame(cell_line = colnames(count_ccle))

mapped_ccle = map_cell_line_names(ccle_df)

conversion_vector = mapped_ccle[[2]] 

colnames(count_ccle) = conversion_vector[colnames(count_ccle)]
ncol(count_ccle)

count_ccle = count_ccle[,!is.na(colnames(count_ccle))]
ncol(count_ccle)

count_ccle[1:10,1:10]

# Map 
count_gdsc = GDSC.expression@assays@data$exprs

colnames(count_gdsc) = sapply(colnames(count_gdsc), sanitize_name, USE.NAMES = FALSE)

gdsc_df = data.frame(cell_line = colnames(count_gdsc))

mapped_gdsc = map_cell_line_names(gdsc_df)

conversion_vector_gdsc = mapped_gdsc[[2]] 

colnames(count_gdsc) = conversion_vector_gdsc[colnames(count_gdsc)]
ncol(count_gdsc)

count_gdsc = count_gdsc[,!is.na(colnames(count_gdsc))]
ncol(count_gdsc)

count_gdsc[1:10, 1:10]

### Overlap

common_cell_lines = intersect(colnames(count_gdsc), colnames(count_ccle))
common_transcript = intersect(rownames(count_gdsc), rownames(count_ccle))

## Filter overlaps

ccle_common = count_ccle[common_transcript, common_cell_lines]
gdsc_common = count_gdsc[common_transcript, common_cell_lines]

write.csv(count_ccle, "CCLe_expression_count-matrix.csv", row.names = TRUE)
write.csv(count_gdsc, "GDSC2_expression_count-matrix.csv", row.names = TRUE)

## Annotation
library(rtracklayer)
library(dplyr)

gtf_path <- "Homo_sapiens.GRCh38.109.gtf" 
gtf <- rtracklayer::import(gtf_path)

gtf_df = data.frame(gtf@elementMetadata)

annot <- gtf_df %>%
  dplyr::select(gene_id, gene_name) %>%
  group_by(gene_id, gene_name) %>%
  distinct()

annot_test = gtf_df %>%
  dplyr::select(gene_id, gene_name) %>%
  distinct()


##

annot_map = annot_test$gene_name
names(annot_map) = annot_test$gene_id
annot_map[is.na(annot_map)] <- names(annot_map)[is.na(annot_map)]


## Map annotations

rownames(ccle_common) = annot_map[rownames(ccle_common)]
rownames(gdsc_common) = annot_map[rownames(gdsc_common)]

ccle_common[1:6,1:6]

write.csv(ccle_common, "CCLE_expression_count-matrix.csv", row.names = TRUE)
write.csv(gdsc_common, "GDSC2_expression_count-matrix.csv", row.names = TRUE)


## Correlation

# Ensure the cell line order is the same in both matrices
ccle_common <- ccle_common[, colnames(gdsc_common)]

# Initialize a vector to store the correlation coefficients
correlations <- numeric(ncol(ccle_common))

# Calculate the correlation for each cell line
for (i in seq_len(ncol(ccle_common))) {
  correlations[i] <- cor(ccle_common[, i], gdsc_common[, i], method = "pearson")
}

# Create a dataframe with the results
results <- data.frame(Cell_Line = colnames(ccle_common), Correlation = correlations)


genes = c("TOP2A", "TOP1", "MLH1", "MSH2", "ERCC1", "MGMT", "GSTP1", "ABCB1", "XRCC1", "ERCC2", "IL10", "ABCG2")
