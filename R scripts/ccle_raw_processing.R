
getCCLErawData <-
  function(path.data=file.path("data", "CCLE"), result.type=c("array", "list")){
    
    
    ccle.raw.drug.sensitivity <- read.csv("CCLE_NP24.2009_Drug_data_2015.02.24.csv", stringsAsFactors=FALSE)
    ccle.raw.drug.sensitivity.list <- do.call(c, apply(ccle.raw.drug.sensitivity, 1, list))
    
    concentrations.no <- max(unlist(lapply(ccle.raw.drug.sensitivity[ , "Doses..uM."], function(x){length(unlist(strsplit(x, split = ",")))})))
    
    if(result.type == "array"){
      ## create the ccle.drug.response object including information viablilities and concentrations for each cell/drug pair
      obj <- array(NA, dim=c(length(unique(ccle.raw.drug.sensitivity[ , "Primary.Cell.Line.Name"])), length(unique(ccle.raw.drug.sensitivity[ , "Compound"])), 2, concentrations.no), dimnames=list(unique(ccle.raw.drug.sensitivity[ , "Primary.Cell.Line.Name"]), unique(ccle.raw.drug.sensitivity[ , "Compound"]), c("concentration", "viability"), 1:concentrations.no))
    }
    fnexperiment <- 
      function(values)  {
        cellline <- values["Primary.Cell.Line.Name"]
        drug <- values["Compound"]
        #doses <- as.numeric(unlist(strsplit(input.matrix["Doses (uM)"], split=", "))) #nature paper raw data
        doses <- as.numeric(unlist(strsplit(values["Doses..uM."], split=","))) # micro molar
        if(concentrations.no > length(doses)) {doses <- c(doses, rep(NA, concentrations.no - length(doses)))}
        
        #responses <- as.numeric(unlist(strsplit(input.matrix["Activity Data\n(raw median data)"], split=",")))  #nature paper raw data
        responses <- as.numeric(unlist(strsplit(values["Activity.Data..median."], split=","))) + 100
        if(concentrations.no > length(responses)) {responses <- c(responses, rep(NA, concentrations.no - length(responses)))}
        
        if(result.type == "array"){
          obj[cellline,drug, "concentration", 1:length(doses)] <<- doses
          obj[cellline,drug, "viability", 1:length(responses)] <<- responses
        }else{
          return(list(cell=cellline, drug=drug, doses=doses, responses=responses))#paste(doses, collapse = ","), responses=paste(responses, collapse = ",")))
        }
      }
    
    ccle.raw.drug.sensitivity.list <- do.call(c, apply(ccle.raw.drug.sensitivity, 1, list))
    ccle.raw.drug.sensitivity.res <- mapply(fnexperiment, values=ccle.raw.drug.sensitivity.list)
    if(result.type == "array"){
      return(list("data"=obj, "concentrations.no"=concentrations.no))
    }else{
      return(list("data"=ccle.raw.drug.sensitivity.res, "concentrations.no"=concentrations.no))
    }
  }

raw.sensitivity <- getCCLErawData(result.type="list")
con_tested <- raw.sensitivity$concentrations.no
raw.sensitivity <- t(raw.sensitivity$data)
raw.sensitivity <- t(apply(raw.sensitivity,1, function(x){unlist(x)}))

## manual curation of drug names
##########################################################################
#raw.sensitivity <- read.csv(file.path(inst("PharmacoGx"), "extdata", "ccle_sensitivity_detail.csv"))
#raw.sensitivity[raw.sensitivity[ ,2]=="PF2341066",2] <- "CRIZOTINIB"
raw.sensitivity[raw.sensitivity[ ,2]=="ZD-6474",2] <- "Vandetanib"
raw.sensitivity[raw.sensitivity[ ,2]=="PF2341066",2] <- "PF-2341066"
##########################################################################

#raw.sensitivity[ ,2] <- gsub(pattern=badchars, replacement="", raw.sensitivity[ ,2])
#raw.sensitivity[ ,2] <- paste("drugid", toupper(raw.sensitivity[ ,2]), sep="_")

rownames(raw.sensitivity)  <- sprintf("drugid_%s_%s",as.character(raw.sensitivity[ ,2]),as.character(raw.sensitivity[ ,1]))
raw.sensitivity <- raw.sensitivity[ ,-c(1,2)]
raw.sensitivity <- array(c(as.matrix(as.numeric(raw.sensitivity[ ,1:con_tested])), as.matrix(as.numeric(raw.sensitivity[ ,(con_tested+1):(2*con_tested)]))), c(nrow(raw.sensitivity), con_tested, 2),
                         dimnames=list(rownames(raw.sensitivity), colnames(raw.sensitivity[ ,1:con_tested]), c("Dose", "Viability")))

head(raw.sensitivity)

save(raw.sensitivity, con_tested, file="./drug_norm_post.RData")


raw.sensitivity.x <- parallel::splitIndices(nrow(raw.sensitivity), floor(nrow(raw.sensitivity)/1000))

dir.create("./slices/")

for(i in seq_along(raw.sensitivity.x)){
  
  slce <- raw.sensitivity[raw.sensitivity.x[[i]],,]
  saveRDS(slce, file=paste0("./slices/ccle_raw_sens_", i, ".rds"))
  
}

slices = list.files("slices/", pattern = "*.rds", full.names = T)

ccle_raw = data.frame()

for (i in 1:length(slices)){
  
  slice = readRDS(slices[i])
  
  slice_df = as.data.frame(slice)
  
  ccle_raw = rbind(ccle_raw, slice_df)
  
}


ccle_df <- ccle_df %>%
  rownames_to_column(var = "rowname") %>%
  separate(rowname, c("label", "drug", "cell_line"), sep = "_") %>%
  select(-label) %>%
  select(cell_line, drug, everything()) 

reshaped_data <- ccle_df %>%
  pivot_longer(cols = -c(cell_line, drug),
               names_to = c("dose_num", "variable"),
               names_sep = "\\.") %>%
  pivot_wider(names_from = variable, values_from = value) %>%
  rename_with(~str_replace_all(.,"doses", ""), starts_with("doses")) %>%
  as.data.frame()



renamed_data = reshaped_data %>%
  rename(dose = Dose,
         response = Viability) %>% 
  select(-dose_num)


summary(renamed_data)

ggplot(renamed_data, aes(x = response)) +
  geom_histogram(binwidth = 0.1, fill = "steelblue", color = "black") +
  theme_minimal() +
  labs(title = "Distribution of Response Variable",
       x = "Response",
       y = "Frequency")

renamed_data$response = renamed_data$response/100

final_df = add_dmso(renamed_data) %>% arrange(cell_line, drug, dose)

final_df = final_df[!is.na(final_df$dose), ]

write.table(final_df, file = "ccle_normalized.tsv", sep = "\t", row.names = FALSE)

