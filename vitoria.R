library(dplyr)
library(httr2)
library(jsonlite)
library(ggplot2)
library(readxl)
library(rvest)
library(readr)
library(stringr)
# Get SMILES from Compound ID
# Param: CID 
# Return (character): SMILES for the compound
get_smiles_from_cid <- function(cid) {
  url <- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/", cid, "/property/CanonicalSMILES/JSON")
  
  response <- request(url) |> req_perform()
  
  # If there was a response, get the data
  if (resp_status(response) == 200) {
    data <- resp_body_json(response)
    if (!is.null(data$PropertyTable$Properties[[1]]$CanonicalSMILES)) {
      return(data$PropertyTable$Properties[[1]]$CanonicalSMILES)
    }
  }
  
  return(NA)
}

# Get SMILES from Compound Name
# Param: name 
# Return (character): SMILES for the compound
get_smiles_from_name <- function(name) {
  # Convert name to URL-encoded format
  encoded_name <- URLencode(name, reserved = TRUE)
  
  # Construct the API URL
  url <- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/", encoded_name, "/property/CanonicalSMILES/JSON")
  
  # Make the request
  response <- request(url) |> req_perform()
  
  # If there was a response, get the data
  if (resp_status(response) == 200) {
    data <- resp_body_json(response)
    if (!is.null(data$PropertyTable$Properties[[1]]$CanonicalSMILES)) {
      return(data$PropertyTable$Properties[[1]]$CanonicalSMILES)
    }
  }
  
  return(NA)
}

# Read the CSV file from gutMGene
gutmgene <- read.csv("data/gutmgene.csv")
cat("N_metabolites raw file:",
    length(gutmgene$Metabolite..ID.))

# # Extract the pubchem_id from the "Metabolite (ID)" column
# gutmgene <- gutmgene %>%
#   mutate(cid = as.numeric(sub(".*\\((\\d+)(?:,.*)?\\)","\\1",
#                                      `Metabolite..ID.`)))
# cat("There are", nrow(gutmgene%>%filter(is.na(cid))), "PubChem ID NAs")

# Extract the metabolite name from the "Metabolite (ID)" column
gutmgene <- gutmgene %>%
  mutate(met_name = str_remove(Metabolite..ID., "\\([^()]*,[^()]*\\)$")) %>%
  mutate(met_name = str_trim(met_name))

# Remove NAs
if (nrow(gutmgene%>%filter(is.na(met_name))) != 0){
  cat("Removing", nrow(gutmgene%>%filter(is.na(met_name))), "met_name NAs")
  gutmgene <- gutmgene %>%
    filter(!is.na(met_name))
} else{
  cat("There are", nrow(gutmgene%>%filter(is.na(met_name))), "met_name NAs")
}
  

# Remove duplicates
gutmgene <- gutmgene %>%
  distinct(met_name, .keep_all=TRUE)
cat("N_metabolites without duplicates:", length(gutmgene$met_name))

# Get SMILES

smiles <- sapply(gutmgene$cid, get_smiles)

# Create dataframes to store SMILES and IDs for target prediction
df_smiles_ids <- data.frame(SMILES=smiles,
                         pubchem_ids=met_vector <- paste0("cid", gutmgene$cid),
                         stringsAsFactors=FALSE)

df_smiles <- data.frame(SMILES=smiles,
                        stringsAsFactors=FALSE)

# Save the dataframes as text files
write.table(df_smiles_ids, file="smiles_ids.txt",
            row.names=FALSE, col.names=FALSE,
            quote=FALSE, sep=" ")
write.table(df_smiles, file="smiles.txt",
            row.names=FALSE, col.names=FALSE,
            quote=FALSE, sep=" ")

# Read SEA results
sea_result <- read.csv("sea-results.xls")

# Get only the human targets 
sea_result <- sea_result[grepl("_HUMAN$", sea_result$Target.ID), ]
sea_targets <- unique(sea_result$Target.ID)
cat("There are", length(sea_targets), "targets predicted with SEA")

# Get STP results
smiles_list <- c("CC(=O)Oc1ccccc1C(=O)O", "C1=CC=CC=C1")  # Your SMILES list

base_url <- "http://www.swisstargetprediction.ch/"
results <- list()

for (smiles in smiles_list) {
  print(paste("Processing:", smiles))
  
  # Submit the form
  response <- POST(base_url, body = list(smiles = smiles), encode = "form")

  if (status_code(response) != 200) {
    print(paste("Failed for", smiles))
    next
  }

  # Parse HTML response
  page <- read_html(response)

  # Extract potential CSV download link
  csv_link <- page %>% html_nodes("a") %>% html_attr("href") %>% 
              .[grepl(".*csv$", .)] %>% 
              {paste0(base_url, .)}

  if (length(csv_link) > 0) {
    csv_filename <- paste0("result_", gsub("[^A-Za-z0-9]", "_", smiles), ".csv")
    
    # Download CSV file properly
    csv_response <- GET(csv_link)
    writeBin(content(csv_response, "raw"), csv_filename)
    
    print(paste("Downloaded:", csv_filename))

    # Optional: Read CSV into R
    csv_data <- read_csv(csv_filename)
    csv_data$SMILES <- smiles
    results[[length(results) + 1]] <- csv_data
  } else {
    print(paste("CSV not found for", smiles))
  }
  
  Sys.sleep(5)  # Avoid getting blocked
}

# Combine all results into a single CSV
if (length(results) > 0) {
  final_df <- do.call(rbind, results)
  write_csv(final_df, "all_results.csv")
  print("Saved all_results.csv")
}

