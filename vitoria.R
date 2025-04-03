library(dplyr)
library(httr2)
library(jsonlite)
library(ggplot2)
library(readxl)
library(rvest)
library(readr)
library(stringr)
library(tidyr)
library(scales)
library(httr)
library(xml2)

get_cid <- function(name) {
  if (is.na(name) || name == "") return(NA)  # Handle missing names
  
  encoded_name <- URLencode(name, reserved = TRUE)
  
  # Try getting CID using compound name endpoint
  url_cid <- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/", 
                    encoded_name, "/cids/JSON")
  response_cid <- tryCatch({
    req <- request(url_cid) |> req_perform()
    if (resp_status(req) == 200) {
      data <- resp_body_json(req)
      if (!is.null(data$IdentifierList$CID)) {
        return(data$IdentifierList$CID[[1]])  # Return first CID
      }
    }
    return(NA)
  }, error = function(e) return(NA))
  
  # If a CID was found, return it immediately
  if (!is.na(response_cid)) {
    return(response_cid)
  }
  
  Sys.sleep(0.5)  # Pause to respect API rate limits
  
  # Try getting SID using compound name endpoint
  url_sid <- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/name/", encoded_name, "/sids/JSON")
  response_sid <- tryCatch({
    req <- request(url_sid) |> req_perform()
    if (resp_status(req) == 200) {
      data <- resp_body_json(req)
      if (!is.null(data$IdentifierList$SID) && length(data$IdentifierList$SID) > 0) {
        data$IdentifierList$SID[[1]]  # Store the SID in a variable instead of returning it
      } else {
        NA
      }
    } else {
      NA
    }
  }, error = function(e) return(NA))
  
  Sys.sleep(0.5)  # Pause to respect API rate limits
  
  # Try getting CID using SID
  url_sid_cid <- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/sid/", 
                        response_sid, "/JSON")
  response_cid <- tryCatch({
    req <- request(url_sid_cid) |> req_perform()
    if (resp_status(req) == 200) {
      data <- resp_body_json(req)
      if (!is.null(data$PC_Substances[[1]]$compound[[2]]$id$id$cid) && length(data$PC_Substances[[1]]$compound[[2]]$id$id$cid) > 0) {
        # print("Found CID using SID")  
        return(data$PC_Substances[[1]]$compound[[2]]$id$id$cid[[1]])  # Return first CID
      } else{
        print(data)
      }
    }
    return(NA)
  }, error = function(e) return(NA))
}

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
  if (is.na(name) || name == "") {
    return(NA)
  }
  
  encoded_name <- URLencode(name, reserved = TRUE)
  url <- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/", encoded_name, "/property/CanonicalSMILES/JSON")
  
  response <- tryCatch(
    {
      req <- request(url) |> req_perform()
      if (resp_status(req) == 200) {
        data <- resp_body_json(req)
        if (!is.null(data$PropertyTable$Properties[[1]]$CanonicalSMILES)) {
          return(data$PropertyTable$Properties[[1]]$CanonicalSMILES)
        }
      }
      return(NA)
    },
    error = function(e) {
      return(NA)
    }
  )
  
  Sys.sleep(0.5)  # Prevent API rate limits
  return(response)
}

# Read the CSV file from gutMGene
gutmgene <- read.csv("data/gutmgene.csv", stringsAsFactors = FALSE)

# Rename columns to avoid special character issues
colnames(gutmgene) <- c("Host_Species", "Gut_Microbe_ID", "Rank", "Metabolite_ID", "Evidence_Type", "Evidence_Amount")

cat("N metabolites raw file:",
    length(gutmgene$Metabolite_ID))
cat("N microbes in raw file:",
    length(unique(gutmgene$Gut_Microbe_ID)))
# cat(length(gutmgene$Rank=='genus'))

get_cid_vec <- Vectorize(get_cid)
# Extract the metabolite name from the "Metabolite_ID" column
gutmgene <- gutmgene %>%
  mutate(
    Metabolite_Name = str_remove(Metabolite_ID, "\\([^()]*\\)$"), # Removes last set of parentheses and their content
    Metabolite_Name = str_trim(Metabolite_Name), # Trims any leading or trailing whitespace
    Gut_Microbe_Name = str_extract(Gut_Microbe_ID, "^[^(]+"), # Extracts text before '('
    Gut_Microbe_ID = as.integer(str_extract(Gut_Microbe_ID, "(?<=\\().*(?=\\))")) # Extracts numbers inside '()'
    # Metabolite_CID = get_cid_vec(Metabolite_Name) # Obtains CID for the metabolite name
  )


# Filter only genus-level microbes
gutmgene_filtered <- gutmgene[gutmgene$Rank == "genus", ]
# Count occurrences of each genus
genus_counts <- table(gutmgene_filtered$Gut_Microbe_Name)
# Filter to only include genera with frequency > 5
genus_to_keep <- names(genus_counts[genus_counts > 5])
gutmgene_filtered <- gutmgene_filtered[gutmgene_filtered$Gut_Microbe_Name %in% genus_to_keep, ]
max_count <- ceiling(max(table(gutmgene_filtered$Gut_Microbe_Name)) / 5) * 5
gutmgene_filtered_count <- gutmgene_filtered %>%
  count(Gut_Microbe_Name) %>%
  mutate(Gut_Microbe_Name = reorder(Gut_Microbe_Name, -n))

# Create the plot
p1 <- ggplot(gutmgene_filtered_count, aes(x = Gut_Microbe_Name, y = n)) +
  geom_bar(stat = "identity", fill = "gray30", width = 0.5) +
  theme_minimal() +  
  theme(
    axis.title.y = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.text.x = element_text(size = 16, angle = 60, hjust = 1),
    panel.grid = element_blank(),  # Remove grid lines
    axis.line.x = element_line(color = "black", linewidth = 0.5),  # Only X-axis line
    axis.line.y = element_line(color = "black", linewidth = 0.5),  # Only Y-axis line
    axis.ticks = element_line(color = "black", linewidth = 1)  # Ensure ticks are visible
  ) +
  labs(x = "Gut Microbe (Genus)", 
       y = "Count") +
  scale_y_continuous(breaks = seq(0, max(gutmgene_filtered_count$n), by = 5)) +  
  coord_cartesian(ylim = c(0, max(gutmgene_filtered_count$n)))


# Save plot as PNG
png("microbe_distribution.png", width = 350, height = 300)
print(p1)
dev.off()

evidence_amount_counts <- gutmgene %>%
  count(Evidence_Amount)  # Count each unique value

p2 <- ggplot(evidence_amount_counts, aes(x = factor(Evidence_Amount), y = n)) +  
  geom_col(fill = "gray30", width = 0.5) +  # Now using geom_col correctly
  geom_text(aes(label = n), vjust = -0.5, size = 5, color = "black") +  # Add text labels
  theme_minimal() +  
  theme(
    axis.title.y = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.text.x = element_text(size = 16, hjust = 1),
    panel.grid = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.line.y = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 1)
  ) +
  scale_y_continuous(breaks = pretty_breaks(n = 7)) +
  labs(x = "Evidence Amount", 
       y = "Count") +  
  coord_cartesian(ylim = c(0, max(evidence_amount_counts)*1.03))

png("evidence_amount.png", width = 350, height = 300)
print(p2)
dev.off()

# Split the 'Evidence' column by commas, then count occurrences of each type of evidence
evidence_count <- gutmgene %>%
  separate_rows(Evidence_Type, sep = ";") %>%  # Split by commas into separate rows
  count(Evidence_Type) %>%  # Count occurrences of each type of evidence
  filter(Evidence_Type %in% c("causally", "correlatively"))  # Only keep causally and correlatively

p3 <- ggplot(evidence_count, aes(x = "", y = n, fill = Evidence_Type)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +  # Converte o gráfico de barras para pizza
  theme_void() +  # Remove fundo e eixos
  scale_fill_manual(values = c("causally" = "gray30", "correlatively" = "gray50")) +  # Personaliza as cores
  geom_text(aes(label = paste(n, "(", round(n / sum(n) * 100, 1), "%)", sep = "")),  # Adiciona contagens e porcentagens
            position = position_stack(vjust = 0.5), size = 6, color = "white") +  # Posiciona os rótulos
  theme(
    plot.title = element_text(size = 14, hjust = 0),
    legend.text = element_text(size = 14),  # Ajusta o tamanho do texto da legenda
    legend.title = element_text(size = 14)  # Ajusta o tamanho do título da legenda (se houver)
  )

png("evidence_types.png", width = 350, height = 300)
print(p3)
dev.off()

# Remove NAs
if (nrow(gutmgene%>%filter(is.na(Metabolite_CID))) != 0){
  cat("Removing", nrow(gutmgene%>%filter(is.na(Metabolite_CID))), "Metabolite_CID NAs")
  gutmgene_na <- gutmgene %>%
    filter(!is.na(Metabolite_CID))
} else{
  cat("There are", nrow(gutmgene%>%filter(is.na(Metabolite_CID))), "Metabolite_CID NAs")
}

# Remove duplicates
gutmgene_clean <- gutmgene_na %>%
  distinct(Metabolite_CID, .keep_all=TRUE)
cat("N_metabolites without duplicates:", length(gutmgene_clean$Metabolite_CID))

smiles <- sapply(gutmgene_clean$Metabolite_CID, get_smiles_from_cid)

# smiles <- sapply(gutmgene$Metabolite_Name, get_smiles_from_name)

# Create dataframes to store SMILES and IDs for target prediction
df_smiles_ids <- data.frame(SMILES=smiles,
                         pubchem_ids=paste0("cid", gutmgene_clean$Metabolite_CID),
                         stringsAsFactors=FALSE)

# df_smiles <- data.frame(SMILES=smiles,
#                         stringsAsFactors=FALSE)

# Save the dataframes as text files
write.table(df_smiles_ids, file="smiles_ids.txt",
            row.names=FALSE, col.names=FALSE,
            quote=FALSE, sep=" ")
# write.table(df_smiles, file="smiles.txt",
#             row.names=FALSE, col.names=FALSE,
#             quote=FALSE, sep=" ")

# Read SEA results
sea_result <- read.csv("data/sea-results.xls")
colnames(sea_result) <- c("Query_ID", "Target_ID", "Affinity_Threshold", "P_Value",
                          "Max_Tc", "Cut_Sum", "Z_Score", "Name", "Description", "Query_Smiles")

# Get only the human targets 
sea_result <- sea_result[grepl("_HUMAN$", sea_result$Target_ID), ]
sea_targets <- sea_result[!duplicated(sea_result$Target_ID), ]

# sea_targets <- unique(sea_result$Target.ID)
cat("There are", length(sea_targets$Target_ID), "targets predicted with SEA")

sea_targets_filtered <- sea_targets[sea_targets$P_Value <= 0.01, ]
cat("There are", length(sea_targets_filtered$Target_ID), "targets predicted with SEA with p <= 0.01")


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

