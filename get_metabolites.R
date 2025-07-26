library(readxl)

get_metabolites <- function(analysis_targets) {
  
  load("data/target_prediction/stpnovo/target_metabolites_stp.RData")
  
  stp_mappings <- read.csv("data/target_prediction/stpnovo/stp_mappings.csv")
  sea_result <- read_xlsx("data/target_prediction/sea/sea-combined.xlsx")
  sea_mappings <- read.csv("data/target_prediction/sea/sea_mappings.csv")
  
  # Rename columns to avoid special character issues
  colnames(sea_result) <- c("Query_ID", "Target_ID", "Affinity_Threshold", "P_Value",
                            "Max_Tc", "Cut_Sum", "Z_Score", "Name", "Description", "Query_Smiles")
  
  # Get only the human targets 
  sea_result <- sea_result[grepl("_HUMAN$", sea_result$Target_ID), ]
  target_metabolites_sea <- split(sea_result$Query_ID, sea_result$Name)
  
  targets_metabolites <- list()
  
  for (target in analysis_targets) {
    all_metabolites <- c()
    
    # Get SEA metabolites
    metabolites <- target_metabolites_sea[[target]]
    metabolites <- gsub("^cid", "", metabolites)
    all_metabolites <- c(all_metabolites, metabolites)
    
    # Get STP metabolites
    all_metabolites <- c(all_metabolites, target_metabolites_stp[[target]])
    
    original_queries <- sea_mappings[!is.na(sea_mappings$HGNC_ID) & sea_mappings$HGNC_ID == target, ]
    for (query in original_queries$Query) {
      metabolites <- target_metabolites_sea[[query]]
      metabolites <- gsub("^cid", "", metabolites)
      all_metabolites <- c(all_metabolites, metabolites)
    }
    
    original_queries <- stp_mappings[!is.na(sea_mappings$HGNC_ID) & stp_mappings$HGNC_ID == target, ]
    for (query in original_queries$Query) {
      all_metabolites <- c(all_metabolites, target_metabolites_stp[[query]])
    }
    
    targets_metabolites[[target]] <- unique(all_metabolites)
  }
  
  return(targets_metabolites)
}
