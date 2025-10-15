library(readxl)

get_metabolites <- function(analysis_targets) {
  
  load("data/target_prediction/stpnovo/target_metabolites_stp.RData")
  stp_mappings <- read.csv("data/target_prediction/stpnovo/stp_mappings.csv")
  
  sea_result <- read_xlsx("data/target_prediction/sea/sea-combined.xlsx")
  sea_mappings <- read.csv("data/target_prediction/sea/sea_mappings.csv")
  colnames(sea_result) <- c("Query_ID", "Target_ID", "Affinity_Threshold", "P_Value",
                            "Max_Tc", "Cut_Sum", "Z_Score", "Name", "Description", "Query_Smiles")
  sea_result <- sea_result[grepl("_HUMAN$", sea_result$Target_ID), ]
  target_metabolites_sea <- split(sea_result$Query_ID, sea_result$Name)
  
  chembl_mappings <- read.csv("data/target_prediction/chemblpred/chemblpred_mappings.csv")
  chembl_mappings <- chembl_mappings[!is.na(chembl_mappings$HGNC_ID),]
  load("data/target_prediction/chemblpred/target_metabolites_chemblpred.RData")
  target_metabolites_chemblpred_chemblid <- target_metabolites_chemblpred
  
  target_metabolites_chemblpred <- list()
  for (target in names(target_metabolites_chemblpred_chemblid)){
    hgnc_id <- chembl_mappings[chembl_mappings$Query == target, "HGNC_ID"]
    
    if (length(hgnc_id) == 1){
      target_metabolites_chemblpred[[hgnc_id]] <- target_metabolites_chemblpred_chemblid[target][[1]]
    }
    else if (length(hgnc_id) == 2){
      for (id in hgnc_id){
        target_metabolites_chemblpred[[id]] <- target_metabolites_chemblpred_chemblid[target][[1]]  
      }
    }
  }
  
  targets_metabolites <- list()
  for (target in analysis_targets){
    print(target)
    all_metabolites <- c()
    
    # Get SEA metabolites
    sea_metabolites <- target_metabolites_sea[[target]]
    sea_metabolites <- gsub("^cid", "", sea_metabolites)

    all_metabolites <- c(sea_metabolites,
                         target_metabolites_stp[[target]],
                         target_metabolites_chemblpred[[target]])
    
    print(length(all_metabolites))
    
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
    
    original_queries <- chembl_mappings[!is.na(chembl_mappings$HGNC_ID) & chembl_mappings$HGNC_ID == target, ]
    for (query in original_queries$Query) {
      all_metabolites <- c(all_metabolites, target_metabolites_chemblpred[[query]])
    }
    print(length(all_metabolites))
    
    targets_metabolites[[target]] <- unique(all_metabolites)
  }
  
  return(targets_metabolites)
}
