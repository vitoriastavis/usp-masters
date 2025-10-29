library(readxl)

get_metabolites <- function(analysis_targets) {
  
  load("data/target_prediction/stpnovo/target_metabolites_stp.RData")
  stp_mappings <- read.csv("data/target_prediction/stpnovo/stp_mappings.csv")
  target_metabolites_stp <- target_metabolites_stp[names(target_metabolites_stp) %in% analysis_targets]
  
  sea_result <- read_xlsx("data/target_prediction/sea/sea-combined.xlsx")
  sea_mappings <- read.csv("data/target_prediction/sea/sea_mappings.csv")
  colnames(sea_result) <- c("Query_ID", "Target_ID", "Affinity_Threshold", "P_Value",
                            "Max_Tc", "Cut_Sum", "Z_Score", "Name", "Description", "Query_Smiles")
  sea_result <- sea_result[grepl("_HUMAN$", sea_result$Target_ID), ]
  sea_result <- sea_result[sea_result$Name %in% analysis_targets, ]
  target_metabolites_sea <- split(sea_result$Query_ID, sea_result$Name)
  
  multitask_mappings <- read.csv("data/target_prediction/multitask/multitask_mappings.csv")
  multitask_mappings <- multitask_mappings[!is.na(multitask_mappings$HGNC_ID),]
  load("data/target_prediction/multitask/target_metabolites_multitask.RData")
  target_metabolites_multitask <- target_metabolites_multitask[names(target_metabolites_multitask) %in% analysis_targets]
  
  ppb_mappings <- read.csv("data/target_prediction/ppb2/ppb_mappings.csv")
  ppb_mappings <- ppb_mappings[!is.na(ppb_mappings$HGNC_ID),]
  load("data/target_prediction/ppb2/target_metabolites_ppb.RData")
  target_metabolites_ppb <- target_metabolites_ppb[names(target_metabolites_ppb) %in% analysis_targets]
  
  targets_metabolites <- list()
  for (target in analysis_targets){
    # print(target)
    all_metabolites <- c()
    
    # Get SEA metabolites
    sea_metabolites <- target_metabolites_sea[[target]]
    sea_metabolites <- gsub("^cid", "", sea_metabolites)

    all_metabolites <- c(sea_metabolites,
                         target_metabolites_stp[[target]],
                         target_metabolites_multitask[[target]],
                         target_metabolites_ppb[[target]])
    
    original_queries <- sea_mappings[!is.na(sea_mappings$HGNC_ID) & sea_mappings$HGNC_ID == target, ]
    for (query in original_queries$Query) {
      metabolites <- target_metabolites_sea[[query]]
      metabolites <- gsub("^cid", "", metabolites)
      all_metabolites <- c(all_metabolites, metabolites)
    }
    
    # print(length(all_metabolites))

    original_queries <- stp_mappings[!is.na(sea_mappings$HGNC_ID) & stp_mappings$HGNC_ID == target, ]
    for (query in original_queries$Query) {
      all_metabolites <- c(all_metabolites, target_metabolites_stp[[query]])
    }
    
    # print(length(all_metabolites))

    original_queries <- multitask_mappings[!is.na(multitask_mappings$HGNC_ID) & multitask_mappings$HGNC_ID == target, ]
    for (query in original_queries$Query) {
      all_metabolites <- c(all_metabolites, target_metabolites_multitask[[query]])
    }
    
    # print(length(all_metabolites))
    
    original_queries <- ppb_mappings[!is.na(ppb_mappings$HGNC_ID) & ppb_mappings$HGNC_ID == target, ]
    for (query in original_queries$Query) {
      all_metabolites <- c(all_metabolites, target_metabolites_ppb[[query]])
    }
    
    # print(length(all_metabolites))
    
    targets_metabolites[[target]] <- unique(all_metabolites)
    
    # print(length(all_metabolites))
  }
  
  return(targets_metabolites)
}
