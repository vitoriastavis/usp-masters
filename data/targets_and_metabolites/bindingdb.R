library(dplyr)
library(stringr)
bindingdb_results <- read.csv("bindingdb_results.tsv", sep="\t")
bindingdb_results <- bindingdb_results %>%
          mutate(Ligand_ID = str_extract(BindingDB.Ligand.Name, "CHEMBL\\d+"))
