load("targets_metabolites.RData")

library(dplyr)

# Number of unique metabolites
length(names(metabolites_to_targets))

# Mean number of metabolites per target
cat(mean(sapply(targets_to_metabolites, length)))

# Distribution of metabolites per target
hist(sapply(targets_to_metabolites, length))

# Targets with more metabolites
sort(sapply(targets_to_metabolites, length), decreasing = TRUE)

# Distribution of targets per metabolite
hist(sapply(metabolites_to_targets, length))

# Number of targets per metabolite
ligand_table <- tibble(
  Target = rep(names(targets_to_metabolites), lengths(targets_to_metabolites)),
  Ligand = unlist(targets_to_metabolites)
)
ligand_counts <- ligand_table %>%
  distinct(Target, Ligand) %>%
  count(Ligand, name = "Target_Count") %>%
  arrange(desc(Target_Count))

# Get the targets for  metabolites with higher target_count 
top_metabolites <- ligand_counts %>%
  filter(Target_Count > 1) %>%
  pull(Ligand)
top_metabolites_targets <- metabolites_to_targets[top_metabolites]

###
targets_to_metabolites_filtered <- targets_to_metabolites
targets_to_metabolites_filtered$S1PR1 <- NULL
targets_to_metabolites_filtered$NR4A1 <- NULL
targets_to_metabolites_filtered$EDNRB<- NULL

# Number of targets per metabolite
ligand_table <- tibble(
  Target = rep(names(targets_to_metabolites_filtered),
               lengths(targets_to_metabolites_filtered)),
  Ligand = unlist(targets_to_metabolites_filtered)
)
ligand_counts <- ligand_table %>%
  distinct(Target, Ligand) %>%
  count(Ligand, name = "Target_Count") %>%
  arrange(desc(Target_Count))

# Get the targets for  metabolites with higher target_count 
top_metabolites <- ligand_counts %>%
  filter(Target_Count > 1) %>%
  pull(Ligand)
top_metabolites_targets <- metabolites_to_targets[top_metabolites]

table(unname(unlist(top_metabolites_targets)))

