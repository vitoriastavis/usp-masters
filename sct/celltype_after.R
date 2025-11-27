# Load necessary libraries
library(dplyr)
library(Matrix)

load("sct/degs_sct_celltype/GSE213982_2_99_2000_celltype.RData")
load("sct/degs_sct_celltype/GSE144136_1_99_2500_celltype.RData")
load("sct/celltypes.RData")

# filter by p and FC
GSE213982_celltype_filtered <- list()
for (ct in cell_types) {
  markers <- GSE213982_celltype[[ct]]
  
  markers_filtered <- markers %>%
    filter(p_val_adj <= 0.05, abs(avg_log2FC) >= 0.6)

  GSE213982_celltype_filtered[[ct]] <- markers_filtered
}

save(GSE213982_celltype_filtered,
file="./sct/degs_sct_celltype/GSE213982_celltype_f.RData")

# filter by p and FC
GSE144136_celltype_filtered <- list()
for (ct in cell_types) {
  markers <- GSE144136_celltype[[ct]]
  
  markers_filtered <- markers %>%
    filter(p_val_adj <= 0.05, abs(avg_log2FC) >= 0.6)

  GSE144136_celltype_filtered[[ct]] <- markers_filtered
}

save(GSE144136_celltype_filtered,
file="./sct/degs_sct_celltype/GSE144136_celltype_f.RData")

# Pasta com os CSVs das interseções
csv_folder <- "./data/target_prediction/intersections"

cell_types <- names(GSE144136_celltype_filtered)
intersection_celltypes <- list()
for (ct in cell_types) {
  intersection_celltypes[[ct]] <- intersect(rownames(GSE144136_celltype_filtered[[ct]]),
                                            rownames(GSE213982_celltype_filtered[[ct]]))
}

# Lista para guardar os resultados
results_sct_celltype <- list()

# Ler cada arquivo CSV e fazer interseção com snrnaseq_intersection
csv_files <- list.files(csv_folder, pattern = "\\.csv$", full.names = TRUE)

for (file in csv_files){
  for (ct in cell_types){
    
    # Prediction methods + cell type degs
    key <- paste0(basename(file), "_", ct)
    
    # Read prediction methods intersections
    df <- read.csv(file, stringsAsFactors = FALSE)
    vec <- as.character(df$target)
    
    # Intersect with DEGs
    inter <- intersect(vec, intersection_celltypes[[ct]])
    
    results_sct_celltype[[key]] <- inter
  }
}

print(results_sct_celltype)
save(results_sct_celltype, file="./sct/degs_sct_celltype/results_sct_celltype.RData")


