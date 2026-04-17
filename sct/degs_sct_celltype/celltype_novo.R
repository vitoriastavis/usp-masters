# Load necessary libraries
library(dplyr)
library(Seurat)
library(Matrix)

load("sct/degs_sct_celltype/GSE213982_2_99_2000_celltype.RData")
load("sct/degs_sct_celltype/GSE144136_1_99_2500_celltype.RData")
load("sct/degs_sct_celltype/celltypes.RData")


pval_adj_threshold <- 0.05
logfc_threshold <- 0.6

# filter by p and FC
GSE213982_celltype_filtered <- list()
for (ct in cell_types) {
  markers <- GSE213982_celltype[[ct]]
  
  markers_filtered <- markers %>%
    filter(p_val_adj <= pval_adj_threshold, abs(avg_log2FC) >= logfc_threshold)
  
  GSE213982_celltype_filtered[[ct]] <- markers_filtered
}

#save(GSE213982_celltype_filtered,
#     file="degs_sct_celltype/GSE213982_celltype_f.RData")

# filter by p and FC
GSE144136_celltype_filtered <- list()
for (ct in cell_types) {
  markers <- GSE144136_celltype[[ct]]
  
  markers_filtered <- markers %>%
    filter(p_val_adj <= pval_adj_threshold, abs(avg_log2FC) >= logfc_threshold)
  
  GSE144136_celltype_filtered[[ct]] <- markers_filtered
}

#save(GSE144136_celltype_filtered,
#     file="degs_sct_celltype/GSE144136_celltype_f.RData")



cell_types <- names(GSE144136_celltype_filtered)
intersection_celltypes <- list()
for (ct in cell_types) {
  intersection_celltypes[[ct]] <- intersect(rownames(GSE144136_celltype_filtered[[ct]]),
                                            rownames(GSE213982_celltype_filtered[[ct]]))
}

# Pasta com os CSVs das interseções
csv_folder <- "./2-targetprediction/data/intersections"
csv_files <- list.files(csv_folder, pattern = "\\.csv$", full.names = TRUE)

# Lista para guardar os resultados
intersection_results_celltype <- list()
for (file in csv_files){
  for (ct in cell_types){
    
    # Prediction methods + cell type degs
    key <- paste0(basename(file), "_", ct)
    
    # Read prediction methods intersections
    df <- read.csv(file, stringsAsFactors = FALSE)
    vec <- as.character(df$target)
    
    # Intersect with DEGs
    inter <- intersect(vec, intersection_celltypes[[ct]])
    intersection_results_celltype[[key]] <- inter
  }
}

print(intersection_results_celltype)
save(intersection_results, file="degs_sct_celltype/intersection_celltype.RData")

all_degs <- c(intersection_results, intersection_results_celltype)
save(all_degs, file="intersection_all_degs.RData")
