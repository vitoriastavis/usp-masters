library(GENIE3)
library(dplyr)
library(Seurat)
library(Matrix)

gse_path <- "data/snrnaseq/GSE144136/"
genie_dir <- paste0(gse_path, "genie/")

if (!dir.exists(genie_dir)) {
  dir.create(genie_dir)
  message("directory created: ", genie_dir)
}

load(paste0(gse_path, "/GSE144136_counts.RData"))
cat("loaded data")

GSE144136 <- FindVariableFeatures(GSE144136,
                                  selection.method = "vst",
                                  nfeatures = 2000)
variable_genes <- VariableFeatures(GSE144136)
cat("found variable features")

GSE144136_matrix <- as.matrix(GetAssayData(GSE144136, assay = "RNA", slot = "data"))
GSE144136_matrix <- GSE144136_matrix[variable_genes, ]

chunk_size <- 50
gene_names <- rownames(GSE144136_matrix)
n_genes <- length(gene_names)

GSE144136 <- NULL
gc()
cat("GENIE3 analysis prepared, starting...")

weight_list <- list()

for (i_reg in seq(1, n_genes, by = chunk_size)) {
  regulators_chunk <- gene_names[i_reg:min(i_reg + chunk_size - 1, n_genes)]
  
  for (i_tgt in seq(1, n_genes, by = chunk_size)) {
    cat("\tprocessing regulators", i_reg, "to", min(i_reg + chunk_size - 1, n_genes),
        "and targets", i_tgt, "to", min(i_tgt + chunk_size - 1, n_genes), "\n")
    
    target_chunk <- gene_names[i_tgt:min(i_tgt + chunk_size - 1, n_genes)]
    
    weights_chunk <- GENIE3(GSE144136_matrix,
                            regulators = regulators_chunk,
                            targets = target_chunk,
                            nTrees = 500,
                            nCores = 1)
    
    weight_list[[paste(i_reg, i_tgt, sep = "_")]] <- weights_chunk
    
    rm(weights_chunk)
    gc()
  }
}

# Combinar resultados 
weight_matrix_combined <- Reduce("+", weight_list) / length(weight_list)
save(weight_matrix_combined, file=paste0(gse_path, "GSE144136_genie_matrix.RData"))

cat("weight matrix built, filtering...")

# Convert to data frame/edge list
edge_list <- as.data.frame(as.table(weight_matrix_combined))
colnames(edge_list) <- c("Regulator", "Target", "Importance")

# Remove zero weights (optional)
edge_list <- edge_list %>%
              filter(Importance > 0)

# Sort by importance
edge_list_sorted <- edge_list %>%
                      arrange(desc(Importance))

# Top N interactions
top_edges <- head(edge_list_sorted, 1000)

# Save results
write.csv(top_edges, paste0(gse_path, "top_genie3_interactions.csv"),
                            row.names = FALSE)

cat("GENIE3 analysis finished")