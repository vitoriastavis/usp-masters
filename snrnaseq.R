library(dplyr)
library(Seurat)
library(Signac)
library(GenomeInfoDb)
library(ggplot2)
library(patchwork)
library(Matrix)
library(BiocParallel)
library(variancePartition)

gse_path <- "data/snrnaseq/GSE144136/"
load(paste0(gse_path, "/GSE144136_counts.RData"))
cat("loaded data")

GSE144136 <- FindVariableFeatures(GSE144136,
                                  selection.method = "vst",
                                  nfeatures = 2000)

cat("found variable features")

# Prepare your data
GSE144136_matrix <- GetAssayData(GSE144136, assay = "RNA", slot = "data")

meta_data <- GSE144136@meta.data

meta_data$condition <- as.numeric(factor(meta_data$condition))
meta_data$cell_type <- as.numeric(factor(meta_data$cell_type))

GSE144136 <- NULL
gc()

# Create a formula describing the variables to partition
form <- ~ (1|orig.ident) + condition + cell_type

param <- SnowParam(workers = 4, type = "SOCK", progressbar = FALSE)

# Define chunk size
chunk_size <- 1000
gene_names <- rownames(GSE144136_matrix)
n_genes <- length(gene_names)

cat("variance analysis prepared, starting...")

# Create empty list to store results
varPart_list <- list()
gc()

# Loop over chunks
for (i in seq(1, n_genes, by = chunk_size)) {
  cat("\tprocessing genes", i, "to", min(i+chunk_size-1, n_genes), "\n")
  
  # Subset matrix
  gene_subset <- gene_names[i:min(i + chunk_size - 1, n_genes)]
  matrix_chunk <- GSE144136_matrix[gene_subset, ]
  
  # Run model
  varPart_chunk <- fitExtractVarPartModel(matrix_chunk, form,
                                          meta_data, 
                                          BPPARAM = SnowParam(1))
  
  # Store result
  varPart_list[[length(varPart_list) + 1]] <- varPart_chunk
  
  # Clean up
  rm(matrix_chunk)
  gc()
}


# Combine results
varPart_combined <- do.call(rbind, varPart_list)
varPart_combined$gene <- rownames(varPart_combined)

cat("variance analysis finished")

top_genes <- head(varPart_combined[order(varPart_combined$condition, decreasing = TRUE),
                                   c("condition","gene")], 100)

load(file='data/target_prediction/targets_intersection.RData')
intersect(targets_intersection, top_genes$gene)

save(top_genes, file=paste0(gse_path, "GSE144136_variance.RData"))