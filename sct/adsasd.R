degs_GSE144136_expr <- read.csv("sct/degs_sct_novo/GSE144136_1_99_2500_results_adj.csv")
degs_GSE144136 <- degs_GSE144136_expr$X

degs_GSE213982_expr <- read.csv("sct/degs_sct_novo/GSE213982_2_99_2000_results_adj.csv")
degs_GSE213982 <- degs_GSE213982_expr$X

pval_adj_threshold <- 0.05
logfc_threshold <- 0.6

degs_GSE213982_f <- degs_GSE213982_expr %>%
  filter(p_val_adj <= pval_adj_threshold, abs(avg_log2FC) >= logfc_threshold)

degs_GSE144136_f <- degs_GSE144136_expr %>%
  filter(p_val_adj <= pval_adj_threshold, abs(avg_log2FC) >= logfc_threshold)

snrnaseq_intersection <- intersect(degs_GSE213982_f$X, degs_GSE144136_f$X)


# Lista fixa para intersectar
snrnaseq_intersection <- as.character(snrnaseq_intersection)  # garantir vetor de caracteres

# Lista para guardar os resultados
results_sct <- list()

# Ler cada arquivo CSV e fazer interseção com snrnaseq_intersection
csv_folder <- "./data/target_prediction/intersections"
csv_files <- list.files(csv_folder, pattern = "\\.csv$", full.names = TRUE)

for (file in csv_files) {
  # Read raw lines first
  raw_lines <- readLines(file, warn = FALSE)

  # Nome base do arquivo para usar como key
  key <- gsub("^intersection_|\\.csv$", "", basename(file))
  
  # Ler CSV
  # print(file)
  df <- read.csv(file, stringsAsFactors = FALSE)
  # print(file)
  # Garantir vetor
  vec <- as.character(df$target)
  
  # Interseção com snrnaseq_intersection
  inter <- intersect(vec, snrnaseq_intersection)
  
  # Salvar na lista
  results_sct[[key]] <- inter
}
print(results_sct)
save(results_sct, file="./sct/degs_sct_novo/results_sct.RData")
