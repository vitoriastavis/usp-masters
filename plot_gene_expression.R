##################################
library(dplyr)
library(purrr)
library(ggplot2)
library(tidyr)

# Função auxiliar: encontra a linha do gene em um data.frame de DEGs
get_gene_row <- function(df, gene) {
  # 1) Se o gene estiver nos rownames
  if (!is.null(rownames(df)) && gene %in% rownames(df)) {
    out <- df[gene, , drop = FALSE]
    out$gene <- gene
    return(out)
  }
  
  # 2) Se houver alguma coluna de gene típica
  candidate_cols <- c("gene", "gene_symbol", "hgnc_symbol", "external_gene_name")
  col_gene <- intersect(candidate_cols, colnames(df))
  
  if (length(col_gene) > 0) {
    col_gene <- col_gene[1]
    df_subset <- df %>% filter(.data[[col_gene]] == gene)
    if (nrow(df_subset) > 0) {
      return(df_subset)
    }
  }
  
  # Se não encontrar, retorna NULL
  return(NULL)
}

# Extrai informações de um gene ao longo dos tipos celulares de um dataset
extract_gene_across_celltypes <- function(deg_list, gene) {
  # deg_list: lista tipo GSE213982_celltype
  # gene: string, ex. "SLC1A2"
  
  map_dfr(
    .x = names(deg_list),
    .f = function(ct) {
      df <- deg_list[[ct]]
      row_gene <- get_gene_row(df, gene)
      
      if (is.null(row_gene)) {
        return(NULL)  # gene não presente nesse tipo celular
      }
      
      # Garante que tenhamos as colunas esperadas (se não existirem, coloca NA)
      tibble(
        celltype   = ct,
        gene       = gene,
        avg_log2FC = row_gene$avg_log2FC %||% NA_real_,
        p_val      = row_gene$p_val %||% NA_real_,
        p_val_adj  = row_gene$p_val_adj %||% NA_real_,
        pct.1      = row_gene$pct.1 %||% NA_real_,
        pct.2      = row_gene$pct.2 %||% NA_real_
      )
    }
  )
}

# Operador "ou" seguro (caso a coluna não exista)
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

# Função de plot para um dataset (uma lista)
plot_gene_for_dataset <- function(deg_list, gene, dataset_name = deparse(substitute(deg_list))) {
  df_gene <- extract_gene_across_celltypes(deg_list, gene)
  
  if (nrow(df_gene) == 0) {
    stop(sprintf("Gene '%s' não encontrado em nenhum tipo celular de %s.", gene, dataset_name))
  }
  
  df_gene <- df_gene %>%
    mutate(
      sig = ifelse(!is.na(p_val_adj) & p_val_adj < 0.05, "FDR < 0.05", "ns"),
      celltype = factor(celltype, levels = names(deg_list))  # mantém ordem da lista
    )
  
  ggplot(df_gene, aes(x = celltype, y = avg_log2FC, fill = sig)) +
    geom_col() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(
      title = paste0("Expressão diferencial de ", gene, " em ", dataset_name),
      subtitle = "avg_log2FC por tipo celular (cores indicam significância FDR)",
      x = "Tipo celular",
      y = "avg_log2FC"
    ) +
    scale_fill_manual(values = c("FDR < 0.05" = "red", "ns" = "grey70")) +
    theme_minimal(base_size = 12)
}

# Exemplo de uso:
# 1) Para o primeiro dataset
p212 <- plot_gene_for_dataset(GSE213982_celltype, gene = "RARB", dataset_name = "GSE213982")
p144 <- plot_gene_for_dataset(GSE144136_celltype, gene = "RARB", dataset_name = "GSE213982")
# 2) Para o segundo dataset (troque pelo nome real, ex.: GSE101521_celltype)
# plot_gene_for_dataset(GSE101521_celltype, gene = "SLC1A2", dataset_name = "GSE101521")
