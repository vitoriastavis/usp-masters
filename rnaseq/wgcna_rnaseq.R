#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(GEOquery)
  library(edgeR)
  library(limma)
  library(WGCNA)
  library(dynamicTreeCut)
})

options(stringsAsFactors = FALSE)
allowWGCNAThreads()

#------------------------------------------------------------
# Função: parse_geo_characteristics
# (extrai age, sex, pmi, brain_ph, rin das characteristics_ch1*)
#------------------------------------------------------------
parse_geo_characteristics <- function(pheno) {
  ph <- pheno %>%
    as.data.frame() %>%
    tibble::rownames_to_column("sample_id")
  
  long <- ph %>%
    tidyr::pivot_longer(
      cols = dplyr::starts_with("characteristics_ch1"),
      names_to  = "char_col",
      values_to = "char"
    ) %>%
    dplyr::filter(!is.na(char), char != "") %>%
    tidyr::separate(
      col  = char,
      into = c("key", "value"),
      sep  = ":\\s*",
      extra = "merge",
      fill  = "right"
    ) %>%
    dplyr::mutate(
      key   = tolower(trimws(key)),
      value = trimws(value)
    ) %>%
    dplyr::mutate(
      std_key = dplyr::case_when(
        grepl("age at death", key) ~ "age",
        grepl("\\bage\\b", key) & !grepl("onset|first episode", key) ~ "age",
        grepl("sex|gender", key) ~ "sex",
        grepl("post.?mortem|pmi", key) ~ "pmi",
        grepl("brain\\s*pH", key, ignore.case = TRUE) ~ "brain_ph",
        grepl("brain\\s*ph", key) ~ "brain_ph",
        grepl("rin|rna integrity", key) ~ "rin",
        TRUE ~ key
      )
    )
  
  wide_char <- long %>%
    dplyr::group_by(sample_id, std_key) %>%
    dplyr::summarise(value = value[1], .groups = "drop") %>%
    tidyr::pivot_wider(
      names_from  = std_key,
      values_from = value
    )
  
  ph_final <- ph %>%
    dplyr::left_join(wide_char, by = "sample_id") %>%
    tibble::column_to_rownames("sample_id")
  
  if ("age" %in% colnames(ph_final)) {
    ph_final$age <- gsub("[^0-9.]", "", ph_final$age)
    ph_final$age <- suppressWarnings(as.numeric(ph_final$age))
  }
  if ("pmi" %in% colnames(ph_final)) {
    ph_final$pmi <- gsub("[^0-9.]", "", ph_final$pmi)
    ph_final$pmi <- suppressWarnings(as.numeric(ph_final$pmi))
  }
  if ("sex" %in% colnames(ph_final)) {
    ph_final$sex <- toupper(trimws(ph_final$sex))
    ph_final$sex[ph_final$sex %in% c("MALE","M","1")]   <- "M"
    ph_final$sex[ph_final$sex %in% c("FEMALE","F","2")] <- "F"
    ph_final$sex <- factor(ph_final$sex)
  }
  if ("rin" %in% colnames(ph_final)) {
    ph_final$rin <- gsub("[^0-9.]", "", ph_final$rin)
    ph_final$rin <- suppressWarnings(as.numeric(ph_final$rin))
  }
  if ("brain_ph" %in% colnames(ph_final)) {
    ph_final$brain_ph <- gsub("[^0-9.]", "", ph_final$brain_ph)
    ph_final$brain_ph <- suppressWarnings(as.numeric(ph_final$brain_ph))
  }
  
  ph_final
}

#------------------------------------------------------------
# Função: prepare_wgcna_matrix
# (voom + remoção de covariáveis → matriz residual para WGCNA)
#------------------------------------------------------------
prepare_wgcna_matrix <- function(counts, meta, sample_col,
                                 covars = c("age","sex","pmi","rin","brain_ph")) {
  
  if (!sample_col %in% colnames(meta)) {
    stop("A coluna '", sample_col, "' não existe na tabela meta.")
  }
  
  meta[[sample_col]] <- make.names(as.character(meta[[sample_col]]))
  colnames(counts)   <- make.names(colnames(counts))
  
  meta   <- meta %>% dplyr::filter(.data[[sample_col]] %in% colnames(counts))
  counts <- counts[, meta[[sample_col]], drop = FALSE]
  
  y <- edgeR::DGEList(counts)
  y <- edgeR::calcNormFactors(y, method = "TMM")
  design0 <- model.matrix(~ 1, data = meta)
  v <- limma::voom(y, design0, plot = FALSE)
  expr <- v$E
  
  meta <- meta[match(colnames(expr), meta[[sample_col]]), , drop = FALSE]
  
  covars_present <- intersect(covars, colnames(meta))
  
  if (length(covars_present) > 0) {
    keep <- stats::complete.cases(meta[, covars_present, drop = FALSE])
    
    if (sum(keep) < 4) {
      stop("Menos de 4 amostras completas após filtrar covariáveis: ",
           paste(covars_present, collapse = ", "))
    }
    
    meta <- meta[keep, , drop = FALSE]
    expr <- expr[, keep, drop = FALSE]
    
    design_cov <- model.matrix(
      as.formula(paste("~", paste(covars_present, collapse = "+"))),
      data = meta
    )
    
    if (nrow(design_cov) != ncol(expr)) {
      stop("Dimensões incompatíveis: nrow(design_cov) = ", nrow(design_cov),
           " vs ncol(expr) = ", ncol(expr))
    }
    
    fit <- limma::lmFit(expr, design_cov)
    expr_resid <- limma::residuals.MArrayLM(fit, expr)
  } else {
    expr_resid <- expr
  }
  
  expr_resid
}

#------------------------------------------------------------
# Função: run_wgcna
#------------------------------------------------------------
run_wgcna <- function(expr_resid, meta, sample_col, trait_col = "group") {
  
  meta[[sample_col]] <- make.names(as.character(meta[[sample_col]]))
  colnames(expr_resid) <- make.names(colnames(expr_resid))
  meta <- meta[match(colnames(expr_resid), meta[[sample_col]]), , drop = FALSE]
  
  trait <- meta[[trait_col]]
  if (is.factor(trait)) trait <- as.character(trait)
  trait_num <- ifelse(trait %in% c("MDD","MDD-S","MDD-NS"), 1, 0)
  
  datExpr <- t(expr_resid)
  
  gsg <- goodSamplesGenes(datExpr, verbose = 0)
  if (!gsg$allOK) {
    datExpr   <- datExpr[gsg$goodSamples, gsg$goodGenes]
    trait_num <- trait_num[gsg$goodSamples]
  }
  
  powers <- 1:20
  sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 0)
  softPower <- ifelse(is.na(sft$powerEstimate), 6, sft$powerEstimate)
  
  adjacency <- adjacency(datExpr, power = softPower)
  TOM <- TOMsimilarity(adjacency)
  dissTOM <- 1 - TOM
  
  geneTree <- hclust(as.dist(dissTOM), method = "average")
  
  modules <- cutreeDynamic(
    dendro = geneTree, distM = dissTOM,
    deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = 30
  )
  moduleColors <- labels2colors(modules)
  
  MEs <- moduleEigengenes(datExpr, colors = moduleColors)$eigengenes
  moduleTraitCor <- cor(MEs, trait_num, use = "pairwise.complete.obs")
  moduleTraitP   <- corPvalueStudent(moduleTraitCor, nrow(datExpr))
  
  list(
    datExpr        = datExpr,
    moduleColors   = moduleColors,
    MEs            = MEs,
    moduleTraitCor = moduleTraitCor,
    moduleTraitP   = moduleTraitP,
    softPower      = softPower,
    geneTree       = geneTree,
    trait_num      = trait_num
  )
}

#------------------------------------------------------------
# Função: extract_module_genes
#------------------------------------------------------------
extract_module_genes <- function(expr_resid, moduleColors, MEs, trait_num, p_thr = 0.05) {
  datExpr <- t(expr_resid)
  
  moduleTraitCor <- cor(MEs, trait_num, use = "pairwise.complete.obs")
  moduleTraitP   <- corPvalueStudent(moduleTraitCor, nrow(datExpr))
  
  modules <- rownames(moduleTraitCor)
  sig_modules <- modules[moduleTraitP[,1] < p_thr]
  
  out <- list()
  if (length(sig_modules) == 0) return(out)
  
  for (mod in sig_modules) {
    genes_mod <- names(moduleColors)[moduleColors == mod]
    if (length(genes_mod) == 0) next
    
    # kME: correlação gene x eigengene do módulo
    kME <- cor(datExpr[, genes_mod, drop = FALSE], MEs[, mod], use = "p")
    
    out[[mod]] <- tibble(
      gene   = genes_mod,
      kME    = as.numeric(kME),
      module = mod
    ) %>% arrange(desc(kME))
  }
  
  out
}

#------------------------------------------------------------
# MAIN: WGCNA para GSE80655 DLPFC (Control vs MDD)
#------------------------------------------------------------

DATA_DIR <- "."
OUT_DIR  <- "wgcna"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# 1) Carregar contagens
load(file.path(DATA_DIR, "80655.RData"))  # deve conter 'cnt_80655'

# 2) GEO + pheno
gse_80655 <- getGEO("GSE80655", GSEMatrix = TRUE)
pheno_80655 <- Biobase::pData(gse_80655[[1]])

pheno_80655 <- pheno_80655[, !duplicated(colnames(pheno_80655))] %>%
  dplyr::mutate(
    group = dplyr::case_when(
      grepl("\\bControl\\b|non-psy|healthy", characteristics_ch1.3, ignore.case = TRUE) ~ "Control",
      grepl("\\bMajor Depression\\b|major\\s*depress", characteristics_ch1.3, ignore.case = TRUE) ~ "MDD",
      grepl("\\bSchizophrenia", characteristics_ch1.3, ignore.case = TRUE) ~ "SCZ",
      grepl("\\bipolar Disorder", characteristics_ch1.3, ignore.case = TRUE) ~ "BPD",
      TRUE ~ NA_character_
    ),
    brain_region = sub("^brain region:\\s*", "", characteristics_ch1.1)
  )

pheno_80655 <- parse_geo_characteristics(pheno_80655)

pheno_80655 <- pheno_80655 %>%
  dplyr::filter(
    brain_region == "DLPFC",
    group %in% c("Control","MDD"),
    description %in% colnames(cnt_80655)
  )

# 3) Matriz residual para WGCNA
expr_resid_80655_DLPFC <- prepare_wgcna_matrix(
  counts    = cnt_80655,
  meta      = pheno_80655,
  sample_col = "description",
  covars    = c("age","sex","pmi","rin","brain_ph")
)

saveRDS(expr_resid_80655_DLPFC,
        file.path(OUT_DIR, "expr_80655_DLPFC.rds"))

# 4) Rodar WGCNA
wgcna_out <- run_wgcna(
  expr_resid = expr_resid_80655_DLPFC,
  meta       = meta_dlpfc,
  sample_col = "description",
  trait_col  = "group"
)

saveRDS(wgcna_out,
        file.path(OUT_DIR, "wgcna_80655_DLPFC.rds"))

# 5) Resumo módulo x MDD
mtc <- wgcna_out$moduleTraitCor
mtp <- wgcna_out$moduleTraitP

df_mtc <- tibble(
  module = rownames(mtc),
  cor_with_MDD = mtc[,1],
  pvalue       = mtp[,1]
)

readr::write_csv(df_mtc,
                 file.path(OUT_DIR, "module_trait_80655_DLPFC.csv"))

# 6) Genes por módulo significativo
trait_vec <- meta_dlpfc$group
trait_num <- ifelse(trait_vec == "MDD", 1, 0)

module_genes <- extract_module_genes(
  expr_resid  = expr_resid_80655_DLPFC,
  moduleColors = wgcna_out$moduleColors,
  MEs         = wgcna_out$MEs,
  trait_num   = trait_num,
  p_thr       = 0.05
)

if (length(module_genes) > 0) {
  for (mod in names(module_genes)) {
    fn <- file.path(OUT_DIR, paste0("module_genes_", mod, "_80655_DLPFC.csv"))
    readr::write_csv(module_genes[[mod]], fn)
  }
}

message("WGCNA 80655 DLPFC concluído. Resultados em: ", OUT_DIR)
