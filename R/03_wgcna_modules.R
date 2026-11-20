# 03_wgcna_modules.R — Signed WGCNA co-expression network, module-trait correlations, and hub genes

suppressPackageStartupMessages({
  library(dplyr)
  library(WGCNA)
})

run_wgcna_analysis <- function(expr_mat, meta, top_mad_genes = 1200, min_module_size = 25,
                               merge_cut_height = 0.25, r2_target = 0.80, fallback_power = 6, seed = 42) {
  set.seed(seed)
  WGCNA::disableWGCNAThreads()

  # Select top variable genes by MAD for scalable, focused co-expression network
  mads <- apply(expr_mat, 1, stats::mad, na.rm = TRUE)
  mads <- mads[is.finite(mads) & mads > 0]
  n_keep <- min(top_mad_genes, length(mads))
  selected_genes <- names(sort(mads, decreasing = TRUE))[seq_len(n_keep)]
  datExpr <- t(expr_mat[selected_genes, meta$sample, drop = FALSE])

  # Soft-thresholding power selection
  powers <- c(seq(2, 12, by = 2), 14, 16)
  sft <- suppressWarnings(
    WGCNA::pickSoftThreshold(datExpr, powerVector = powers, networkType = "signed", verbose = 0)
  )
  fit_df <- sft$fitIndices
  fit_df$signed_R2 <- -sign(fit_df$slope) * fit_df$SFT.R.sq

  qualifying <- fit_df$Power[fit_df$signed_R2 >= r2_target]
  soft_power <- if (length(qualifying) > 0) qualifying[1] else fallback_power
  message(sprintf("[INFO] WGCNA selected soft-thresholding power = %d (%d genes x %d samples).",
                  soft_power, ncol(datExpr), nrow(datExpr)))

  # Guard against WGCNA::cor namespace clash
  old_cor <- stats::cor
  cor <- WGCNA::cor
  on.exit({ cor <- old_cor }, add = TRUE)

  net <- WGCNA::blockwiseModules(
    datExpr,
    power = soft_power,
    networkType = "signed",
    TOMType = "signed",
    minModuleSize = min_module_size,
    reassignThreshold = 0,
    mergeCutHeight = merge_cut_height,
    numericLabels = FALSE,
    pamRespectsDendro = FALSE,
    saveTOMs = FALSE,
    maxBlockSize = max(ncol(datExpr) + 100, 5000),
    verbose = 0
  )

  module_colors <- net$colors
  names(module_colors) <- colnames(datExpr)
  MEs <- WGCNA::orderMEs(net$MEs)
  rownames(MEs) <- rownames(datExpr)

  # Module-trait correlation against OS event, OS time (years), and Advanced stage
  traits <- data.frame(
    OS_event = as.numeric(meta$OS),
    OS_years = as.numeric(meta$OS.time) / 365.25,
    Advanced_stage = ifelse(is.na(meta$stage_group), NA_real_, as.numeric(meta$stage_group == "Advanced")),
    row.names = meta$sample
  )
  # Drop all-NA or zero-variance traits
  valid_traits <- vapply(traits, function(col) sum(!is.na(col)) >= 10 && stats::sd(col, na.rm = TRUE) > 0, logical(1))
  traits <- traits[, valid_traits, drop = FALSE]

  cor_mat <- stats::cor(MEs, traits, use = "pairwise.complete.obs")
  pval_mat <- WGCNA::corPvalueStudent(cor_mat, nSamples = nrow(datExpr))

  mod_trait_rows <- list()
  for (m in rownames(cor_mat)) {
    for (tr in colnames(cor_mat)) {
      mod_trait_rows[[length(mod_trait_rows) + 1]] <- data.frame(
        module = m,
        trait = tr,
        cor = as.numeric(cor_mat[m, tr]),
        pvalue = as.numeric(pval_mat[m, tr]),
        n_genes = sum(paste0("ME", module_colors) == m),
        stringsAsFactors = FALSE
      )
    }
  }
  module_traits_df <- bind_rows(mod_trait_rows) %>%
    mutate(fdr_bh = stats::p.adjust(pvalue, method = "BH"))

  # Intra-modular connectivity (kME) per gene
  kme_mat <- stats::cor(datExpr, MEs, use = "pairwise.complete.obs")
  own_me <- paste0("ME", module_colors)
  kme_own <- vapply(seq_along(module_colors), function(i) {
    me_name <- own_me[i]
    if (me_name %in% colnames(kme_mat)) as.numeric(kme_mat[i, me_name]) else NA_real_
  }, numeric(1))

  membership_df <- data.frame(
    gene = names(module_colors),
    module = paste0("ME", module_colors),
    module_color = as.character(module_colors),
    kME = kme_own,
    stringsAsFactors = FALSE
  ) %>%
    arrange(module, desc(abs(kME)))

  list(
    MEs = MEs,
    membership = membership_df,
    module_traits = module_traits_df,
    sft_fit = fit_df %>% select(Power, signed_R2, mean.k.),
    soft_power = soft_power
  )
}
