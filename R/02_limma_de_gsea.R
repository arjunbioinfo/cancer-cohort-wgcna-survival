# 02_limma_de_gsea.R — limma empirical Bayes differential expression and fgsea cancer pathway enrichment

suppressPackageStartupMessages({
  library(dplyr)
  library(fgsea)
  library(limma)
})

CANCER_HALLMARK_SETS <- list(
  Cell_Cycle_Mitotic = c(
    "MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "MCM7", "CDC45", "KIF11", "KIF2C",
    "CDK1", "CCNB1", "CCNB2", "CCNA2", "MKI67", "TOP2A", "BUB1", "BUB1B", "AURKA", "AURKB", "PLK1", "TPX2", "PCNA"
  ),
  EMT_Extracellular_Matrix = c(
    "FN1", "CDH2", "VIM", "SNAI1", "SNAI2", "ZEB1", "ZEB2", "TWIST1",
    "COL1A1", "COL1A2", "COL3A1", "MMP2", "MMP9", "MMP14", "LOX", "POSTN", "SPARC", "TAGLN", "ACTA2", "ITGA5"
  ),
  DNA_Replication_Repair = c(
    "BRCA1", "BRCA2", "RAD51", "CHEK1", "CHEK2", "ATR", "ATM", "RFC4", "FEN1",
    "EXO1", "POLE", "POLD1", "RRM1", "RRM2", "TYMS", "UNG", "MSH2", "MSH6", "FANCD2", "BARD1"
  ),
  MYC_Proliferation_Targets = c(
    "MYC", "ODC1", "NPM1", "NCL", "NOP56", "FBL", "DKC1", "PA2G4", "LDHA",
    "ENO1", "PGK1", "PKM", "SLC2A1", "SRM", "CAD", "EIF4E", "EIF4A1", "CDK4", "CCND2", "HSP90AA1"
  ),
  Immune_Interferon_Infiltration = c(
    "CD8A", "CD8B", "CD3D", "CD3E", "GZMA", "GZMB", "PRF1", "NKG7", "IFNG",
    "CXCL9", "CXCL10", "STAT1", "ISG15", "IFIT1", "IFIT3", "MX1", "HLA-DRA", "HLA-DRB1", "CD74", "PDCD1"
  ),
  Glycolysis_Hypoxia_Metabolism = c(
    "VEGFA", "SLC2A1", "HK2", "PFKP", "ALDOA", "GAPDH", "PGK1", "ENO1",
    "LDHA", "PDK1", "BNIP3", "CA9", "ADM", "NDRG1", "P4HA1", "EGLN3", "HIF1A", "TPI1", "GPI"
  ),
  p53_Apoptosis_Signaling = c(
    "TP53", "CDKN1A", "MDM2", "BAX", "BBC3", "PMAIP1", "GADD45A", "FAS",
    "TNFRSF10B", "CASP3", "CASP8", "CASP9", "BID", "BCL2", "MCL1", "BIRC5", "APAF1", "BTG2"
  ),
  Keratinization_Squamous_Diff = c(
    "KRT5", "KRT6A", "KRT6B", "KRT14", "KRT16", "KRT17", "TP63", "SOX2",
    "IVL", "SPRR1A", "SPRR1B", "SPRR2A", "S100A7", "S100A8", "S100A9", "DSG3", "DSC3", "PKP1"
  )
)

run_limma_and_gsea <- function(expr_mat, meta, min_mean_expr = 1.0, padj_thresh = 0.05, lfc_thresh = 0.5, seed = 42) {
  # Prefer Advanced vs Early stage; fall back to OS Deceased vs Alive if stage has too few samples
  has_stage <- !is.na(meta$stage_group)
  if (sum(meta$stage_group == "Advanced", na.rm = TRUE) >= 10 &&
      sum(meta$stage_group == "Early", na.rm = TRUE) >= 10) {
    sub_meta <- meta[has_stage, , drop = FALSE]
    group <- factor(sub_meta$stage_group, levels = c("Early", "Advanced"))
    contrast_label <- "Advanced_vs_Early"
  } else {
    sub_meta <- meta
    group <- factor(ifelse(sub_meta$OS == 1, "Deceased", "Alive"), levels = c("Alive", "Deceased"))
    contrast_label <- "Deceased_vs_Alive"
  }

  sub_expr <- expr_mat[, sub_meta$sample, drop = FALSE]
  keep_genes <- rowMeans(sub_expr, na.rm = TRUE) >= min_mean_expr &
    apply(sub_expr, 1, stats::sd, na.rm = TRUE) > 1e-6
  sub_expr <- sub_expr[keep_genes, , drop = FALSE]

  design <- stats::model.matrix(~ group)
  fit <- limma::lmFit(sub_expr, design)
  fit <- limma::eBayes(fit)
  tt <- limma::topTable(fit, coef = 2, number = Inf, sort.by = "none")
  tt$gene <- rownames(tt)
  tt$contrast <- contrast_label
  tt$significant <- (tt$adj.P.Val < padj_thresh) & (abs(tt$logFC) >= lfc_thresh)
  de_df <- tt %>%
    select(gene, contrast, logFC, AveExpr, t, P.Value, adj.P.Val, significant) %>%
    arrange(adj.P.Val, desc(abs(logFC)))

  # Run fgsea using moderated t-statistics
  ranks <- stats::setNames(de_df$t, de_df$gene)
  ranks <- ranks[is.finite(ranks)]
  ranks <- sort(ranks, decreasing = TRUE)

  set.seed(seed)
  fg_res <- suppressWarnings(
    fgsea::fgseaMultilevel(
      pathways = CANCER_HALLMARK_SETS,
      stats = ranks,
      minSize = 5,
      maxSize = 500
    )
  )
  gsea_df <- as.data.frame(fg_res) %>%
    mutate(leadingEdge = vapply(leadingEdge, function(x) paste(head(x, 8), collapse = ","), character(1))) %>%
    select(pathway, pval, padj, ES, NES, size, leadingEdge) %>%
    arrange(padj, desc(abs(NES)))

  message(sprintf(
    "[INFO] limma (%s): %d genes tested, %d significant DEGs; fgsea tested %d pathways.",
    contrast_label, nrow(de_df), sum(de_df$significant), nrow(gsea_df)
  ))

  list(de = de_df, gsea = gsea_df, contrast = contrast_label)
}
