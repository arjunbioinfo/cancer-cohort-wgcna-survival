library(testthat)

repo_root <- normalizePath(file.path(getwd(), "..", ".."), mustWork = FALSE)
if (!file.exists(file.path(repo_root, "R", "01_download_xena.R"))) {
  repo_root <- getwd()
}

source(file.path(repo_root, "R", "01_download_xena.R"))
source(file.path(repo_root, "R", "02_limma_de_gsea.R"))
source(file.path(repo_root, "R", "03_wgcna_modules.R"))
source(file.path(repo_root, "R", "04_survival_cox_km.R"))
source(file.path(repo_root, "R", "05_plot_summary.R"))

make_synthetic_cancer_cohort <- function(seed = 42) {
  set.seed(seed)
  n_samples <- 60
  samples <- sprintf("TCGA-00-%04d-01", seq_len(n_samples))
  stage_group <- rep(c("Early", "Advanced"), each = n_samples / 2)

  # Latent biological factors driving two co-expressed modules + survival risk
  f_prolif <- rnorm(n_samples) + ifelse(stage_group == "Advanced", 1.2, -1.2)
  f_immune <- rnorm(n_samples) + ifelse(stage_group == "Early", 0.9, -0.9)

  # Survival risk increases with f_prolif and decreases with f_immune
  hazard <- exp(0.9 * f_prolif - 0.6 * f_immune)
  os_time <- pmax(40, round(rexp(n_samples, rate = 0.0012 * hazard)))
  os_event <- ifelse(os_time < 1500, 1, 0)

  meta <- data.frame(
    sample = samples,
    OS = os_event,
    OS.time = os_time,
    stage_group = stage_group,
    stringsAsFactors = FALSE
  )

  prolif_genes <- unique(c(CANCER_HALLMARK_SETS$Cell_Cycle_Mitotic, paste0("PROL_", seq_len(40))))
  immune_genes <- unique(c(CANCER_HALLMARK_SETS$Immune_Interferon_Infiltration, paste0("IMMU_", seq_len(40))))
  emt_genes <- CANCER_HALLMARK_SETS$EMT_Extracellular_Matrix
  bg_genes <- paste0("BG_", seq_len(60))
  all_genes <- unique(c(prolif_genes, immune_genes, emt_genes, bg_genes))

  expr <- matrix(rnorm(length(all_genes) * n_samples, mean = 6, sd = 0.6),
                 nrow = length(all_genes), ncol = n_samples,
                 dimnames = list(all_genes, samples))

  for (g in prolif_genes) {
    expr[g, ] <- 7 + 1.4 * f_prolif + rnorm(n_samples, sd = 0.35)
  }
  for (g in immune_genes) {
    expr[g, ] <- 7 + 1.3 * f_immune + rnorm(n_samples, sd = 0.35)
  }
  for (g in emt_genes) {
    expr[g, ] <- 6.5 + 0.8 * f_prolif + rnorm(n_samples, sd = 0.5)
  }

  list(expr = expr, meta = meta)
}

test_that("parse_stage_group maps FIGO/AJCC stage strings cleanly", {
  parsed <- parse_stage_group(c("Stage IB1", "Stage IIB", "Stage IIIB", "Stage IVA", "Not Reported"))
  expect_equal(parsed, c("Early", "Early", "Advanced", "Advanced", NA_character_))
})

test_that("end-to-end limma, fgsea, WGCNA, Cox/KM survival, and plotting succeed on synthetic cohort", {
  cohort <- make_synthetic_cancer_cohort(seed = 42)

  de_res <- run_limma_and_gsea(cohort$expr, cohort$meta, min_mean_expr = 1.0, padj_thresh = 0.05, lfc_thresh = 0.5, seed = 42)
  expect_true(nrow(de_res$de) > 100)
  expect_true(sum(de_res$de$significant) >= 20)
  expect_true(nrow(de_res$gsea) >= 2)

  wgcna_res <- run_wgcna_analysis(
    cohort$expr, cohort$meta,
    top_mad_genes = nrow(cohort$expr),
    min_module_size = 15,
    merge_cut_height = 0.25,
    r2_target = 0.80,
    fallback_power = 6,
    seed = 42
  )
  expect_true(nrow(wgcna_res$membership) == nrow(cohort$expr))
  expect_true(nrow(wgcna_res$module_traits) >= 2)

  surv_res <- run_survival_analysis(cohort$expr, cohort$meta, wgcna_res, top_hub_per_module = 5)
  expect_true(nrow(surv_res$module_cox) >= 1)
  expect_true(nrow(surv_res$hub_cox) >= 3)

  tmp_dir <- tempdir()
  f1 <- file.path(tmp_dir, "volc_gsea.png")
  f2 <- file.path(tmp_dir, "wgcna_mt.png")
  f3 <- file.path(tmp_dir, "surv_km.png")
  f4 <- file.path(tmp_dir, "hub_genes.png")

  plot_limma_volcano_gsea(de_res$de, de_res$gsea, 0.05, 0.5, c("MCM5", "CD8A"), f1)
  plot_wgcna_module_traits(wgcna_res, f2)
  plot_survival_forest_km(cohort$meta, wgcna_res, surv_res, f3)
  plot_hub_genes_summary(surv_res$hub_cox, f4)

  for (f in c(f1, f2, f3, f4)) {
    expect_true(file.exists(f) && file.info(f)$size > 1000)
  }
})
