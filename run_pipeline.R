#!/usr/bin/env Rscript
# run_pipeline.R — End-to-end TCGA cancer cohort DE, WGCNA, GSEA, and survival pipeline

suppressPackageStartupMessages({
  library(yaml)
  library(data.table)
})

parse_config_arg <- function(args) {
  idx <- which(args == "--config")
  if (length(idx) == 1 && length(args) >= idx + 1) {
    return(args[idx + 1])
  }
  "config/config.yaml"
}

args <- commandArgs(trailingOnly = TRUE)
config_path <- parse_config_arg(args)
message(sprintf("[INFO] Loading config: %s", config_path))
cfg <- yaml::read_yaml(config_path)

dir.create(cfg$output_dirs$data, recursive = TRUE, showWarnings = FALSE)
dir.create(cfg$output_dirs$results, recursive = TRUE, showWarnings = FALSE)
dir.create(cfg$output_dirs$figures, recursive = TRUE, showWarnings = FALSE)

source("R/01_download_xena.R")
source("R/02_limma_de_gsea.R")
source("R/03_wgcna_modules.R")
source("R/04_survival_cox_km.R")
source("R/05_plot_summary.R")

# 1. Download & align cohort
paths <- download_xena_cohort(cfg)
cohort <- load_and_align_xena_cohort(
  paths$expr_path, paths$surv_path, paths$clin_path,
  min_os_days = cfg$survival$min_os_days
)

# 2. limma DE + fgsea
de_res <- run_limma_and_gsea(
  cohort$expr, cohort$meta,
  min_mean_expr = cfg$de$min_mean_expr,
  padj_thresh = cfg$de$padj_threshold,
  lfc_thresh = cfg$de$lfc_threshold,
  seed = cfg$random_seed
)

# 3. Signed WGCNA co-expression modules
wgcna_res <- run_wgcna_analysis(
  cohort$expr, cohort$meta,
  top_mad_genes = cfg$wgcna$top_mad_genes,
  min_module_size = cfg$wgcna$min_module_size,
  merge_cut_height = cfg$wgcna$merge_cut_height,
  r2_target = cfg$wgcna$r2_target,
  fallback_power = cfg$wgcna$fallback_power,
  seed = cfg$random_seed
)

# 4. Survival Cox PH + Kaplan-Meier
surv_res <- run_survival_analysis(
  cohort$expr, cohort$meta, wgcna_res,
  top_hub_per_module = cfg$survival$top_hub_genes_per_module
)

# 5. Export tables
res_dir <- cfg$output_dirs$results
fwrite(de_res$de, file.path(res_dir, "limma_de_results.tsv"), sep = "\t")
fwrite(de_res$gsea, file.path(res_dir, "fgsea_pathways.tsv"), sep = "\t")
fwrite(wgcna_res$membership, file.path(res_dir, "wgcna_module_membership.tsv"), sep = "\t")
fwrite(wgcna_res$module_traits, file.path(res_dir, "wgcna_module_traits.tsv"), sep = "\t")
fwrite(surv_res$module_cox, file.path(res_dir, "survival_module_cox.tsv"), sep = "\t")
fwrite(surv_res$hub_cox, file.path(res_dir, "survival_hub_genes.tsv"), sep = "\t")

# 6. Generate summary figures
fig_dir <- cfg$output_dirs$figures
plot_limma_volcano_gsea(
  de_res$de, de_res$gsea,
  padj_thresh = cfg$de$padj_threshold,
  lfc_thresh = cfg$de$lfc_threshold,
  highlight_genes = cfg$plotting$highlight_genes,
  out_path = file.path(fig_dir, "limma_volcano_gsea.png")
)
plot_wgcna_module_traits(wgcna_res, file.path(fig_dir, "wgcna_module_trait_heatmap.png"))
plot_survival_forest_km(cohort$meta, wgcna_res, surv_res, file.path(fig_dir, "survival_forest_km.png"))
plot_hub_genes_summary(surv_res$hub_cox, file.path(fig_dir, "hub_genes_summary.png"))

message("[INFO] Pipeline complete. All tables and figures saved.")
