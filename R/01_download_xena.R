# 01_download_xena.R — Download and align TCGA cohort expression, survival, and clinical data from UCSC Xena

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(jsonlite)
})

download_if_missing <- function(url, dest) {
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(dest) && file.info(dest)$size > 0) {
    message(sprintf("[INFO] Using cached file: %s (%d bytes)", dest, file.info(dest)$size))
    return(file.info(dest)$size)
  }
  message(sprintf("[INFO] Downloading %s -> %s", url, dest))
  old_timeout <- getOption("timeout")
  on.exit(options(timeout = old_timeout), add = TRUE)
  options(timeout = 300)
  status <- utils::download.file(url, destfile = dest, mode = "wb", quiet = TRUE)
  if (status != 0 || !file.exists(dest) || file.info(dest)$size == 0) {
    stop(sprintf("Failed to download %s", url))
  }
  file.info(dest)$size
}

parse_stage_group <- function(stage_vec) {
  s <- toupper(trimws(as.character(stage_vec)))
  out <- rep(NA_character_, length(s))
  out[grepl("STAGE\\s*(I($|[^IV])|II($|[^I]))", s)] <- "Early"
  out[grepl("STAGE\\s*(III|IV)", s)] <- "Advanced"
  out
}

download_xena_cohort <- function(cfg) {
  raw_dir <- cfg$output_dirs$data
  cohort_tag <- gsub("[^A-Za-z0-9_-]", "_", cfg$cohort$name)

  expr_path <- file.path(raw_dir, paste0(cohort_tag, "_HiSeqV2.gz"))
  surv_path <- file.path(raw_dir, paste0(cohort_tag, "_survival.txt"))
  clin_path <- file.path(raw_dir, paste0(cohort_tag, "_clinicalMatrix.tsv"))

  s_expr <- download_if_missing(cfg$cohort$expr_url, expr_path)
  s_surv <- download_if_missing(cfg$cohort$survival_url, surv_path)
  s_clin <- download_if_missing(cfg$cohort$clinical_url, clin_path)

  manifest <- list(
    cohort = cfg$cohort$name,
    files = list(
      expression = list(url = cfg$cohort$expr_url, path = expr_path, bytes = s_expr),
      survival = list(url = cfg$cohort$survival_url, path = surv_path, bytes = s_surv),
      clinical = list(url = cfg$cohort$clinical_url, path = clin_path, bytes = s_clin)
    )
  )
  write_json(manifest, file.path(raw_dir, "download_manifest.json"), pretty = TRUE, auto_unbox = TRUE)

  list(expr_path = expr_path, surv_path = surv_path, clin_path = clin_path)
}

load_and_align_xena_cohort <- function(expr_path, surv_path, clin_path, min_os_days = 30) {
  # Expression: rows = HGNC genes, cols = TCGA sample barcodes
  expr_dt <- fread(expr_path, data.table = FALSE)
  gene_col <- colnames(expr_dt)[1]
  genes <- make.unique(as.character(expr_dt[[gene_col]]))
  expr_mat <- as.matrix(expr_dt[, -1, drop = FALSE])
  rownames(expr_mat) <- genes
  storage.mode(expr_mat) <- "numeric"

  # Keep primary tumor samples (-01)
  is_primary <- grepl("-01$", colnames(expr_mat))
  if (sum(is_primary) >= 20) {
    expr_mat <- expr_mat[, is_primary, drop = FALSE]
  }

  # Survival table
  surv_df <- fread(surv_path, data.table = FALSE)
  sample_col <- intersect(c("sample", "Sample", "_PATIENT", "bcr_patient_barcode"), colnames(surv_df))[1]
  if (is.na(sample_col)) stop("No sample identifier column found in survival table.")
  surv_clean <- data.frame(
    sample = as.character(surv_df[[sample_col]]),
    OS = suppressWarnings(as.numeric(surv_df[["OS"]])),
    OS.time = suppressWarnings(as.numeric(surv_df[["OS.time"]])),
    stringsAsFactors = FALSE
  ) %>%
    filter(!is.na(OS), !is.na(OS.time), OS.time >= min_os_days) %>%
    distinct(sample, .keep_all = TRUE)

  # Clinical metadata (optional stage extraction)
  clin_df <- fread(clin_path, data.table = FALSE)
  clin_id_col <- colnames(clin_df)[1]
  stage_candidates <- c("clinical_stage", "pathologic_stage", "ajcc_pathologic_tumor_stage", "tumor_stage")
  stage_col <- intersect(stage_candidates, colnames(clin_df))[1]
  clin_clean <- data.frame(
    sample = as.character(clin_df[[clin_id_col]]),
    stage_group = if (!is.na(stage_col)) parse_stage_group(clin_df[[stage_col]]) else NA_character_,
    stringsAsFactors = FALSE
  ) %>%
    distinct(sample, .keep_all = TRUE)

  meta <- surv_clean %>%
    left_join(clin_clean, by = "sample") %>%
    filter(sample %in% colnames(expr_mat)) %>%
    arrange(sample)

  if (nrow(meta) < 20) {
    stop(sprintf("Too few aligned tumor samples with survival data (%d).", nrow(meta)))
  }

  expr_aligned <- expr_mat[, meta$sample, drop = FALSE]
  # Replace any residual NAs with row medians
  if (anyNA(expr_aligned)) {
    row_meds <- apply(expr_aligned, 1, median, na.rm = TRUE)
    na_idx <- which(is.na(expr_aligned), arr.ind = TRUE)
    expr_aligned[na_idx] <- row_meds[na_idx[, 1]]
  }

  message(sprintf(
    "[INFO] Aligned %d genes x %d primary tumors (OS events=%d; Early=%d, Advanced=%d).",
    nrow(expr_aligned), ncol(expr_aligned), sum(meta$OS == 1),
    sum(meta$stage_group == "Early", na.rm = TRUE),
    sum(meta$stage_group == "Advanced", na.rm = TRUE)
  ))

  list(expr = expr_aligned, meta = meta)
}
