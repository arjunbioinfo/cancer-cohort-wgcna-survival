# 04_survival_cox_km.R — Cox proportional hazards and Kaplan-Meier overall survival modeling

suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
})

fit_single_cox_km <- function(feature_vec, os_time, os_event, stage_vec = NULL) {
  df <- data.frame(
    time = as.numeric(os_time) / 365.25,
    status = as.numeric(os_event),
    x = as.numeric(scale(feature_vec)),
    stringsAsFactors = FALSE
  )
  if (!is.null(stage_vec)) {
    df$stage <- factor(stage_vec, levels = c("Early", "Advanced"))
  }
  df <- df[is.finite(df$time) & is.finite(df$status) & is.finite(df$x), , drop = FALSE]
  if (nrow(df) < 20 || sum(df$status == 1) < 5) {
    return(NULL)
  }

  cox_uni <- tryCatch(survival::coxph(survival::Surv(time, status) ~ x, data = df), error = function(e) NULL)
  if (is.null(cox_uni)) return(NULL)
  s_uni <- summary(cox_uni)
  hr <- as.numeric(s_uni$conf.int[1, "exp(coef)"])
  ci_low <- as.numeric(s_uni$conf.int[1, "lower .95"])
  ci_high <- as.numeric(s_uni$conf.int[1, "upper .95"])
  cox_p <- as.numeric(s_uni$coefficients[1, "Pr(>|z|)"])

  # Stage-adjusted multivariate Cox when stage is available
  adj_hr <- NA_real_
  adj_p <- NA_real_
  if (!is.null(stage_vec) && sum(!is.na(df$stage)) >= 20 && length(unique(stats::na.omit(df$stage))) == 2) {
    cox_mv <- tryCatch(
      survival::coxph(survival::Surv(time, status) ~ x + stage, data = df[!is.na(df$stage), ]),
      error = function(e) NULL
    )
    if (!is.null(cox_mv)) {
      s_mv <- summary(cox_mv)
      adj_hr <- as.numeric(s_mv$conf.int["x", "exp(coef)"])
      adj_p <- as.numeric(s_mv$coefficients["x", "Pr(>|z|)"])
    }
  }

  # Median-split Kaplan-Meier log-rank test
  med_x <- stats::median(df$x, na.rm = TRUE)
  df$group <- ifelse(df$x >= med_x, "High", "Low")
  lr <- survival::survdiff(survival::Surv(time, status) ~ group, data = df)
  km_p <- stats::pchisq(lr$chisq, df = 1, lower.tail = FALSE)

  data.frame(
    n_samples = nrow(df),
    n_events = sum(df$status == 1),
    hr = hr,
    ci_low = ci_low,
    ci_high = ci_high,
    cox_pvalue = cox_p,
    adj_hr = adj_hr,
    adj_pvalue = adj_p,
    km_logrank_p = km_p,
    stringsAsFactors = FALSE
  )
}

run_survival_analysis <- function(expr_mat, meta, wgcna_res, top_hub_per_module = 10) {
  MEs <- wgcna_res$MEs[meta$sample, , drop = FALSE]
  me_names <- setdiff(colnames(MEs), "MEgrey")
  if (length(me_names) == 0) me_names <- colnames(MEs)

  # 1. Module-level Cox & KM
  mod_rows <- list()
  for (m in me_names) {
    res <- fit_single_cox_km(MEs[[m]], meta$OS.time, meta$OS, meta$stage_group)
    if (!is.null(res)) {
      res$module <- m
      res$n_genes <- sum(wgcna_res$membership$module == m)
      mod_rows[[length(mod_rows) + 1]] <- res
    }
  }
  module_cox_df <- bind_rows(mod_rows) %>%
    mutate(fdr_bh = stats::p.adjust(cox_pvalue, method = "BH")) %>%
    select(module, n_genes, n_samples, n_events, hr, ci_low, ci_high, cox_pvalue, fdr_bh, adj_hr, adj_pvalue, km_logrank_p) %>%
    arrange(cox_pvalue)

  # 2. Hub-gene-level Cox & KM for top 3 modules by Cox p-value
  top_mods <- head(module_cox_df$module, 3)
  hub_candidates <- wgcna_res$membership %>%
    filter(module %in% top_mods) %>%
    group_by(module) %>%
    slice_max(order_by = abs(kME), n = top_hub_per_module, with_ties = FALSE) %>%
    ungroup()

  hub_rows <- list()
  for (i in seq_len(nrow(hub_candidates))) {
    g <- hub_candidates$gene[i]
    if (!(g %in% rownames(expr_mat))) next
    res <- fit_single_cox_km(expr_mat[g, meta$sample], meta$OS.time, meta$OS, meta$stage_group)
    if (!is.null(res)) {
      res$gene <- g
      res$module <- hub_candidates$module[i]
      res$kME <- hub_candidates$kME[i]
      hub_rows[[length(hub_rows) + 1]] <- res
    }
  }
  hub_cox_df <- bind_rows(hub_rows) %>%
    mutate(fdr_bh = stats::p.adjust(cox_pvalue, method = "BH")) %>%
    select(gene, module, kME, hr, ci_low, ci_high, cox_pvalue, fdr_bh, adj_hr, adj_pvalue, km_logrank_p) %>%
    arrange(cox_pvalue)

  message(sprintf(
    "[INFO] Survival Cox/KM evaluated %d WGCNA modules and %d hub genes (top module: %s, HR=%.2f, p=%.3g).",
    nrow(module_cox_df), nrow(hub_cox_df),
    module_cox_df$module[1], module_cox_df$hr[1], module_cox_df$cox_pvalue[1]
  ))

  list(module_cox = module_cox_df, hub_cox = hub_cox_df)
}
