# 05_plot_summary.R — Publication-quality ggplot2/patchwork summary figures (300 DPI)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(survival)
})

plot_limma_volcano_gsea <- function(de_df, gsea_df, padj_thresh, lfc_thresh, highlight_genes, out_path) {
  df <- de_df %>%
    mutate(
      neg_log10_p = -log10(pmax(adj.P.Val, 1e-300)),
      status = case_when(
        adj.P.Val < padj_thresh & logFC >= lfc_thresh ~ "Up",
        adj.P.Val < padj_thresh & logFC <= -lfc_thresh ~ "Down",
        TRUE ~ "NS"
      )
    )
  contrast_lbl <- df$contrast[1]

  hl <- df %>% filter(gene %in% highlight_genes)
  if (nrow(hl) == 0) {
    hl <- df %>% arrange(adj.P.Val) %>% head(6)
  }

  p1 <- ggplot(df, aes(x = logFC, y = neg_log10_p, color = status)) +
    geom_point(alpha = 0.65, size = 1.3) +
    scale_color_manual(values = c(Up = "#D95F02", Down = "#1B9E77", NS = "#BDBDBD")) +
    geom_hline(yintercept = -log10(padj_thresh), linetype = "dashed", color = "#555555") +
    geom_vline(xintercept = c(-lfc_thresh, lfc_thresh), linetype = "dashed", color = "#555555") +
    geom_text(data = hl, aes(label = gene), color = "#111111", size = 2.8, vjust = -0.7, fontface = "bold", check_overlap = TRUE) +
    labs(
      title = paste0("limma DE (", gsub("_", " ", contrast_lbl), ")"),
      x = "log2 Fold Change",
      y = "-log10(FDR)",
      color = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "top", plot.title = element_text(face = "bold"))

  g_plot <- gsea_df %>%
    mutate(
      pathway_clean = gsub("_", " ", pathway),
      sig_label = ifelse(padj < 0.05, "FDR < 0.05", "NS")
    )

  p2 <- ggplot(g_plot, aes(x = stats::reorder(pathway_clean, NES), y = NES, fill = sig_label)) +
    geom_col(width = 0.68) +
    coord_flip() +
    geom_hline(yintercept = 0, color = "#333333") +
    scale_fill_manual(values = c("FDR < 0.05" = "#386CB0", "NS" = "#BDBDBD")) +
    labs(
      title = "fgsea Cancer Pathway Enrichment",
      x = NULL,
      y = "Normalized Enrichment Score (NES)",
      fill = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "top", plot.title = element_text(face = "bold"))

  combined <- p1 + p2 + plot_layout(widths = c(1.1, 1.0))
  ggsave(out_path, combined, width = 12, height = 5.0, dpi = 300, bg = "white")
}

plot_wgcna_module_traits <- function(wgcna_res, out_path) {
  sft <- wgcna_res$sft_fit
  p1 <- ggplot(sft, aes(x = Power, y = signed_R2)) +
    geom_line(color = "#386CB0", linewidth = 0.9) +
    geom_point(size = 2.5, color = "#386CB0") +
    geom_hline(yintercept = 0.80, linetype = "dashed", color = "#D95F02") +
    geom_vline(xintercept = wgcna_res$soft_power, linetype = "dotted", color = "#1B9E77") +
    labs(
      title = sprintf("WGCNA Soft-Thresholding (selected \u03b2 = %d)", wgcna_res$soft_power),
      x = "Soft-thresholding power (\u03b2)",
      y = "Signed Scale-Free Topology R\u00b2"
    ) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))

  mt <- wgcna_res$module_traits %>%
    mutate(
      label = sprintf("%.2f\n(p=%.2g)", cor, pvalue),
      trait_clean = gsub("_", " ", trait)
    )

  p2 <- ggplot(mt, aes(x = trait_clean, y = module, fill = cor)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = label), size = 2.7, color = "#111111") +
    scale_fill_gradient2(low = "#1B9E77", mid = "white", high = "#D95F02", midpoint = 0, limits = c(-1, 1)) +
    labs(
      title = "Module Eigengene\u2013Clinical Trait Correlations",
      x = NULL,
      y = NULL,
      fill = "Pearson r"
    ) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))

  combined <- p1 + p2 + plot_layout(widths = c(0.85, 1.15))
  ggsave(out_path, combined, width = 12, height = 5.4, dpi = 300, bg = "white")
}

plot_survival_forest_km <- function(meta, wgcna_res, surv_res, out_path) {
  mod_cox <- surv_res$module_cox %>%
    head(12) %>%
    mutate(
      ci_low_c = pmax(ci_low, 0.1),
      ci_high_c = pmin(ci_high, 10),
      sig = ifelse(cox_pvalue < 0.05, "p < 0.05", "NS")
    )

  p1 <- ggplot(mod_cox, aes(x = hr, y = stats::reorder(module, hr), color = sig)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "#555555") +
    geom_errorbarh(aes(xmin = ci_low_c, xmax = ci_high_c), height = 0.25, linewidth = 0.8) +
    geom_point(size = 2.8) +
    scale_x_log10() +
    scale_color_manual(values = c("p < 0.05" = "#D95F02", "NS" = "#636363")) +
    labs(
      title = "Univariate Cox PH: WGCNA Module Eigengenes",
      x = "Hazard Ratio (log scale, 95% CI)",
      y = NULL,
      color = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "top", plot.title = element_text(face = "bold"))

  # Kaplan-Meier step curve for top prognostic module
  top_mod <- surv_res$module_cox$module[1]
  km_p <- surv_res$module_cox$km_logrank_p[1]
  me_vals <- wgcna_res$MEs[meta$sample, top_mod]
  km_df <- data.frame(
    time = meta$OS.time / 365.25,
    status = meta$OS,
    group = ifelse(me_vals >= stats::median(me_vals, na.rm = TRUE), "High ME", "Low ME")
  )
  fit_km <- survival::survfit(survival::Surv(time, status) ~ group, data = km_df)

  # Convert survfit to tidy data frame for ggplot2
  strata_names <- names(fit_km$strata)
  km_plot_df <- data.frame(
    time = fit_km$time,
    surv = fit_km$surv,
    lower = fit_km$lower,
    upper = fit_km$upper,
    group = rep(gsub("^group=", "", strata_names), fit_km$strata),
    stringsAsFactors = FALSE
  )

  p2 <- ggplot(km_plot_df, aes(x = time, y = surv, color = group, fill = group)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.14, color = NA) +
    geom_step(linewidth = 1.0) +
    scale_color_manual(values = c("High ME" = "#D95F02", "Low ME" = "#1B9E77")) +
    scale_fill_manual(values = c("High ME" = "#D95F02", "Low ME" = "#1B9E77")) +
    coord_cartesian(ylim = c(0, 1.02)) +
    annotate("text", x = max(km_plot_df$time) * 0.65, y = 0.18,
             label = sprintf("Log-rank p = %.3g", km_p), fontface = "bold", size = 3.6) +
    labs(
      title = sprintf("Kaplan\u2013Meier Overall Survival (%s)", top_mod),
      x = "Overall Survival (years)",
      y = "Survival Probability",
      color = NULL, fill = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "top", plot.title = element_text(face = "bold"))

  combined <- p1 + p2 + plot_layout(widths = c(1.0, 1.05))
  ggsave(out_path, combined, width = 12, height = 5.0, dpi = 300, bg = "white")
}

plot_hub_genes_summary <- function(hub_cox_df, out_path) {
  df <- hub_cox_df %>%
    head(15) %>%
    mutate(
      direction = ifelse(hr >= 1, "Risk (HR > 1)", "Protective (HR < 1)"),
      neg_log10_p = -log10(pmax(cox_pvalue, 1e-100))
    )

  p <- ggplot(df, aes(x = abs(kME), y = hr, size = neg_log10_p, color = direction)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "#555555") +
    geom_point(alpha = 0.85) +
    geom_text(aes(label = gene), size = 3.1, vjust = -0.9, show.legend = FALSE, check_overlap = TRUE) +
    scale_color_manual(values = c("Risk (HR > 1)" = "#D95F02", "Protective (HR < 1)" = "#1B9E77")) +
    labs(
      title = "Top Prognostic WGCNA Hub Genes (Intramodular Connectivity vs. Hazard Ratio)",
      x = "Intramodular Connectivity |kME|",
      y = "Univariate Cox Hazard Ratio (per SD)",
      color = NULL,
      size = "-log10(p)"
    ) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))

  ggsave(out_path, p, width = 8.4, height = 5.4, dpi = 300, bg = "white")
}
