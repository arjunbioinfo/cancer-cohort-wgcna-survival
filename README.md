# cancer-cohort-wgcna-survival

[![CI](https://github.com/arjunbioinfo/cancer-cohort-wgcna-survival/actions/workflows/ci.yml/badge.svg)](https://github.com/arjunbioinfo/cancer-cohort-wgcna-survival/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

A reproducible, config-driven **R / Bioconductor** workflow for bulk cancer cohort transcriptomics combining **limma** empirical Bayes differential expression, **WGCNA** signed co-expression modules, **fgsea** cancer pathway enrichment, and **Cox proportional hazards + Kaplan–Meier** overall survival modeling.

Every push runs offline `testthat` unit tests and executes the full workflow on the public **TCGA-CESC** (Cervical Squamous Cell Carcinoma) cohort from **UCSC Xena** in GitHub Actions CI.

## Pipeline Overview

1. **Cohort Download & Alignment (`R/01_download_xena.R`)**: Fetches `TCGA-CESC` log2-normalized gene expression (`HiSeqV2.gz`), overall survival (`CESC_survival.txt`), and clinical stage metadata (`CESC_clinicalMatrix`) from UCSC Xena ([UCSC Xena, TCGA-CESC cohort](https://xenabrowser.net/datapages/?cohort=TCGA%20Cervical%20Cancer%20(CESC))), filters to primary tumors (`-01`) with `OS.time >= 30` days, and writes `data/raw/download_manifest.json`.
2. **Differential Expression & GSEA (`R/02_limma_de_gsea.R`)**: Fits `limma` empirical Bayes linear models comparing Advanced (Stage III–IV) vs. Early (Stage I–II) tumors and runs `fgseaMultilevel` across curated cancer hallmark pathways (cell cycle/mitotic, EMT/ECM, DNA repair, MYC targets, immune/interferon infiltration, glycolysis/hypoxia, p53/apoptosis, squamous differentiation).
3. **Signed WGCNA Co-Expression Network (`R/03_wgcna_modules.R`)**: Selects top MAD genes, chooses a scale-free soft-thresholding power (`pickSoftThreshold`), detects signed co-expression modules (`blockwiseModules`), correlates module eigengenes with stage and survival, and computes intramodular connectivity (`kME`) to identify hub genes.
4. **Survival Modeling (`R/04_survival_cox_km.R`)**: Runs univariate and stage-adjusted multivariate Cox proportional hazards regression (`survival::coxph`) plus median-split Kaplan–Meier log-rank tests (`survival::survdiff`) on every WGCNA module eigengene and top intramodular hub genes.
5. **Visualization (`R/05_plot_summary.R`)**: Exports four 300 DPI `ggplot2` / `patchwork` summary figures in `figures/`.

## Quick Start

```bash
# 1. Install dependencies (R 4.4+ with pak)
Rscript -e "install.packages('pak'); pak::local_install_deps(dependencies = TRUE)"

# 2. Run offline unit tests
make test

# 3. Fast end-to-end run on real TCGA-CESC data (CI profile, top 1,200 MAD genes)
make ci

# 4. Full cohort run (top 3,000 MAD genes)
make run

# Or run via Docker
docker build -t cancer-cohort-wgcna-survival .
docker run --rm -v "$(pwd)/results:/app/results" -v "$(pwd)/figures:/app/figures" cancer-cohort-wgcna-survival
```

## Outputs

* `results/limma_de_results.tsv` — Full `limma` moderated t-test differential expression table (`logFC`, `AveExpr`, `t`, `P.Value`, `adj.P.Val`).
* `results/fgsea_pathways.tsv` — `fgsea` normalized enrichment scores (`NES`), p-values, BH-FDR, and leading-edge genes.
* `results/wgcna_module_membership.tsv` — Gene-to-module assignments and intramodular connectivity (`kME`).
* `results/wgcna_module_traits.tsv` — Module eigengene correlations and p-values against OS event, OS time, and clinical stage.
* `results/survival_module_cox.tsv` — Univariate & stage-adjusted Cox HR, 95% CI, and Kaplan–Meier log-rank p-values per WGCNA module.
* `results/survival_hub_genes.tsv` — Prognostic Cox/KM statistics for top intramodular hub genes.
* `figures/` — `limma_volcano_gsea.png`, `wgcna_module_trait_heatmap.png`, `survival_forest_km.png`, and `hub_genes_summary.png`.

## Adapting to Another TCGA Cohort

Swap the three UCSC Xena URLs in `config/config.yaml` (e.g., `TCGA.LUSC`, `TCGA.LUAD`, or `TCGA.BRCA`) and re-run `make run`.

## Repository Layout

* `run_pipeline.R` — Single-command workflow runner
* `R/` — Modular scripts (`01_download_xena.R` through `05_plot_summary.R`)
* `config/config.yaml` / `config/ci.yaml` — Full and fast CI TCGA-CESC configurations
* `DESCRIPTION` — R/Bioconductor dependency manifest for `pak` and GitHub Actions
* `tests/testthat/test_pipeline.R` — Offline synthetic cohort unit + smoke tests
* `Makefile`, `Dockerfile`, `.dockerignore`, `CITATION.cff`
* `.github/workflows/ci.yml` — GitHub Actions CI (R 4.4 + Posit Package Manager binaries + real TCGA-CESC run)

## License

MIT — see [LICENSE](LICENSE).

## Author

**Mallikarjuna Thippana, PhD** — Computational Biologist, Karolinska Institutet  
[Portfolio](https://arjunbioinfo.github.io/) · [ORCID](https://orcid.org/0000-0001-8485-8626) · [Google Scholar](https://scholar.google.com/citations?user=NpLUj68AAAAJ) · [GitHub](https://github.com/arjunbioinfo)
