FROM rocker/r-ver:4.4.1

LABEL org.opencontainers.image.title="cancer-cohort-wgcna-survival" \
      org.opencontainers.image.description="Reproducible R/Bioconductor workflow for TCGA cancer cohort limma DE, WGCNA, fgsea, and Cox/Kaplan-Meier survival analysis" \
      org.opencontainers.image.authors="Mallikarjuna Thippana" \
      org.opencontainers.image.licenses="MIT"

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
      make libcurl4-openssl-dev libssl-dev libxml2-dev zlib1g-dev && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY DESCRIPTION ./
RUN Rscript -e "install.packages('pak', repos = 'https://cloud.r-project.org'); pak::local_install_deps(dependencies = TRUE)"

COPY run_pipeline.R Makefile README.md LICENSE CITATION.cff ./
COPY config ./config
COPY R ./R
COPY tests ./tests
RUN mkdir -p data/raw results figures

CMD ["Rscript", "run_pipeline.R", "--config", "config/ci.yaml"]
