FROM bioconductor/bioconductor_docker:3.23

RUN apt-get update && \
    apt-get install -y --no-install-recommends zip libcairo2 libfontconfig1 fonts-dejavu-core && \
    rm -rf /var/lib/apt/lists/*

# CRAN + Bioconductor packages required by the app
RUN R -e "install.packages(c('shiny', 'ggplot2', 'viridis', 'plotly', 'dplyr', 'here'))" && \
    R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager'); \
          BiocManager::install('SingleCellExperiment', ask = FALSE, update = FALSE)"

WORKDIR /srv/app
RUN mkdir ./data/
COPY app.R ./

EXPOSE 3838

CMD ["R", "-e", "shiny::runApp('.', host = '0.0.0.0', port = 3838)"]
