FROM rocker/rstudio:4.5.0



ARG DIAREPORT_VERSION=0.9.1
ENV DIAREPORT_VERSION=$DIAREPORT_VERSION
LABEL version=$DIAREPORT_VERSION

LABEL org.opencontainers.image.title="Gevaert-Lab diareport"
LABEL org.opencontainers.image.authors="Andrea Argentini <aargentini@gmail.com>"
LABEL org.opencontainers.image.description="RStudio environment for diareport analysis with Chrome/Quarto support"
LABEL org.opencontainers.image.url="https://github.com/Gevaert-Lab/diareport"

# --------------------------------------------------
# Basic utilities & Chromium/Chrome dependencies
# --------------------------------------------------
# Note: libasound2t64 is the modern replacement for libasound2
RUN apt-get update && apt-get install -y --fix-missing \
    wget \
    unzip \
    default-jre \
    gdebi-core \
    libfontconfig1 \
    libnss3 \
    libatk-bridge2.0-0 \
    libgtk-3-0 \
    libxss1 \
    libasound2t64 \
    && rm -rf /var/lib/apt/lists/*

# --------------------------------------------------
# Development libraries
# --------------------------------------------------
RUN apt-get update && apt-get install -y --fix-missing \
    libcurl4-openssl-dev \
    libssl-dev \
    libssl3 \
    libxml2-dev \
    zlib1g-dev \
    libpng-dev \
    libjpeg-dev \
    libfreetype6-dev \
    libnetcdf-dev \
    libhdf5-dev \
    && rm -rf /var/lib/apt/lists/*

# --------------------------------------------------
# Numerical and optimization libraries
# --------------------------------------------------
RUN apt-get update && apt-get install -y --fix-missing \
    libglpk40 \
    libglpk-dev \
    libboost-all-dev \
    && rm -rf /var/lib/apt/lists/*

# --------------------------------------------------
# Node.js and npm
# --------------------------------------------------
RUN apt-get update && apt-get install -y --fix-missing \
    nodejs \
    npm \
    && rm -rf /var/lib/apt/lists/*

# --------------------------------------------------
# Quarto CLI v1.8.26
# --------------------------------------------------
RUN wget -qO quarto.deb \
    https://github.com/quarto-dev/quarto-cli/releases/download/v1.8.26/quarto-1.8.26-linux-amd64.deb && \
    apt-get update && apt-get install -y ./quarto.deb && \
    rm quarto.deb

# --------------------------------------------------
# Google Chrome (Bypassing Snap issues with Chromium)
# --------------------------------------------------
RUN wget -q https://dl.google.com/linux/direct/google-chrome-stable_current_amd64.deb && \
    apt-get update && \
    apt-get install -y ./google-chrome-stable_current_amd64.deb && \
    rm google-chrome-stable_current_amd64.deb

# --------------------------------------------------
# Environment variables
# --------------------------------------------------
ENV OPENSSL_CONF=/etc/ssl/openssl.cnf
ENV R_LIBS_SITE=/usr/local/lib/R/site-library
# Point R packages to the correct Chrome binary
ENV CHROMOTE_CHROME=/usr/bin/google-chrome

# --------------------------------------------------
# PhantomJS (for some R packages)
# --------------------------------------------------
RUN npm install -g phantomjs-prebuilt

# --------------------------------------------------
# R packages
# --------------------------------------------------
# Added missing closing quote and parenthesis at the end
RUN R -e "install.packages(c('BiocManager','remotes','ggrepel','tidyr','reactable','stringr','htmltools','purrr','ggVennDiagram','upsetjs','heatmaply','plotly','yaml','UpSetR','GGally'), repos='https://cloud.r-project.org', lib='/usr/local/lib/R/site-library'); \
          BiocManager::install(c('QFeatures','msqrob2','MSnbase'), lib='/usr/local/lib/R/site-library'); \
          remotes::install_github('Gevaert-Lab/diareport@${DIAREPORT_VERSION}', dependencies=TRUE)"

# --------------------------------------------------
# Verify Installations
# --------------------------------------------------
RUN quarto --version && google-chrome --version