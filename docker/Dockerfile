# =============================================================================
#  Dockerfile  --  pinned, reproducible tidymarkers image with the MCP server
# =============================================================================
#  Base is pinned to the R the lockfile targets (R 4.6.1).  Dependencies are
#  installed with `pak`, which pulls Posit Package Manager Linux binaries, so the
#  build does not compile the Bioconductor stack from source.  Only `Imports` are
#  installed: the heavy `Suggests` (annotation, rBLAST, Rsubread, rprimer,
#  crisprVerse) stay out of the image, and the package reports `unknown` for the
#  capabilities that need them instead of pretending they exist.
#
#  Build:  docker build -t tidymarkers:0.4.0 .
#  Run:    see docker/README.md
# =============================================================================
FROM rocker/r-ver:4.6.1

ARG DEBIAN_FRONTEND=noninteractive

# Runtime libraries for the Imports (curl/ssl/xml/icu/PNG/fonts) and the
# toolchain R CMD INSTALL needs.
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libicu-dev \
    libpng-dev \
    libfreetype6-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libbz2-dev \
    liblzma-dev \
    zlib1g-dev \
    libgit2-dev \
    make \
    g++ \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/tidymarkers
COPY . /opt/tidymarkers

# pak defaults to PPM and picks Linux binaries for this distro.
RUN R -q -e 'install.packages("pak", repos = "https://cloud.r-project.org")' \
    && R -q -e 'pak::local_install_dev_deps("/opt/tidymarkers", dependencies = c("Depends", "Imports", "LinkingTo"), ask = FALSE)'

RUN R CMD INSTALL --no-multiarch --no-test-load /opt/tidymarkers

RUN chmod +x docker/entrypoint.sh

# Stdio only: the server reads JSON-RPC on stdin and writes it on stdout, so the
# container opens no port.  Mount data read-only at /data (see docker/README.md).
ENTRYPOINT ["/opt/tidymarkers/docker/entrypoint.sh"]
