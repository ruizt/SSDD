## ssdd CV-sweep image
##
## R 4.5 + the packages needed by compute.R (dplyr, readr, sf, blockCV,
## ranger, yardstick). The script is NOT baked in — it's mounted at
## /scripts/compute.R via a Kubernetes ConfigMap at job-submit time, so a
## tweak to the R logic doesn't require an image rebuild.
##
## Build and push from the repo root:
##   docker buildx build --platform linux/amd64 \
##     -f scripts/tide/sweep_cv/Dockerfile \
##     -t ghcr.io/ruizt/ssdd-r:latest --push .
##
## Run locally:
##   docker run --rm \
##     -e SSDD_R_D=200 -e SSDD_R_S=50 \
##     -v $(pwd)/_data/processed:/data:ro \
##     -v $(pwd)/_tmp/cv:/jobs/output \
##     -v $(pwd)/scripts/tide/sweep_cv/compute.R:/scripts/compute.R:ro \
##     ghcr.io/ruizt/ssdd-r:latest

FROM rocker/r-ver:4.5.0

# System libraries required by sf / ranger / yardstick
RUN apt-get update && apt-get install -y --no-install-recommends \
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev \
        libgdal-dev \
        libgeos-dev \
        libproj-dev \
        libudunits2-dev \
        libfontconfig1-dev \
        && rm -rf /var/lib/apt/lists/*

# Install R packages from Posit Public Package Manager (Linux binaries —
# much faster + more reliable than source installs).
#
# PPM's binary repo is keyed on base-OS codename: packages built for
# Jammy link against libproj.so.22, packages built for Noble against
# libproj.so.25, etc. Hardcoding the codename breaks the image whenever
# rocker bumps its base, so source /etc/os-release at build time and let
# the shell substitute the right codename in.
RUN . /etc/os-release && \
    Rscript -e "install.packages( \
        c('dplyr', 'readr', 'sf', 'blockCV', 'ranger', 'yardstick'), \
        repos = 'https://packagemanager.posit.co/cran/__linux__/${VERSION_CODENAME}/latest')"

# Default: run whatever script is mounted at /scripts/compute.R.
CMD ["Rscript", "/scripts/compute.R"]
