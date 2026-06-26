## ssdd radius-sweep image
##
## Contains Python 3.11 + the ssdd package (Imports only — geopandas, shapely,
## numpy, pandas, tqdm). Scripts and input data are NOT baked in:
##   - compute.py is mounted at /scripts/compute.py via a Kubernetes ConfigMap
##   - input shapefiles + DINS GeoJSONs are mounted at /data from a PVC
##   - per-job outputs are written to /jobs/output (also a PVC)
##
## Build and push from the repo root:
##   docker buildx build --platform linux/amd64 \
##     -f scripts/tide/sweep_process/Dockerfile \
##     -t ghcr.io/ruizt/ssdd:latest --push .
##
## Run locally (mount a script + data + output):
##   docker run --rm \
##     -e SSDD_FIRE=palisades -e SSDD_R_D=100 -e SSDD_R_S=50 \
##     -v $(pwd)/_data/raw:/data:ro \
##     -v $(pwd)/_tmp/sweep:/jobs/output \
##     -v $(pwd)/scripts/tide/sweep_process/compute.py:/scripts/compute.py:ro \
##     ghcr.io/ruizt/ssdd:latest

FROM python:3.11-slim

# Build-essential is needed only if any wheel falls back to source.
# geopandas/shapely/pyogrio ship manylinux wheels with bundled GEOS/GDAL/PROJ
# in their current 1.x / 2.x releases, so no system geospatial libs needed.
RUN apt-get update && apt-get install -y --no-install-recommends \
        gcc \
        && rm -rf /var/lib/apt/lists/*

# Copy the package source and install with its declared dependencies.
WORKDIR /build/ssdd
COPY . /build/ssdd
RUN pip install --no-cache-dir -e .

# Default: run whatever script is mounted at /scripts/compute.py.
CMD ["python", "/scripts/compute.py"]
