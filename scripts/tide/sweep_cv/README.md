# SSDD CV sweep — Tide batch infrastructure

Phase 2 of the sweep. Reads `sweep_all.csv` (per-`(fire, r_D)` metrics from
`sweep_process`), NN-joins terrain covariates, and for each `r_D` evaluates every
metric **form-setting** (4 SD × 3 SS = 12, one `{SD, SS} + terrain` RF each) at
**two scales**:

- **structure** — RF classifier → `destroyed`
- **neighbourhood** — aggregate-then-fit: fine-block means → block damage-rate
  (RF regression)

CV is **one pooled spatial-block scheme** (no LOFO): each fold holds out a block
from every fire, so held-out performance is **disaggregated by fire**. Coarse
blocks define folds (`SSDD_CV_BLOCK`, default 1500 m); a separate fine grid
(`SSDD_NBHD_BLOCK`, default 400 m) is the neighbourhood unit. Per-fire within-fire
CV runs for the shipped form only, as a fire-heterogeneity diagnostic.

**Outer parallelism** = one Kubernetes Job per `r_D`. **Inner** = `mclapply`
forks across the fits. Output: `rD<r_D>_cv.csv` keyed
`(r_D, scale, cv, form, fire, fold)`.

## Files

| File | Role |
|---|---|
| `Dockerfile` | Build `ghcr.io/ruizt/ssdd-r:latest` (R 4.5 + dplyr/sf/blockCV/ranger/yardstick). |
| `compute.R` | Job entrypoint — reads `SSDD_R_D`, `SSDD_R_S`, `SSDD_CORES`; NN-joins Kenny terrain; runs CV in parallel. |
| `submit.sh` | Unified driver. Subcommands: `upload`, `submit`, `wait`, `fetch`, `all`, `clean`. |
| `collect.R` | Assemble per-job CSVs; also writes per-setting summary and per-setting optimum. |

## One-time setup

Build and push the R image:

```bash
docker buildx build --platform linux/amd64 \
  -f scripts/tide/sweep_cv/Dockerfile \
  -t ghcr.io/ruizt/ssdd-r:latest --push .
```

Make `ghcr.io/ruizt/ssdd-r` public on GitHub Packages so the cluster can
pull without a secret.

## Prerequisites for any run

The CV sweep consumes outputs from the radius sweep:

- `_data/processed/sweep/sweep_all.csv` (produced by `sweep_process/collect.py`)
- `_data/processed/{eaton,palisades}/covariates/*` (already in repo)

If you haven't run the radius sweep + collect yet:

```bash
./scripts/tide/sweep_process/submit.sh all
python scripts/tide/sweep_process/collect.py
```

## Running the CV sweep

End-to-end:

```bash
./scripts/tide/sweep_cv/submit.sh all
```

That uploads inputs, submits **6 × 3 = 18 Jobs** (matching the radius sweep
grid), polls until done, and fetches per-job CSVs to
`_data/processed/sweep_cv/`. Then:

```bash
Rscript scripts/tide/sweep_cv/collect.R
```

assembles three CSVs in `_data/processed/sweep_cv/`:

| File | Contents |
|---|---|
| `sweep_results.csv` | Long-format: one row per (r_D, r_S, setting, fold). |
| `sweep_summary.csv` | Mean ± SD per (r_D, r_S, setting). |
| `sweep_optima.csv`  | Best (r_D, r_S) per setting on mean AUC. |

The `sweep_optima.csv` is the small summary you usually want for a quick
answer about "which radii give the best held-out skill in each CV regime."

## Resources per Job

- **4 CPU** / **6 Gi memory** (request = limit, per the cluster's
  `limit/request ≤ 1.2` policy)
- `mclapply(..., mc.cores = 4)` forks across the 32 CV fits
- ranger inside each fork is single-threaded

## Adjusting the sweep grid

Three arrays near the top of `submit.sh`:

```bash
FIRES=(eaton palisades)
R_D_VALUES=(150 175 200 250 300 400)
R_S_VALUES=(10 25 50)
```

Should match whatever `sweep_process/submit.sh` produced. If you re-ran
`sweep_process` with a different grid, update both before launching the CV
sweep.

## Test locally

Without the cluster, single combination via R:

```bash
SSDD_R_D=200 SSDD_R_S=50 SSDD_CORES=4 \
  SSDD_DATA_DIR=$(pwd)/_data/processed SSDD_OUT_DIR=$(pwd)/_tmp/cv \
  Rscript scripts/tide/sweep_cv/compute.R
```

Or via Docker (good for confirming the image's R env):

```bash
docker run --rm --platform linux/amd64 \
  -e SSDD_R_D=200 -e SSDD_R_S=50 -e SSDD_CORES=4 \
  -v "$(pwd)/_data/processed":/data:ro \
  -v "$(pwd)/_tmp/cv":/jobs/output \
  -v "$(pwd)/scripts/tide/sweep_cv/compute.R":/scripts/compute.R:ro \
  ghcr.io/ruizt/ssdd-r:latest
```
