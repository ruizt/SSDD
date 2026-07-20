# SSDD radius sweep — Tide batch infrastructure

Phase 1 of the sweep. For each `(fire, r_D)` it computes the **candidate metric
forms** the CV sweep compares — 4 SD (`{uniform, quartic}` kernel × `{root_area,
unit}` weight) and 3 SS (nearest-neighbour, `orient ∈ {flat, gauss, cos2}`) — in
one vectorised pass, and writes a per-job CSV with the 7 metric columns plus
`ssdd_id`, `cent_x/y`, `DAMAGE`. Fires: `eaton palisades mountain`.

`r_S` is gone: SS uses the true nearest neighbour (no isolation cutoff), so the
only metric-radius dimension is `r_D`. Input is each fire's **processed** gpkg
(`_data/processed/<fire>/<fire>_buildings.gpkg`, DINS already joined), staged to
`/data/<fire>/`. Each `(fire, r_D)` is one independent Kubernetes Job.

## Files

| File | Role |
|---|---|
| `Dockerfile` | Build the `ghcr.io/ruizt/ssdd:latest` image (Python + the ssdd package). |
| `compute.py` | Job entrypoint — mounted via ConfigMap; reads `SSDD_*` env vars and calls `compute_raw_metrics`. |
| `submit.sh` | Unified driver. Subcommands: `upload`, `submit`, `wait`, `fetch`, `all`, `clean`. All PVC, accessor-pod, and Job manifests are inline heredocs. |
| `collect.py` | Concatenate per-job CSVs into one long-format `sweep_all.csv`. |

## One-time setup

Build and push the image (only needed when the package or its dependencies
change):

```bash
docker buildx build --platform linux/amd64 \
  -f scripts/tide/sweep_process/Dockerfile \
  -t ghcr.io/ruizt/ssdd:latest --push .
```

After the first push, set the package visibility to public on GitHub so the
cluster can pull without a registry secret (Settings → Packages → Change
visibility).

## Running the sweep

End-to-end:

```bash
./scripts/tide/sweep_process/submit.sh all
```

That runs `upload` → `submit` → `wait` → `fetch` in sequence. Or run the
subcommands individually:

```bash
./scripts/tide/sweep_process/submit.sh upload   # creates PVCs, kubectl cp's _data/raw into the data PVC
./scripts/tide/sweep_process/submit.sh submit   # creates ConfigMap, submits 32 Jobs
./scripts/tide/sweep_process/submit.sh wait     # polls every 30 s until all Jobs are 1/1 complete
./scripts/tide/sweep_process/submit.sh fetch    # kubectl cp's results into _data/processed/sweep/
```

Assemble per-job CSVs into a single table:

```bash
python scripts/tide/sweep_process/collect.py
```

Cleanup when done (PVCs preserved so reruns are quick):

```bash
./scripts/tide/sweep_process/submit.sh clean    # deletes only the Jobs
```

To wipe the cluster state entirely, also remove the PVCs by name:

```bash
kubectl delete pvc -n cal-poly-ruiz ssdd-sweep-data ssdd-sweep-output
```

## Adjusting the sweep grid

Three arrays near the top of `submit.sh`:

```bash
FIRES=(eaton palisades)
R_D_VALUES=(50 100 150 200)
R_S_VALUES=(25 50 75 100)
```

Defaults yield 2 × 4 × 4 = 32 independent Jobs. The cluster runs them in
parallel up to your namespace's resource quota.

## Test locally before pushing to the cluster

A single combination, no Kubernetes involved:

```bash
SSDD_FIRE=palisades SSDD_R_D=100 SSDD_R_S=50 \
  SSDD_DATA_DIR=$(pwd)/_data/raw SSDD_OUT_DIR=$(pwd)/_tmp \
  python scripts/tide/sweep_process/compute.py
```

Or with the Docker image:

```bash
docker run --rm --platform linux/amd64 \
  -e SSDD_FIRE=palisades -e SSDD_R_D=100 -e SSDD_R_S=50 \
  -v "$(pwd)/_data/raw":/data:ro \
  -v "$(pwd)/_tmp/sweep":/jobs/output \
  -v "$(pwd)/scripts/tide/sweep_process/compute.py":/scripts/compute.py:ro \
  ghcr.io/ruizt/ssdd:latest
```

(Note: Docker Desktop on macOS sometimes raises EDEADLK when reading
shapefile sidecars over the bind mount. That's a host-side quirk and does
not happen on the cluster, where data is read from a real PVC.)

## Configuration reference

Constants at the top of `submit.sh` if you want to retarget:

| Variable | Default | Meaning |
|---|---|---|
| `NAMESPACE` | `cal-poly-ruiz` | Kubernetes namespace |
| `IMAGE` | `ghcr.io/ruizt/ssdd:latest` | Image jobs pull |
| `CONFIGMAP` | `ssdd-sweep-script` | Name of the ConfigMap holding `compute.py` |
| `DATA_PVC` | `ssdd-sweep-data` | PVC for input shapefiles + DINS |
| `OUTPUT_PVC` | `ssdd-sweep-output` | PVC for per-job CSVs |
| `ACCESSOR_POD` | `ssdd-sweep-accessor` | Sleep-pod name used by `upload` and `fetch` |
