#!/bin/bash
# submit.sh — driver for the CV sweep on Tide.
#
# Usage:
#   ./submit.sh upload   Create PVCs and upload sweep_all.csv + Kenny covariates
#   ./submit.sh submit   Create ConfigMap and submit one Job per (r_D, r_S)
#   ./submit.sh wait     Poll until all sweep jobs are complete
#   ./submit.sh fetch    Copy per-job outputs to _data/processed/sweep_cv/
#   ./submit.sh all      upload + submit + wait + fetch (full pipeline)
#   ./submit.sh clean    Delete sweep Jobs (PVCs preserved for reruns)
#
# This is the SECOND-PHASE sweep: it reads the per-(r_D, r_S) raw-metric
# CSVs that the sweep_process produces (collected into sweep_all.csv by
# scripts/tide/sweep_process/collect.py), runs spatial-block + LOFO CV for
# each combination in parallel across cores, and writes per-job CV-score
# CSVs that collect.R then assembles into one results table.
#
# Edit the sweep grid below to match the sweep_process grid.

set -euo pipefail

# ----- Configuration ----------------------------------------------------------

NAMESPACE="cal-poly-ruiz"
IMAGE="ghcr.io/ruizt/ssdd-r:latest"
CONFIGMAP="ssdd-cv-script"
DATA_PVC="ssdd-cv-data"
OUTPUT_PVC="ssdd-cv-output"
ACCESSOR_POD="ssdd-cv-accessor"

# Sweep grid — should match scripts/tide/sweep_process/submit.sh
FIRES=(eaton palisades)            # kept for the data layout
R_D_VALUES=(150 175 200 250 300 400)
R_S_VALUES=(10 25 50)

# Resources per job. mclapply uses CORES forks; ranger inside each is
# single-threaded, so request matches the inner parallelism.
CORES="4"
MEM="6Gi"

# Resolve paths from script location so this works from any CWD.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
COMPUTE_R="${SCRIPT_DIR}/compute.R"
LOCAL_SWEEP_CSV="${REPO_ROOT}/_data/processed/sweep/sweep_all.csv"
LOCAL_PROCESSED="${REPO_ROOT}/_data/processed"
LOCAL_OUT_DIR="${REPO_ROOT}/_data/processed/sweep_cv"

# ----- Helpers ----------------------------------------------------------------

ensure_pvcs() {
    kubectl apply -n "${NAMESPACE}" -f - <<EOF
apiVersion: v1
kind: PersistentVolumeClaim
metadata:
  name: ${DATA_PVC}
  namespace: ${NAMESPACE}
spec:
  storageClassName: rook-cephfs-tide
  accessModes: [ReadWriteMany]
  resources:
    requests: {storage: 2Gi}
---
apiVersion: v1
kind: PersistentVolumeClaim
metadata:
  name: ${OUTPUT_PVC}
  namespace: ${NAMESPACE}
spec:
  storageClassName: rook-cephfs-tide
  accessModes: [ReadWriteMany]
  resources:
    requests: {storage: 1Gi}
EOF
}

start_accessor() {
    kubectl apply -n "${NAMESPACE}" -f - <<EOF
apiVersion: v1
kind: Pod
metadata:
  name: ${ACCESSOR_POD}
  namespace: ${NAMESPACE}
spec:
  restartPolicy: Never
  containers:
    - name: accessor
      image: busybox
      command: ["sleep", "600"]
      volumeMounts:
        - {name: data,   mountPath: /data}
        - {name: output, mountPath: /jobs/output}
  volumes:
    - name: data
      persistentVolumeClaim: {claimName: ${DATA_PVC}}
    - name: output
      persistentVolumeClaim: {claimName: ${OUTPUT_PVC}}
EOF
    kubectl wait -n "${NAMESPACE}" --for=condition=Ready "pod/${ACCESSOR_POD}" --timeout=120s
}

stop_accessor() {
    kubectl delete -n "${NAMESPACE}" pod "${ACCESSOR_POD}" --ignore-not-found
}

# ----- Subcommands ------------------------------------------------------------

cmd_upload() {
    echo "=== upload: PVCs + sweep_all.csv + Kenny covariates ==="
    if [[ ! -f "${LOCAL_SWEEP_CSV}" ]]; then
        echo "ERROR: ${LOCAL_SWEEP_CSV} not found." >&2
        echo "  Run sweep_process first (./scripts/tide/sweep_process/submit.sh all)," >&2
        echo "  then assemble (python scripts/tide/sweep_process/collect.py)." >&2
        exit 1
    fi
    for fire in "${FIRES[@]}"; do
        local d="${LOCAL_PROCESSED}/${fire}/covariates"
        if [[ ! -d "${d}" ]]; then
            echo "ERROR: ${d} not found." >&2
            exit 1
        fi
    done

    ensure_pvcs
    start_accessor

    # Layout we want inside the PVC, mirroring local _data/processed shape:
    #   /data/sweep/sweep_all.csv
    #   /data/<fire>/covariates/...
    kubectl exec -n "${NAMESPACE}" "${ACCESSOR_POD}" -- mkdir -p \
        /data/sweep /data/eaton/covariates /data/palisades/covariates

    kubectl cp -n "${NAMESPACE}" "${LOCAL_SWEEP_CSV}" \
        "${ACCESSOR_POD}:/data/sweep/sweep_all.csv"

    for fire in "${FIRES[@]}"; do
        kubectl cp -n "${NAMESPACE}" \
            "${LOCAL_PROCESSED}/${fire}/covariates/." \
            "${ACCESSOR_POD}:/data/${fire}/covariates/"
    done

    kubectl exec -n "${NAMESPACE}" "${ACCESSOR_POD}" -- ls -la \
        /data/sweep /data/eaton/covariates /data/palisades/covariates

    stop_accessor
}

cmd_submit() {
    echo "=== submit: ConfigMap + Jobs ==="
    ensure_pvcs

    kubectl create configmap "${CONFIGMAP}" -n "${NAMESPACE}" \
        --from-file=compute.R="${COMPUTE_R}" \
        --dry-run=client -o yaml | kubectl apply -f -

    local n=0
    for r_d in "${R_D_VALUES[@]}"; do
        for r_s in "${R_S_VALUES[@]}"; do
            local job_name="ssdd-cv-rd${r_d}-rs${r_s}"
            kubectl apply -n "${NAMESPACE}" -f - <<EOF
apiVersion: batch/v1
kind: Job
metadata:
  name: ${job_name}
  namespace: ${NAMESPACE}
  labels: {app: ssdd-cv, rd: "${r_d}", rs: "${r_s}"}
spec:
  backoffLimit: 3
  template:
    spec:
      restartPolicy: Never
      containers:
        - name: ssdd-cv
          image: ${IMAGE}
          imagePullPolicy: IfNotPresent
          resources:
            # Cluster policy requires limit/request <= 1.2 — keep them equal.
            requests: {cpu: "${CORES}", memory: "${MEM}"}
            limits:   {cpu: "${CORES}", memory: "${MEM}"}
          env:
            - {name: SSDD_R_D,   value: "${r_d}"}
            - {name: SSDD_R_S,   value: "${r_s}"}
            - {name: SSDD_CORES, value: "${CORES}"}
          volumeMounts:
            - {name: script, mountPath: /scripts}
            - {name: data,   mountPath: /data, readOnly: true}
            - {name: output, mountPath: /jobs/output}
      volumes:
        - name: script
          configMap: {name: ${CONFIGMAP}}
        - name: data
          persistentVolumeClaim: {claimName: ${DATA_PVC}}
        - name: output
          persistentVolumeClaim: {claimName: ${OUTPUT_PVC}}
EOF
            n=$((n + 1))
        done
    done
    echo "Submitted ${n} jobs."
}

cmd_wait() {
    echo "=== wait: polling for completion ==="
    while true; do
        local total complete
        total=$(kubectl get jobs -n "${NAMESPACE}" -l app=ssdd-cv --no-headers 2>/dev/null | wc -l | tr -d ' ')
        complete=$(kubectl get jobs -n "${NAMESPACE}" -l app=ssdd-cv --no-headers 2>/dev/null | grep -c "1/1" || true)
        echo "  ${complete}/${total} complete"
        if [[ "${complete}" -eq "${total}" && "${total}" -gt 0 ]]; then
            break
        fi
        sleep 30
    done
}

cmd_fetch() {
    echo "=== fetch: outputs -> ${LOCAL_OUT_DIR} ==="
    mkdir -p "${LOCAL_OUT_DIR}"
    start_accessor
    kubectl cp -n "${NAMESPACE}" "${ACCESSOR_POD}:/jobs/output/." "${LOCAL_OUT_DIR}/"
    stop_accessor
    echo "Done. Assemble with: Rscript ${SCRIPT_DIR}/collect.R"
}

cmd_clean() {
    echo "=== clean: deleting CV-sweep Jobs (PVCs preserved) ==="
    kubectl delete jobs -n "${NAMESPACE}" -l app=ssdd-cv --ignore-not-found
}

cmd_all() {
    cmd_upload
    cmd_submit
    cmd_wait
    cmd_fetch
}

# ----- Dispatch ---------------------------------------------------------------

case "${1:-}" in
    upload) cmd_upload ;;
    submit) cmd_submit ;;
    wait)   cmd_wait   ;;
    fetch)  cmd_fetch  ;;
    clean)  cmd_clean  ;;
    all)    cmd_all    ;;
    *)
        echo "Usage: $0 {upload|submit|wait|fetch|all|clean}" >&2
        exit 1
        ;;
esac
