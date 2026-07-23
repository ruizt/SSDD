#!/bin/bash
# submit.sh — driver for the SSDD radius sweep on Tide.
#
# Usage:
#   ./submit.sh upload    Create PVCs and upload _data/raw inputs to the data PVC
#   ./submit.sh submit    Create ConfigMap and submit one Job per (fire, r_D)
#   ./submit.sh wait      Poll until all sweep jobs are complete
#   ./submit.sh fetch     Copy per-job outputs to _data/processed/sweep/
#   ./submit.sh all       upload + submit + wait + fetch (full pipeline)
#   ./submit.sh clean     Delete sweep Jobs (PVCs preserved for reruns)
#
# Edit the sweep grid below to change which (fire, r_D) combinations run.
# Paths resolve from the script's location, so it works from any CWD.

set -euo pipefail

# ----- Configuration ----------------------------------------------------------

NAMESPACE="cal-poly-ruiz"
# Versioned, never :latest — jobs pull with IfNotPresent, so a floating tag
# would let a node that cached an older :latest serve stale package code.
# Bump this whenever the package source changes and you rebuild.
IMAGE="ghcr.io/ruizt/ssdd:v0.2"
CONFIGMAP="ssdd-sweep-script"
DATA_PVC="ssdd-sweep-data"
OUTPUT_PVC="ssdd-sweep-output"
ACCESSOR_POD="ssdd-sweep-accessor"

# Sweep grid — edit to taste. r_S is gone (SS uses the true nearest neighbour).
FIRES=(eaton palisades mountain)
R_D_VALUES=(50 100 150 200 250 300)

# Resolve paths from script location so this works from any CWD.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
COMPUTE_PY="${SCRIPT_DIR}/compute.py"
LOCAL_DATA_DIR="${REPO_ROOT}/_data/processed"
LOCAL_OUT_DIR="${REPO_ROOT}/_data/processed/sweep"

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
    requests: {storage: 1Gi}
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
    requests: {storage: 5Gi}
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
    echo "=== upload: PVCs + processed buildings gpkgs ==="
    for fire in "${FIRES[@]}"; do
        if [[ ! -f "${LOCAL_DATA_DIR}/${fire}/${fire}_buildings.gpkg" ]]; then
            echo "ERROR: expected ${LOCAL_DATA_DIR}/${fire}/${fire}_buildings.gpkg" >&2
            exit 1
        fi
    done
    ensure_pvcs
    start_accessor
    # compute.py reads /data/<fire>/<fire>_buildings.gpkg
    for fire in "${FIRES[@]}"; do
        kubectl exec -n "${NAMESPACE}" "${ACCESSOR_POD}" -- mkdir -p "/data/${fire}"
        kubectl cp -n "${NAMESPACE}" \
            "${LOCAL_DATA_DIR}/${fire}/${fire}_buildings.gpkg" \
            "${ACCESSOR_POD}:/data/${fire}/${fire}_buildings.gpkg"
    done
    kubectl exec -n "${NAMESPACE}" "${ACCESSOR_POD}" -- ls -la /data
    stop_accessor
}

cmd_submit() {
    echo "=== submit: ConfigMap + Jobs ==="
    ensure_pvcs

    kubectl create configmap "${CONFIGMAP}" -n "${NAMESPACE}" \
        --from-file=compute.py="${COMPUTE_PY}" \
        --dry-run=client -o yaml | kubectl apply -f -

    local n=0
    for fire in "${FIRES[@]}"; do
        for r_d in "${R_D_VALUES[@]}"; do
                local job_name="ssdd-sweep-${fire}-rd${r_d}"
                kubectl apply -n "${NAMESPACE}" -f - <<EOF
apiVersion: batch/v1
kind: Job
metadata:
  name: ${job_name}
  namespace: ${NAMESPACE}
  labels: {app: ssdd-sweep, fire: "${fire}", rd: "${r_d}"}
spec:
  # Allow a few retries — ghcr.io occasionally rate-limits anonymous pulls
  # and a stuck pull-pod can otherwise burn the only retry attempt.
  backoffLimit: 3
  template:
    spec:
      restartPolicy: Never
      containers:
        - name: ssdd-sweep
          image: ${IMAGE}
          # IfNotPresent: once a node has the image cached, subsequent pods
          # on that node skip the pull. Cuts registry pressure dramatically
          # across a multi-job sweep. (Switch back to Always if you push a
          # new image while a sweep is mid-flight.)
          imagePullPolicy: IfNotPresent
          resources:
            # Cluster policy requires limit/request <= 1.2 — keep them equal.
            requests: {cpu: "1", memory: "3Gi"}
            limits:   {cpu: "1", memory: "3Gi"}
          env:
            - {name: SSDD_FIRE, value: "${fire}"}
            - {name: SSDD_R_D,  value: "${r_d}"}
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
        total=$(kubectl get jobs -n "${NAMESPACE}" -l app=ssdd-sweep --no-headers 2>/dev/null | wc -l | tr -d ' ')
        complete=$(kubectl get jobs -n "${NAMESPACE}" -l app=ssdd-sweep --no-headers 2>/dev/null | grep -c "1/1" || true)
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
    echo "Done. Assemble with: python ${SCRIPT_DIR}/collect.py"
}

cmd_clean() {
    echo "=== clean: deleting sweep Jobs (PVCs preserved) ==="
    kubectl delete jobs -n "${NAMESPACE}" -l app=ssdd-sweep --ignore-not-found
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
