#!/bin/bash
# submit.sh — driver for the SSDD radius sweep on Tide.
#
# Usage:
#   ./submit.sh upload    Create PVCs and upload _data/raw inputs to the data PVC
#   ./submit.sh submit    Create ConfigMap and submit one Job per (fire, r_D, r_S)
#   ./submit.sh wait      Poll until all sweep jobs are complete
#   ./submit.sh fetch     Copy per-job outputs to _data/processed/sweep/
#   ./submit.sh all       upload + submit + wait + fetch (full pipeline)
#   ./submit.sh clean     Delete sweep Jobs (PVCs preserved for reruns)
#
# Edit the sweep grid below to change which (fire, r_D, r_S) combinations run.
# Paths resolve from the script's location, so it works from any CWD.

set -euo pipefail

# ----- Configuration ----------------------------------------------------------

NAMESPACE="cal-poly-ruiz"
IMAGE="ghcr.io/ruizt/ssdd:latest"
CONFIGMAP="ssdd-sweep-script"
DATA_PVC="ssdd-sweep-data"
OUTPUT_PVC="ssdd-sweep-output"
ACCESSOR_POD="ssdd-sweep-accessor"

# Sweep grid — edit to taste.
FIRES=(eaton palisades)
R_D_VALUES=(50 100 150 200)
R_S_VALUES=(25 50 75 100)
R_NN="200"

# Resolve paths from script location so this works from any CWD.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
COMPUTE_PY="${SCRIPT_DIR}/compute.py"
LOCAL_DATA_DIR="${REPO_ROOT}/_data/raw"
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
    echo "=== upload: PVCs + input data ==="
    if [[ ! -d "${LOCAL_DATA_DIR}/buildings" || ! -d "${LOCAL_DATA_DIR}/dins" ]]; then
        echo "ERROR: expected ${LOCAL_DATA_DIR}/{buildings,dins}" >&2
        exit 1
    fi
    ensure_pvcs
    start_accessor
    kubectl cp -n "${NAMESPACE}" "${LOCAL_DATA_DIR}/." "${ACCESSOR_POD}:/data/"
    kubectl exec -n "${NAMESPACE}" "${ACCESSOR_POD}" -- ls -la /data /data/buildings /data/dins
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
            for r_s in "${R_S_VALUES[@]}"; do
                local job_name="ssdd-sweep-${fire}-rd${r_d}-rs${r_s}"
                kubectl apply -n "${NAMESPACE}" -f - <<EOF
apiVersion: batch/v1
kind: Job
metadata:
  name: ${job_name}
  namespace: ${NAMESPACE}
  labels: {app: ssdd-sweep, fire: "${fire}", rd: "${r_d}", rs: "${r_s}"}
spec:
  backoffLimit: 1
  template:
    spec:
      restartPolicy: Never
      containers:
        - name: ssdd-sweep
          image: ${IMAGE}
          imagePullPolicy: Always
          resources:
            # Cluster policy requires limit/request <= 1.2 — keep them equal.
            requests: {cpu: "1", memory: "3Gi"}
            limits:   {cpu: "1", memory: "3Gi"}
          env:
            - {name: SSDD_FIRE, value: "${fire}"}
            - {name: SSDD_R_D,  value: "${r_d}"}
            - {name: SSDD_R_S,  value: "${r_s}"}
            - {name: SSDD_R_NN, value: "${R_NN}"}
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
