#!/usr/bin/env bash
set -euo pipefail

module use /work/modules
module load miniforge
conda-init

conda activate snakebase
cd /work/repos/af3_screening

RUN_NAME="example"
TMP_ROOT="/tmp/af3_pipeline"
TMP_RUN="${TMP_ROOT}/${RUN_NAME}"

WORKDIR="/work/datadir/af3_runs/${RUN_NAME}"
CONFIG="${WORKDIR}/pipeline_inputs/ucloud_config.yaml"
PROFILE="/work/repos/af3_screening/profiles/ucloud"
LOG="${WORKDIR}/pipeline_inputs/pipeline.log"

cleanup() {
    local status=$?

    # Prevent the EXIT trap from running recursively.
    trap - EXIT

    echo "Cleaning temporary directory: $TMP_RUN"
    rm -rf "$TMP_RUN" || true

    exit "$status"
}
trap cleanup EXIT

# Preserve TMP_ROOT even if Snakemake cleans empty parent directories.
mkdir -p "$TMP_ROOT"
touch "$TMP_ROOT/.keep"

# Remove stale data from a previous run.
rm -rf "$TMP_RUN"

SNAKEMAKE_ARGS=(
    --cores 192
    --profile "$PROFILE"
    --configfile "$CONFIG"
    --directory "$WORKDIR"
    --logger snkmt
    --logger-snkmt-db "${WORKDIR}/../snkmt.db"
)

# Unlock only when explicitly requested:
#
#     UNLOCK=1 ./run_pipeline.sh
#
if [[ "${UNLOCK:-1}" == "1" ]]; then
    snakemake "${SNAKEMAKE_ARGS[@]}" --unlock
fi

: >"$LOG"

snakemake \
    --dry-run \
    "${SNAKEMAKE_ARGS[@]}" \
    2>&1 | tee -a "$LOG"

snakemake \
    "${SNAKEMAKE_ARGS[@]}" \
    2>&1 | tee -a "$LOG"
