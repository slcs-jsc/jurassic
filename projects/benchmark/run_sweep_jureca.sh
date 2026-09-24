#!/bin/bash
#SBATCH --account=slmet
#SBATCH --partition=dc-cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=03:00:00
#SBATCH --exclusive
#SBATCH --disable-perfparanoid
#SBATCH --job-name=e3_sweep
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/base.sh"

bench_init "problem_size"
bench_analyse_topology

CHANNEL_COUNTS_FILE="${JR_CONFIG_DIR}/channel_counts.txt"
GAS_SETS_DIR="${JR_CONFIG_DIR}/gas_sets"

THREADS="${THREADS:-${JR_PHYS_PER_SOCKET}}"
CORES="${CORES:-E:S0:${THREADS}}"
REP="${REP:-3}"
BATCH_SIZE="${BATCH_SIZE:-$(get_batch_size 2>/dev/null || echo 1)}"

# ND/NG are compile-time array bounds 
BIN_CACHE_DIR="${JR_BENCH_DIR}/bin_cache"
mkdir -p "${BIN_CACHE_DIR}"

# build_or_reuse ND NG
# Echoes the path to a formod binary compiled with -DND=<nd> -DNG=<ng>,
# building it once per distinct (nd, ng) pair and reusing it thereafter.
build_or_reuse() {
    local nd="$1" ng="$2"
    local key="nd${nd}_ng${ng}"
    local variant_dir="${BIN_CACHE_DIR}/${key}"
    local variant_bin="${variant_dir}/formod"

    if [[ -x "${variant_bin}" ]]; then
        echo "[e3] reusing cached build for ND=${nd} NG=${ng}" >&2
        echo "${variant_bin}"
        return
    fi

    echo "[e3] building ND=${nd} NG=${ng} -> ${variant_dir}" >&2
    mkdir -p "${variant_dir}"
    (
        cd "${JR_REPO_ROOT}/src"
        make clean
        make DEFINES="-DND=${nd} -DNG=${ng}"
    ) 1>&2
    cp "${JR_REPO_ROOT}/src/formod" "${variant_bin}"

    echo "${variant_bin}"
}

run_point() {
    local label="$1" nd="$2" ng="$3" gas_file="$4"

    local ctl_out="${JR_BENCH_DIR}/ctl/${label}.ctl"
    mkdir -p "$(dirname "${ctl_out}")"

    local gen_args=(--in "${JR_ACTIVE_CTL_BASE}" --out "${ctl_out}" --nd "${nd}")
    if [[ -n "${gas_file}" ]]; then
        gen_args+=(--gas-file "${gas_file}")
    fi
    python3 "${SCRIPT_DIR}/generate_ctl.py" "${gen_args[@]}"

    JR_ACTIVE_CTL="${ctl_out}"

    local formod_bin
    formod_bin="$(build_or_reuse "${nd}" "${ng}")"
    JR_FORMOD_BIN="${formod_bin}"

    bench_prepare_inputs
    bench_run_forward "${label}" "${THREADS}" "FLOPS_DP" "${BATCH_SIZE}" "${REP}" "${CORES}" ""
    bench_run_forward "${label}" "${THREADS}" "MEM_DP"   "${BATCH_SIZE}" "${REP}" "${CORES}" ""
}

JR_ACTIVE_CTL_BASE="${JR_ACTIVE_CTL}"

# NG held fixed at the baseline gas count (7, from core.txt) while ND varies;
# ND held fixed at the baseline channel count (32) while NG varies.
BASELINE_ND="${BASELINE_ND:-32}"
BASELINE_NG="${BASELINE_NG:-7}"

echo "=== channel scaling ==="
while read -r nd; do
    [[ -z "${nd}" ]] && continue
    run_point "channels_${nd}" "${nd}" "${BASELINE_NG}" ""
done < "${CHANNEL_COUNTS_FILE}"

echo "=== gas set scaling ==="
for gas_file in "${GAS_SETS_DIR}"/*.txt; do
    [[ -e "${gas_file}" ]] || continue
    name="$(basename "${gas_file}" .txt)"
    ng="$(grep -c . "${gas_file}")"
    run_point "gases_${name}" "${BASELINE_ND}" "${ng}" "${gas_file}"
done

bench_finish