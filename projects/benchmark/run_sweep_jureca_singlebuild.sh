#!/bin/bash
#SBATCH --account=slmet
#SBATCH --partition=dc-cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=02:00:00
#SBATCH --exclusive
#SBATCH --disable-perfparanoid
#SBATCH --job-name=e3_sweep_singlebuild
set -euo pipefail
shopt -s inherit_errexit

# Single-build counterpart of run_sweep_jureca.sh: that script rebuilds
# formod/climatology/<geometry> once per (ND,NG) point so every point's
# compile-time array bounds exactly match its .ctl, at the cost of one
# private "make clean && make" per distinct (ND,NG) pair in the sweep. This
# script instead builds exactly once, sized for the LARGEST ND and NG
# anywhere in the sweep's config files (channel_counts.txt and gas_sets/),
# and reuses that one binary set for every point.
#
# Tradeoff: fewer builds (one instead of one per distinct ND/NG), but every
# point other than the single largest one runs with compile-time array
# bounds bigger than its own .ctl needs -- los_t/atm_t (src/jurassic.h) are
# sized [[ND]]/[NLOS][ND] etc., so this inflates per-ray scratch memory
# (cache footprint, memset/fill cost) uniformly across all points, which
# will shift absolute timings somewhat relative to run_sweep_jureca.sh's
# per-point-sized binaries. Fine for a quick/cheap pass; use
# run_sweep_jureca.sh instead when the per-point memory footprint itself is
# part of what's being measured.

JR_EXPERIMENT=e3_sweep_singlebuild
RUN_ID=${RUN_ID:-e3_sweep_singlebuild_${SLURM_JOB_ID:-manual}}

if [ -n "${SLURM_SUBMIT_DIR:-}" ] && [ -f "$SLURM_SUBMIT_DIR/base.sh" ]; then
  jr_scripts_dir="$SLURM_SUBMIT_DIR"
else
  jr_scripts_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi

BATCH_SIZE=${BATCH_SIZE:-64}

export JR_SCRIPTS_DIR_OVERRIDE="$jr_scripts_dir"
source "$jr_scripts_dir/base.sh"

SCRIPT_DIR="$jr_scripts_dir"

bench_init

CONFIG_DIR="$JR_REPO_ROOT/projects/benchmark/configs"
CHANNEL_COUNTS_FILE="$CONFIG_DIR/channel_counts.txt"
GAS_SETS_DIR="$CONFIG_DIR/gas_sets"

THREADS="${THREADS:-${JR_PHYS_PER_SOCKET}}"
CORES="${CORES:-$(cpus_phys "$THREADS")}"
REP="${REP:-3}"

# NG held fixed at the selected baseline case's gas count while ND varies;
# ND held fixed at the baseline channel count while NG varies -- same
# baseline semantics as run_sweep_jureca.sh.
BASELINE_ND="${BASELINE_ND:-32}"
BASELINE_NG="${BASELINE_NG:-7}"

# Largest ND across the channel-scaling sweep (or BASELINE_ND if that's bigger).
MAX_ND=$(awk -v base="$BASELINE_ND" '
  /./ { if ($1 > max || NR == 1) max = $1 }
  END { if (base > max) max = base; print max }
' "$CHANNEL_COUNTS_FILE")

# Largest NG across the gas-set-scaling sweep (or BASELINE_NG if that's bigger).
MAX_NG="$BASELINE_NG"
for gas_file in "$GAS_SETS_DIR"/*.txt; do
  [[ -e "$gas_file" ]] || continue
  ng="$(grep -c . "$gas_file")"
  (( ng > MAX_NG )) && MAX_NG="$ng"
done

echo "=== single build sized for ND=${MAX_ND} NG=${MAX_NG} (covers every point below) ===" >&2
bench_build_isolated MPI="$JR_MPI" MPICC="$JR_MPICC" COMPILER="$JR_COMPILER" \
  GPU=0 LIKWID=1 DEFINES="-DND=${MAX_ND} -DNG=${MAX_NG}" || exit 1
echo "max_nd=$MAX_ND" >> "$JR_RUN_DIR/config.txt"
echo "max_ng=$MAX_NG" >> "$JR_RUN_DIR/config.txt"

# JR_BIN_DIR is now set (by bench_build_isolated) to the one build everything
# below uses. JR_FORMOD_BIN is deliberately left unset so bench_run_forward
# falls back to "${JR_BIN_DIR:-$JR_SRC_DIR}/formod" -- the single build.
run_point() {
    local label="$1" nd="$2" ng="$3" gas_file="$4"

    local ctl_out="$JR_WORK_DIR/ctl/${label}.ctl"
    mkdir -p "$(dirname "$ctl_out")"

    local gen_args=(--nd "$nd")
    if [[ -n "$gas_file" ]]; then
        gen_args+=(--gas-file "$gas_file")
    fi
    python3 "$SCRIPT_DIR/generate_ctl.py" "${gen_args[@]}" "$JR_ACTIVE_CTL_BASE" "$ctl_out"

    JR_ACTIVE_CTL="$ctl_out"

    bench_prepare_inputs
    bench_run_forward "$label" "$THREADS" "FLOPS_DP" "$BATCH_SIZE" "$REP" "$CORES" ""
    bench_run_forward "$label" "$THREADS" "MEM_DP"   "$BATCH_SIZE" "$REP" "$CORES" ""
}

JR_ACTIVE_CTL_BASE="$JR_ACTIVE_CTL"

echo "=== channel scaling ==="
while read -r nd; do
    [[ -z "$nd" ]] && continue
    run_point "channels_${nd}" "$nd" "$BASELINE_NG" ""
done < "$CHANNEL_COUNTS_FILE"

echo "=== gas set scaling ==="
for gas_file in "$GAS_SETS_DIR"/*.txt; do
    [[ -e "$gas_file" ]] || continue
    name="$(basename "$gas_file" .txt)"
    ng="$(grep -c . "$gas_file")"
    run_point "gases_${name}" "$BASELINE_ND" "$ng" "$gas_file"
done

cp -a "$JR_WORK_DIR/ctl" "$JR_RUN_DIR/" 2>/dev/null || true
bench_finish
