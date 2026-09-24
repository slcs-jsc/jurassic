#!/bin/bash
#SBATCH --account=slmet
#SBATCH --partition=dc-cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=03:00:00
#SBATCH --exclusive
#SBATCH --disable-perfparanoid
#SBATCH --job-name=e1_sweep
set -euo pipefail

JR_EXPERIMENT=e3_sweep
RUN_ID=${RUN_ID:-e3_sweep_${SLURM_JOB_ID:-manual}}

if [ -n "${SLURM_SUBMIT_DIR:-}" ] && [ -f "$SLURM_SUBMIT_DIR/base.sh" ]; then
  jr_scripts_dir="$SLURM_SUBMIT_DIR"
else
  jr_scripts_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi

BATCH_SIZE=${BATCH_SIZE:-1024}

export JR_SCRIPTS_DIR_OVERRIDE="$jr_scripts_dir"
source "$jr_scripts_dir/base.sh"

SCRIPT_DIR="$jr_scripts_dir"

if ! grep -q 'JR_FORMOD_BIN' "$jr_scripts_dir/base.sh"; then
  echo "ERROR: base.sh does not support JR_FORMOD_BIN yet." >&2
  echo "Apply base.sh.patch.txt (bench_run_forward: use \${JR_FORMOD_BIN:-\$JR_SRC_DIR/formod})" >&2
  echo "before running this sweep, otherwise every point silently reuses whichever" >&2
  echo "binary is currently at \$JR_SRC_DIR/formod." >&2
  exit 1
fi

bench_init

CONFIG_DIR="$JR_REPO_ROOT/projects/benchmark/configs"
CHANNEL_COUNTS_FILE="$CONFIG_DIR/channel_counts.txt"
GAS_SETS_DIR="$CONFIG_DIR/gas_sets"

THREADS="${THREADS:-${JR_PHYS_PER_SOCKET}}"
CORES="${CORES:-E:S0:${THREADS}}"
REP="${REP:-3}"

# ND/NG are compile-time array bounds (see src/jurassic.h)
BIN_CACHE_DIR="$JR_REPO_ROOT/projects/benchmark/bin_cache"
BUILD_SCRATCH_DIR="$BIN_CACHE_DIR/_build"
mkdir -p "$BIN_CACHE_DIR" "$BUILD_SCRATCH_DIR"

# build_or_reuse ND NG
build_or_reuse() {
    local nd="$1" ng="$2"
    local key="nd${nd}_ng${ng}"
    local variant_dir="$BIN_CACHE_DIR/${key}"
    local variant_bin="$variant_dir/formod"

    if [[ -x "$variant_bin" ]]; then
        echo "[e3] reusing cached build for ND=${nd} NG=${ng}" >&2
        echo "$variant_bin"
        return
    fi

    local lock_dir="$BUILD_SCRATCH_DIR/${key}.lock"
    local waited=0
    while ! mkdir "$lock_dir" 2>/dev/null; do
        sleep 5
        waited=$((waited + 5))
        if [[ -x "$variant_bin" ]]; then
            echo "[e3] variant ND=${nd} NG=${ng} built by another process" >&2
            echo "$variant_bin"
            return
        fi
        if (( waited > 1800 )); then
            echo "[e3] ERROR: timed out waiting for lock $lock_dir" >&2
            exit 1
        fi
    done
    trap 'rmdir "'"$lock_dir"'" 2>/dev/null' RETURN EXIT

    if [[ -x "$variant_bin" ]]; then
        echo "[e3] reusing cached build for ND=${nd} NG=${ng}" >&2
        echo "$variant_bin"
        return
    fi

    local build_dir="$BUILD_SCRATCH_DIR/${key}"
    echo "[e3] building ND=${nd} NG=${ng} in private copy -> $build_dir" >&2
    rm -rf "$build_dir"
    cp -r "$JR_SRC_DIR" "$build_dir"
    (
        cd "$build_dir"
        make clean
        make -j MPI="$JR_MPI" MPICC="$JR_MPICC" COMPILER="$JR_COMPILER" \
            GPU=0 LIKWID=1 DEFINES="-DND=${nd} -DNG=${ng}"
    ) 1>&2

    mkdir -p "$variant_dir"
    cp "$build_dir/formod" "$variant_bin"
    rm -rf "$build_dir"

    echo "$variant_bin"
}

run_point() {
    local label="$1" nd="$2" ng="$3" gas_file="$4"

    local ctl_out="$JR_WORK_DIR/ctl/${label}.ctl"
    mkdir -p "$(dirname "$ctl_out")"

    local gen_args=(--in "$JR_ACTIVE_CTL_BASE" --out "$ctl_out" --nd "$nd")
    if [[ -n "$gas_file" ]]; then
        gen_args+=(--gas-file "$gas_file")
    fi
    python3 "$SCRIPT_DIR/generate_ctl.py" "${gen_args[@]}"

    JR_ACTIVE_CTL="$ctl_out"
    JR_FORMOD_BIN="$(build_or_reuse "$nd" "$ng")"
    export JR_FORMOD_BIN

    bench_prepare_inputs
    bench_run_forward "$label" "$THREADS" "FLOPS_DP" "$BATCH_SIZE" "$REP" "$CORES" ""
    bench_run_forward "$label" "$THREADS" "MEM_DP"   "$BATCH_SIZE" "$REP" "$CORES" ""
}

JR_ACTIVE_CTL_BASE="$JR_ACTIVE_CTL"

# NG held fixed at the selected baseline case's gas count while ND varies;
# ND held fixed at the baseline channel count while NG varies. 
BASELINE_ND="${BASELINE_ND:-32}"
BASELINE_NG="${BASELINE_NG:-7}"

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