#!/bin/bash
#SBATCH --account=slmet
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --time=02:00:00
#SBATCH --exclusive
#SBATCH --disable-perfparanoid
#SBATCH --job-name=e3_sweep
set -euo pipefail
shopt -s inherit_errexit

JR_EXPERIMENT=e3_sweep
RUN_ID=${RUN_ID:-e3_sweep_${SLURM_JOB_ID:-manual}}

if [ -n "${SLURM_SUBMIT_DIR:-}" ] && [ -f "$SLURM_SUBMIT_DIR/base.sh" ]; then
  jr_scripts_dir="$SLURM_SUBMIT_DIR"
else
  jr_scripts_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi

BATCH_SIZE=${BATCH_SIZE:-48}

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

# bench_prepare_inputs() (base.sh) generates atm/obs tables via
# "${JR_BIN_DIR:-$JR_SRC_DIR}/climatology" and the geometry binary. Without
# JR_BIN_DIR being set, that falls back to $JR_SRC_DIR -- i.e. whatever
# happens to currently be sitting in the shared repo checkout's src/, which
# can be (and has been) a GPU/OpenACC build left behind by an unrelated
# run_ncu_profile*.sh job (missing CPU-job modules like nvidia-compilers,
# failing with "libacchost.so: cannot open shared object file"), or a build
# compiled for a different ND/NG than this sweep point needs (failing in
# read_ctl() with "Set 0 <= ND <= MAX!"). build_or_reuse() below builds
# climatology/geometry alongside formod, per (ND,NG), precisely to avoid
# both failure modes -- see run_point(), which sets JR_BIN_DIR per point.
# No one-time build is needed here: every point sets its own JR_BIN_DIR
# before bench_prepare_inputs runs.

CONFIG_DIR="$JR_REPO_ROOT/projects/benchmark/configs"
CHANNEL_FILE="$CONFIG_DIR/channels_alt3.tsv"
CHANNEL_LIST="${CHANNEL_LIST:-8 16 32 64 128}"
GAS_SETS_DIR="$CONFIG_DIR/gas_sets"
 
THREADS="${THREADS:-${JR_PHYS_PER_SOCKET}}"
CORES="${CORES:-$(cpus_phys "$THREADS")}"
REP="${REP:-3}"

# ND/NG are compile-time array bounds (see src/jurassic.h)
BIN_CACHE_DIR="$JR_REPO_ROOT/projects/benchmark/bin_cache"
BUILD_SCRATCH_DIR="$BIN_CACHE_DIR/_build"
mkdir -p "$BIN_CACHE_DIR" "$BUILD_SCRATCH_DIR"

# build_or_reuse ND NG
# formod's statically-sized arrays (los_t/atm_t etc., see src/jurassic.h) are
# bounded by the compile-time ND/NG macros, so climatology and the geometry
# binary (zenith/nadir/limb) need the SAME per-point -DND/-DNG as formod --
# not just formod itself. Building only formod here (as this used to) while
# leaving climatology/geometry at whatever ND/NG they were last compiled with
# means any sweep point past that bound fails in read_ctl() with
# "Set 0 <= ND <= MAX!" as soon as the .ctl's runtime ND exceeds the binary's
# compiled-in array size. All three binaries are cached together per (nd,ng)
# key below so every sweep point is internally consistent.
build_or_reuse() {
    local nd="$1" ng="$2"
    local key="nd${nd}_ng${ng}_${JR_MARCH_TAG}"
    local variant_dir="$BIN_CACHE_DIR/${key}"
    local variant_bin="$variant_dir/formod"
    local variant_climatology="$variant_dir/climatology"
    local variant_geom="$variant_dir/$JR_GEOMETRY"

    if [[ -x "$variant_bin" && -x "$variant_climatology" && -x "$variant_geom" ]]; then
        echo "[e3] reusing cached build for ND=${nd} NG=${ng}" >&2
        echo "$variant_bin"
        return
    fi

    local lock_dir="$BUILD_SCRATCH_DIR/${key}.lock"
    local waited=0
    while ! mkdir "$lock_dir" 2>/dev/null; do
        sleep 5
        waited=$((waited + 5))
        if [[ -x "$variant_bin" && -x "$variant_climatology" && -x "$variant_geom" ]]; then
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

    if [[ -x "$variant_bin" && -x "$variant_climatology" && -x "$variant_geom" ]]; then
        echo "[e3] reusing cached build for ND=${nd} NG=${ng}" >&2
        echo "$variant_bin"
        return
    fi

    local root_dir="$BUILD_SCRATCH_DIR/${key}"
    local build_dir="$root_dir/src"
    echo "[e3] building ND=${nd} NG=${ng} (formod, climatology, $JR_GEOMETRY) in private copy -> $build_dir" >&2
    rm -rf "$root_dir"
    mkdir -p "$root_dir"
    for entry in "$JR_REPO_ROOT"/*; do
        local base
        base="$(basename "$entry")"
        if [[ "$base" == "src" ]]; then
            cp -r "$entry" "$root_dir/src"
        else
            ln -s "$entry" "$root_dir/$base"
        fi
    done

    if ! ( cd "$build_dir" \
        && make clean \
        && make -j formod climatology "$JR_GEOMETRY" MARCH_NATIVE="$JR_MARCH_NATIVE" MPI="$JR_MPI" MPICC="$JR_MPICC" COMPILER="$JR_COMPILER" \
             GPU=0 LIKWID=1 DEFINES="-DND=${nd} -DNG=${ng}" ) 1>&2; then
        echo "[e3] ERROR: build failed for ND=${nd} NG=${ng} -> $build_dir" >&2
        rm -rf "$root_dir"
        exit 1
    fi

    if [[ ! -x "$build_dir/formod" || ! -x "$build_dir/climatology" || ! -x "$build_dir/$JR_GEOMETRY" ]]; then
        echo "[e3] ERROR: build for ND=${nd} NG=${ng} did not produce formod/climatology/$JR_GEOMETRY -> $build_dir" >&2
        rm -rf "$root_dir"
        exit 1
    fi

    mkdir -p "$variant_dir"
    cp "$build_dir/formod" "$variant_bin"
    cp "$build_dir/climatology" "$variant_climatology"
    cp "$build_dir/$JR_GEOMETRY" "$variant_geom"
    rm -rf "$root_dir"

    echo "$variant_bin"
}


run_point() {
    local label="$1" nd="$2" ng="$3" gas_file="$4"

    local ctl_out="$JR_WORK_DIR/ctl/${label}.ctl"
    mkdir -p "$(dirname "$ctl_out")"

    python3 "$SCRIPT_DIR/generate_ctl.py" --channels "$CHANNEL_FILE" --nd "$nd" \
        --gas-file "${gas_file:-$BASE_GAS_FILE}" "$JR_ACTIVE_CTL_BASE" "$ctl_out" >/dev/null

    JR_ACTIVE_CTL="$ctl_out"
    JR_FORMOD_BIN="$(build_or_reuse "$nd" "$ng")"
    export JR_FORMOD_BIN
    # bench_prepare_inputs() (base.sh) falls back to "${JR_BIN_DIR:-$JR_SRC_DIR}"
    # for climatology/geometry; point it at this point's own ND/NG-sized
    # build so it never reaches for the shared, uncontrolled $JR_SRC_DIR.
    JR_BIN_DIR="$(dirname "$JR_FORMOD_BIN")"
    export JR_BIN_DIR

    bench_prepare_inputs
    bench_run_forward "$label" "$THREADS" "FLOPS_DP" "$BATCH_SIZE" "$REP" "$CORES" ""
    bench_run_forward "$label" "$THREADS" "MEM_DP"   "$BATCH_SIZE" "$REP" "$CORES" ""
}

JR_ACTIVE_CTL_BASE="$JR_ACTIVE_CTL"

# NG held fixed at the selected baseline case's gas count while ND varies;
# ND held fixed at the baseline channel count while NG varies.
BASELINE_ND="${BASELINE_ND:-32}"
BASELINE_NG="${BASELINE_NG:-18}"
BASE_GAS_FILE="$GAS_SETS_DIR/ng18.txt"

echo "=== pre-building formod for all ND/NG variants ==="
for nd in $CHANNEL_LIST; do
    build_or_reuse "$nd" "$BASELINE_NG" >/dev/null
done

for gas_file in "$GAS_SETS_DIR"/*.txt; do
    [[ -e "$gas_file" ]] || continue
    ng="$(grep -c '^[^#[:space:]]' "$gas_file")"
    build_or_reuse "$BASELINE_ND" "$ng" >/dev/null
done
echo "=== pre-build complete ==="

echo "=== channel scaling ==="
for nd in $CHANNEL_LIST; do
    run_point "channels_${nd}" "$nd" "$BASELINE_NG" ""
done

echo "=== gas set scaling ==="
for gas_file in "$GAS_SETS_DIR"/*.txt; do
    [[ -e "$gas_file" ]] || continue
    name="$(basename "$gas_file" .txt)"
    ng="$(grep -c '^[^#[:space:]]' "$gas_file")"
    run_point "gases_${name}" "$BASELINE_ND" "$ng" "$gas_file"
done

cp -a "$JR_WORK_DIR/ctl" "$JR_RUN_DIR/" 2>/dev/null || true
bench_finish