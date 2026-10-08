#!/bin/bash
#SBATCH --account=slmet
#SBATCH --partition=booster
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --gpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --exclusive
#SBATCH --disable-dcgm
#SBATCH --disable-perfparanoid

# GPU counterpart of run_hermes_profile.sh: sweeps Nsight Compute (ncu) over
# BATCH_SIZE instead of sweeping LIKWID over OMP_NUM_THREADS. 


set -euo pipefail
set -x
trap 'echo "FAILED at line $LINENO: $BASH_COMMAND" >&2' ERR

# Resolve repository-relative paths once so Slurm can stage the run from temporary launch directories.
script_source=${BASH_SOURCE[0]:-$0}
script_dir=$(cd "$(dirname "$script_source")" && pwd)
repo_root=$(cd "$script_dir/../../.." && pwd)
if [ ! -f "$repo_root/projects/benchmark/configs/baseline_cases.tsv" ] && [ -n "${SLURM_SUBMIT_DIR:-}" ]; then
  submit_repo_root=$(cd "$SLURM_SUBMIT_DIR/../../.." && pwd)
  if [ -f "$submit_repo_root/projects/benchmark/configs/baseline_cases.tsv" ]; then
    repo_root=$submit_repo_root
    script_dir="$repo_root/projects/benchmark/scripts"
  fi
fi
src_dir="$repo_root/src"
runs_root="$repo_root/projects/benchmark/runs"
run_id=${RUN_ID:-ncu_profile_report_${SLURM_JOB_ID:-manual}}
run_dir="$runs_root/$run_id"
work_dir="$run_dir/work"

# Select the benchmark case from the shared baseline matrix (same lookup as run_hermes_profile.sh).
case_name=${CASE_NAME:-${GEOMETRY:-zenith}_baseline}
baseline_cases="$repo_root/projects/benchmark/configs/baseline_cases.tsv"
case_row=$(awk -F'	' -v key="$case_name" 'NR > 1 && $1 == key { print; exit }' "$baseline_cases")
if [ -z "$case_row" ]; then
  echo "Unknown benchmark baseline case: $case_name" >&2
  exit 1
fi

geometry=$(printf '%s\n' "$case_row" | awk -F'	' '{print $2}')
ctl_rel=$(printf '%s\n' "$case_row" | awk -F'	' '{print $3}')
ctl_template=${CTLFILE:-$repo_root/$ctl_rel}
bench_tblbase=${BENCH_TBLBASE:-/p/data1/slmet/model_data/jurassic/tab/tria_1cm/nc_1e-6/tria}

compiler_gpu=${COMPILER_GPU:-nvc}
mpicc=${MPICC:-mpicc}
mpi=${MPI:-0}
gpu_pin=${GPU_PIN:-1}
gpu_target=${GPU_TARGET:-gpu}
info=${INFO:-0}
flat_arrays=${FLAT_ARRAYS:-1}
rebuild=${REBUILD:-1}

ncu_batch_sizes=${NCU_BATCH_SIZES:-"8 32 128 256 1024 2048"}

# DP FLOPs = dadd + dmul + 2*dfma; bytes = DRAM read + write; occupancy,
# compute/memory throughput %, and L1/L2 cache traffic + hit rate 
ncu_metrics=${NCU_METRICS:-"gpu__time_duration.sum,smsp__sass_thread_inst_executed_op_dadd_pred_on.sum,smsp__sass_thread_inst_executed_op_dmul_pred_on.sum,smsp__sass_thread_inst_executed_op_dfma_pred_on.sum,dram__bytes_read.sum,dram__bytes_write.sum,sm__warps_active.avg.pct_of_peak_sustained_active,sm__throughput.avg.pct_of_peak_sustained_elapsed,dram__throughput.avg.pct_of_peak_sustained_elapsed,l1tex__t_bytes.sum,l1tex__t_sector_hit_rate.pct,lts__t_bytes.sum,lts__t_sector_hit_rate.pct"}

# Launch 1 of formod_batch is always the one-element reference call; launch 2
# is the first timed benchmark iteration.
ncu_launch_skip=${NCU_LAUNCH_SKIP:-1}
ncu_launch_count=${NCU_LAUNCH_COUNT:-1}
ncu_max_iter=$((ncu_launch_skip + ncu_launch_count))

mkdir -p "$work_dir"

echo "perf_event_paranoid: $(cat /proc/sys/kernel/perf_event_paranoid 2>/dev/null || echo unavailable)" > "$run_dir/perf_paranoid_status.txt"

# Validate the selected control template and LUT base path before staging work files.
if [ ! -f "$ctl_template" ]; then
  echo "Control file not found: $ctl_template" >&2
  exit 1
fi
bench_tbl_dir=$(dirname "$bench_tblbase")
if [ ! -d "$bench_tbl_dir" ]; then
  echo "Benchmark LUT directory not found: $bench_tbl_dir" >&2
  exit 1
fi

case "$geometry" in
  zenith|nadir|limb)
    ;;
  *)
    echo "Unsupported geometry: $geometry" >&2
    exit 1
    ;;
esac

cd "$work_dir"
export LANG=C
export LC_ALL=C

# Load the GPU compiler, Nsight Compute, and the plotting/validation stack
# expected on JUWELS Booster.
if command -v ml >/dev/null 2>&1; then
  ml Stages/2026 GCCcore/14.3.0
  ml CMake/4.0.3
  ml ecBuild
  ml nvidia-compilers ParaStationMPI
  ml Nsight-Compute/2025.3.1
  ml SciPy-bundle/2025.07
  ml netcdf4-python/1.7.2
fi

if ! command -v ncu >/dev/null 2>&1; then
  echo "ncu not found on PATH after 'ml Nsight-Compute'." >&2
  echo "Check module availability with: ml spider Nsight-Compute" >&2
  exit 1
fi

ncu --version | tee "$run_dir/ncu_version.txt"
ncu --query-metrics > "$run_dir/ncu_available_metrics.txt" 2>&1 || true

export LD_LIBRARY_PATH="$repo_root/libs/build/lib:$repo_root/libs/build/lib64:${LD_LIBRARY_PATH:-}"

active_ctl="$work_dir/${case_name}.ctl"
awk -v tblbase="$bench_tblbase" '{ if ($1 == "TBLBASE") print "TBLBASE = " tblbase; else print $0; }' "$ctl_template" > "$active_ctl"

collect_hardware_info() {
  if command -v lscpu >/dev/null 2>&1; then
    lscpu > "$run_dir/lscpu.txt"
  fi
  if command -v numactl >/dev/null 2>&1; then
    numactl --hardware > "$run_dir/numactl_hardware.txt"
  fi
  if command -v nvidia-smi >/dev/null 2>&1; then
    nvidia-smi --query-gpu=index,name,driver_version,compute_cap,memory.total,memory.free \
      --format=csv > "$run_dir/gpu_info.csv" 2>&1 || true
    nvidia-smi topo -m > "$run_dir/gpu_topology.txt" 2>&1 || true
  fi
}

collect_hardware_info

record_git_info() {
  local git_info="$run_dir/git_info.txt"
  if git -C "$repo_root" rev-parse --is-inside-work-tree >/dev/null 2>&1; then
    {
      echo "commit=$(git -C "$repo_root" rev-parse HEAD)"
      echo "branch=$(git -C "$repo_root" rev-parse --abbrev-ref HEAD)"
      if [ -z "$(git -C "$repo_root" status --porcelain)" ]; then
        echo "dirty=0"
      else
        echo "dirty=1"
      fi
    } > "$git_info"
    git -C "$repo_root" status --porcelain > "$run_dir/git_status.txt" 2>/dev/null || true
  else
    echo "repo_root is not a git work tree: $repo_root" > "$git_info"
  fi
}

record_git_info

printf 'case_name=%s\ngeometry=%s\nctl_template=%s\nactive_ctl=%s\nbench_tblbase=%s\ncompiler_gpu=%s\nmpicc=%s\nmpi=%s\ngpu_pin=%s\ngpu_target=%s\nflat_arrays=%s\nrebuild=%s\nncu_batch_sizes=%s\nncu_metrics=%s\nncu_launch_skip=%s\nncu_launch_count=%s\n' \
  "$case_name" \
  "$geometry" \
  "$ctl_template" \
  "$active_ctl" \
  "$bench_tblbase" \
  "$compiler_gpu" \
  "$mpicc" \
  "$mpi" \
  "$gpu_pin" \
  "$gpu_target" \
  "$flat_arrays" \
  "$rebuild" \
  "$ncu_batch_sizes" \
  "$ncu_metrics" \
  "$ncu_launch_skip" \
  "$ncu_launch_count" \
  > "$run_dir/config.txt"

build_gpu() {
  cd "$src_dir" || return 1
  make clean || return 1
  make -j DEFINES="-DNG=18" MPI="$mpi" MPICC="$mpicc" COMPILER="$compiler_gpu" GPU=1 GPU_TARGET="$gpu_target" \
    GPU_PIN="$gpu_pin" INFO="$info" FLAT_ARRAYS="$flat_arrays" LIKWID=0 || return 1
  # Return to work_dir (may have been entered via Slurm's temporary launch dir)
  cd "$work_dir" 2>/dev/null || true
  return 0
}

# Create atmospheric and observation inputs for the selected geometry.
prepare_inputs() {
  rm -rf data
  mkdir -p data
  "$src_dir/climatology" "$active_ctl" data/atm.tab
  "$src_dir/$geometry" "$active_ctl" data/obs.tab
}

if [ "$rebuild" = 1 ]; then
  build_gpu
fi

validation_status="$run_dir/validation_status.txt"
run_profiling=1

if [ "${SKIP_VALIDATION:-0}" != "1" ]; then
  trap - ERR
  set +e
  ( cd "$repo_root/projects/validation" && \
    VALIDATION_TBLBASE="$bench_tblbase" scripts/run_validation.py \
      > "$run_dir/validation.log" 2>&1 )
  validation_rc=$?
  set -e
  trap 'echo "FAILED at line $LINENO: $BASH_COMMAND" >&2' ERR

  echo "exit_code=$validation_rc" > "$validation_status"

  latest_validation_run=$(ls -td "$repo_root/projects/validation/runs/validation_"*/ 2>/dev/null | head -n1 || true)
  if [ -n "$latest_validation_run" ]; then
    cp -a "${latest_validation_run}summary.tsv" "$run_dir/validation_summary.tsv" 2>/dev/null || true
  fi
  if [ "$validation_rc" -ne 0 ]; then
    echo "Validation FAILED (exit $validation_rc) -- skipping Nsight Compute sweep for this candidate." >&2
    run_profiling=0
  fi
else
  echo "exit_code=skipped" > "$validation_status"
fi

mkdir -p ncu
mkdir -p "$run_dir/ncu"
cd ncu
prepare_inputs

# Validate requested metrics against what this GPU/ncu build actually supports,
# mirroring the LIKWID-group check in run_hermes_profile.sh.
if [ "$run_profiling" = 1 ]; then
  IFS=',' read -r -a requested_metrics <<< "$ncu_metrics"
  missing_metrics=()
  for metric in "${requested_metrics[@]}"; do
    base_metric=${metric%%.*}
    if ! grep -q -w "$base_metric" "$run_dir/ncu_available_metrics.txt" 2>/dev/null; then
      missing_metrics+=("$metric")
    fi
  done
  if [ ${#missing_metrics[@]} -gt 0 ]; then
    echo "WARNING: metric(s) not found in ncu --query-metrics output, continuing anyway: ${missing_metrics[*]}" >&2
    echo "See $run_dir/ncu_available_metrics.txt for the supported list on this node." >&2
  fi
fi

if [ "$run_profiling" = 1 ]; then
  for batch in $ncu_batch_sizes; do
    ncu_output="$run_dir/ncu/formod_batch${batch}"
    log_txt="log.batch${batch}.txt"
    out_tab="/tmp/jurassic_ncu_${run_id}_b${batch}.tab"

    echo "Running Nsight Compute BATCH_SIZE=$batch ..."

    trap - ERR
    set +e
    JURASSIC_MAX_ITER=$ncu_max_iter srun -n1 -N1 ncu \
      --target-processes all \
      --kernel-name-base function \
      -k regex:formod_batch \
      --launch-skip "$ncu_launch_skip" --launch-count "$ncu_launch_count" \
      --metrics "$ncu_metrics" \
      --force-overwrite \
      --export "$ncu_output" \
      "$src_dir/formod" "$active_ctl" data/obs.tab data/atm.tab "$out_tab" \
      TASK time BATCH_SIZE "$batch" \
      > "$log_txt" 2>&1
    ncu_rc=$?
    set -e
    trap 'echo "FAILED at line $LINENO: $BASH_COMMAND" >&2' ERR

    printf 'BATCH_SIZE=%s\nLAUNCH_SKIP=%s\nLAUNCH_COUNT=%s\nEXIT_CODE=%s\n' \
      "$batch" "$ncu_launch_skip" "$ncu_launch_count" "$ncu_rc" >> "$log_txt"

    if [ "$ncu_rc" -ne 0 ]; then
      echo "WARNING: ncu run failed for BATCH_SIZE=$batch (exit $ncu_rc), see $log_txt" >&2
    fi

    csv_out="log.batch${batch}.csv"
    ncu_rep="${ncu_output}.ncu-rep"
    if [ -f "$ncu_rep" ]; then
      trap - ERR
      set +e
      ncu --import "$ncu_rep" --csv --page raw > "$csv_out" 2>"log.batch${batch}.decode_err.txt"
      decode_rc=$?
      set -e
      trap 'echo "FAILED at line $LINENO: $BASH_COMMAND" >&2' ERR
      if [ "$decode_rc" -ne 0 ]; then
        echo "WARNING: 'ncu --import' failed for BATCH_SIZE=$batch (exit $decode_rc)," \
             "see log.batch${batch}.decode_err.txt -- $csv_out will be missing/incomplete." >&2
      elif [ ! -s "$csv_out" ]; then
        echo "WARNING: 'ncu --import' produced an empty $csv_out for BATCH_SIZE=$batch." >&2
      fi
    else
      echo "WARNING: expected $ncu_rep not found, skipping CSV decode for BATCH_SIZE=$batch." >&2
    fi

    rm -f "$out_tab"
  done
else
  echo "Skipped Nsight Compute sweep -- validation failed for this candidate, or ncu is unavailable on this node." > "$run_dir/skipped_profiling.txt"
fi

cp -a data "$run_dir/data.ncu"
cp -a log.batch*.txt "$run_dir/" 2>/dev/null || true
cp -a log.batch*.csv "$run_dir/" 2>/dev/null || true
cp -a "$active_ctl" "$run_dir/"

echo "Nsight Compute run directory: $run_dir"
echo "Raw per-config output: $run_dir/log.batch<N>.txt, $run_dir/log.batch<N>.csv (decoded metrics), and $run_dir/ncu/formod_batch<N>.ncu-rep"
echo "Available metrics on this node: $run_dir/ncu_available_metrics.txt"
echo "GPU info: $run_dir/gpu_info.csv, $run_dir/gpu_topology.txt"
echo "Code provenance: $run_dir/git_info.txt"

# Fail the job itself when the candidate didn't validate
if [ "${validation_rc:-0}" -ne 0 ]; then
  echo "Exiting non-zero: validation failed for this candidate (see $validation_status)." >&2
  exit 3
fi
