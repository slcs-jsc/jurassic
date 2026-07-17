#!/bin/bash
#SBATCH --account=slmet
#SBATCH --partition=booster
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --gpus-per-task=1
# The default Slurm request matches one quarter of a Booster node: one A100
# plus 12 CPU cores. Use sbatch --cpus-per-task=48 for CPU-only node tests.
#SBATCH --time=01:30:00

set -euo pipefail

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
run_id=${RUN_ID:-juwels_booster_${SLURM_JOB_ID:-manual}}
run_dir="$runs_root/$run_id"
work_dir="$run_dir/work"

# Select the benchmark case from the shared baseline matrix.
case_name=${CASE_NAME:-${GEOMETRY:-zenith}_baseline}
baseline_cases="$repo_root/projects/benchmark/configs/baseline_cases.tsv"
case_row=$(awk -F'	' -v key="$case_name" 'NR > 1 && $1 == key { print; exit }' "$baseline_cases")
if [ -z "$case_row" ]; then
  echo "Unknown benchmark baseline case: $case_name" >&2
  exit 1
fi
geometry=$(printf '%s
' "$case_row" | awk -F'	' '{print $2}')
ctl_rel=$(printf '%s
' "$case_row" | awk -F'	' '{print $3}')
ctl_template=${CTLFILE:-$repo_root/$ctl_rel}
bench_tblbase=${BENCH_TBLBASE:-/p/data1/slmet/model_data/jurassic/tab/tria_1cm/nc_1e-6/tria}
target=${TARGET:-gpu}
slurm_cpus_per_task=${SLURM_CPUS_PER_TASK:-12}
full_node_cpu_cores=${FULL_NODE_CPU_CORES:-48}
if [ -n "${THREADS:-}" ]; then
  threads=$THREADS
elif [ "$target" = gpu ]; then
  threads="1 2 4 8 12"
else
  threads="1 2 4 8 12 24 48"
fi
cpu_batch_size=${CPU_BATCH_SIZE:-48}
batches=${BATCHES:-"1 8 64 256"}
compiler_cpu=${COMPILER_CPU:-gcc}
compiler_gpu=${COMPILER_GPU:-nvc}
mpicc=${MPICC:-mpicc}
mpi=${MPI:-1}
gpu_pin=${GPU_PIN:-1}
info=${INFO:-0}
acc_time=${ACC_TIME:-0}
acc_notify=${ACC_NOTIFY:-0}
rebuild=${REBUILD:-1}
omp_places=${OMP_PLACES:-cores}
omp_proc_bind=${OMP_PROC_BIND:-close}

mkdir -p "$work_dir"

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

case "$target" in
  cpu|gpu|both)
    ;;
  *)
    echo "Unsupported target: $target" >&2
    exit 1
    ;;
esac

if [ "$target" = both ] && [ "$rebuild" = 0 ]; then
  echo "TARGET=both requires REBUILD=1 so CPU and GPU binaries are rebuilt separately" >&2
  exit 1
fi

if [ "$target" = cpu ] || [ "$target" = both ]; then
  if [ -z "${THREADS:-}" ] && [ "$slurm_cpus_per_task" -lt "$full_node_cpu_cores" ]; then
    echo "CPU benchmark target needs --cpus-per-task=$full_node_cpu_cores for the default full-node thread sweep" >&2
    echo "Override THREADS explicitly if you want a smaller CPU sweep inside a GPU-sized allocation" >&2
    exit 1
  fi
fi

# Run in a clean working directory with deterministic locale settings.
cd "$work_dir"
export LANG=C
export LC_ALL=C
export OMP_PLACES="$omp_places"
export OMP_PROC_BIND="$omp_proc_bind"
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-$slurm_cpus_per_task}

# Load the compiler and plotting stack expected on JUWELS Booster.
if command -v ml >/dev/null 2>&1; then
  ml Stages/2026 GCCcore/14.3.0
  ml CMake/4.0.3
  ml ecBuild
  ml matplotlib/3.10.5
  ml nvidia-compilers ParaStationMPI
fi

export LD_LIBRARY_PATH="$repo_root/libs/build/lib:$repo_root/libs/build/lib64:${LD_LIBRARY_PATH:-}"

# Materialize a run-local control file with the chosen LUT base name.
active_ctl="$work_dir/${case_name}.ctl"
awk -v tblbase="$bench_tblbase" '{ if ($1 == "TBLBASE") print "TBLBASE = " tblbase; else print $0; }' "$ctl_template" > "$active_ctl"

# Record the effective benchmark configuration for later inspection.
collect_hardware_info() {
  if command -v lscpu >/dev/null 2>&1; then
    lscpu > "$run_dir/lscpu.txt"
  fi
  if command -v numactl >/dev/null 2>&1; then
    numactl --hardware > "$run_dir/numactl_hardware.txt"
  fi
}

collect_hardware_info

printf 'case_name=%s
geometry=%s
ctl_template=%s
active_ctl=%s
bench_tblbase=%s
target=%s
slurm_cpus_per_task=%s
full_node_cpu_cores=%s
threads=%s
cpu_batch_size=%s
batches=%s
compiler_cpu=%s
compiler_gpu=%s
mpicc=%s
mpi=%s
gpu_pin=%s
info=%s
acc_time=%s
acc_notify=%s
rebuild=%s
omp_places=%s
omp_proc_bind=%s
'   "$case_name" "$geometry" "$ctl_template" "$active_ctl" "$bench_tblbase" "$target" "$slurm_cpus_per_task" "$full_node_cpu_cores" "$threads" "$cpu_batch_size" "$batches" "$compiler_cpu" "$compiler_gpu" "$mpicc" "$mpi" "$gpu_pin" "$info" "$acc_time" "$acc_notify" "$rebuild" "$omp_places" "$omp_proc_bind" > "$run_dir/config.txt"

# Rebuild a CPU-only binary when the run requests it.
build_cpu() {
  cd "$src_dir"
  make clean
  make -j MPI="$mpi" MPICC="$mpicc" COMPILER="$compiler_cpu" GPU=0
  cd "$work_dir"
}

# Rebuild a GPU-enabled binary with the requested compiler settings.
build_gpu() {
  cd "$src_dir"
  make clean
  make -j MPI="$mpi" MPICC="$mpicc" COMPILER="$compiler_gpu" GPU=1 GPU_PIN="$gpu_pin" INFO="$info"
  cd "$work_dir"
}

# Keep summary generation working even on systems without matplotlib.
maybe_plot() {
  local summary_tsv=$1
  local output_png=$2
  local title=$3

  if python3 -c "import matplotlib" >/dev/null 2>&1; then
    python3 "$script_dir/plot_benchmark_results.py" "$summary_tsv" --output "$output_png" --title "$title"
  else
    echo "WARNING: matplotlib not available; skipping plot generation for $output_png" >&2
  fi
}

# Create atmospheric and observation inputs for the selected geometry.
prepare_inputs() {
  rm -rf data
  mkdir -p data
  "$src_dir/climatology" "$active_ctl" data/atm.tab
  "$src_dir/$geometry" "$active_ctl" data/obs.tab
}

# Execute the CPU benchmark sweep over OpenMP thread counts.
run_cpu() {
  mkdir -p cpu
  cd cpu
  prepare_inputs
  # Sweep OpenMP thread counts while keeping the synthetic workload fixed.
  for omp in $threads; do
    log="log.omp${omp}"
    out="/tmp/jurassic_bench_${run_id}_cpu_omp${omp}.tab"
    OMP_NUM_THREADS=$omp "$src_dir/formod" "$active_ctl" data/obs.tab data/atm.tab "$out" TASK time BATCH_SIZE "$cpu_batch_size" > "$log" 2>&1
    printf 'OMP_NUM_THREADS=%s
CPU_BATCH_SIZE=%s
OMP_PLACES=%s
OMP_PROC_BIND=%s
' "$omp" "$cpu_batch_size" "$OMP_PLACES" "$OMP_PROC_BIND" >> "$log"
  done
  python3 "$script_dir/summarize_time_logs.py" omp 'log.omp*' --tsv-out "$run_dir/summary.cpu.tsv" | tee "$run_dir/summary.cpu.txt"
  maybe_plot "$run_dir/summary.cpu.tsv" "$run_dir/plot_cpu_scaling.png" "JURASSIC Booster CPU baseline"
  cp -a data "$run_dir/data.cpu"
  cp -a log.omp* "$run_dir/"
  cd "$work_dir"
}

# Execute the GPU benchmark sweep over configured batch sizes.
run_gpu() {
  mkdir -p gpu
  cd gpu
  prepare_inputs
  for batch in $batches; do
    log="log.batch${batch}"
    out="/tmp/jurassic_bench_${run_id}_gpu_batch${batch}.tab"
    if [ "$acc_time" = 1 ] || [ "$acc_notify" != 0 ]; then
      env_cmd=(env)
      if [ "$acc_time" = 1 ]; then
        env_cmd+=(NVCOMPILER_ACC_TIME=1)
      fi
      if [ "$acc_notify" != 0 ]; then
        env_cmd+=("NVCOMPILER_ACC_NOTIFY=$acc_notify")
      fi
      "${env_cmd[@]}" "$src_dir/formod" "$active_ctl" data/obs.tab data/atm.tab "$out" TASK time BATCH_SIZE "$batch" > "$log" 2>&1
    else
      "$src_dir/formod" "$active_ctl" data/obs.tab data/atm.tab "$out" TASK time BATCH_SIZE "$batch" > "$log" 2>&1
    fi
    printf 'BATCH_SIZE=%s
ACC_TIME=%s
ACC_NOTIFY=%s
' "$batch" "$acc_time" "$acc_notify" >> "$log"
  done
  python3 "$script_dir/summarize_time_logs.py" batch 'log.batch*' --tsv-out "$run_dir/summary.gpu.tsv" | tee "$run_dir/summary.gpu.txt"
  maybe_plot "$run_dir/summary.gpu.tsv" "$run_dir/plot_gpu_batch.png" "JURASSIC Booster GPU baseline"
  cp -a data "$run_dir/data.gpu"
  cp -a log.batch* "$run_dir/"
  cd "$work_dir"
}

# Run the requested CPU and/or GPU benchmark sections.
if [ "$target" = cpu ] || [ "$target" = both ]; then
  if [ "$rebuild" = 1 ]; then
    build_cpu
  fi
  run_cpu
fi

if [ "$target" = gpu ] || [ "$target" = both ]; then
  if [ "$rebuild" = 1 ]; then
    build_gpu
  fi
  run_gpu
fi

cp -a "$active_ctl" "$run_dir/"

echo "Run directory: $run_dir"
