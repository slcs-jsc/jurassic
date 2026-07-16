#!/bin/bash
set -euo pipefail

# Compare CPU and GPU formod outputs for one benchmark baseline case and flag
# brightness-temperature NaNs in the generated logs.
script_source=${BASH_SOURCE[0]:-$0}
script_dir=$(cd "$(dirname "$script_source")" && pwd)
repo_root=$(cd "$script_dir/../../.." && pwd)
src_dir="$repo_root/src"
runs_root="$repo_root/projects/benchmark/runs"
run_id=${RUN_ID:-bt_nan_$(date +%Y%m%d_%H%M%S)}
run_dir="$runs_root/$run_id"
work_dir="$run_dir/work"

# Select the benchmark baseline case and resolve the matching control template.
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
mpi_cpu=${MPI_CPU:-0}
mpi_gpu=${MPI_GPU:-1}
mpicc_gpu=${MPICC_GPU:-mpicc}
compiler_cpu=${COMPILER_CPU:-gcc}
compiler_gpu=${COMPILER_GPU:-nvc}
gpu_pin=${GPU_PIN:-1}
gpu_batch_size=${GPU_BATCH_SIZE:-1}
acc_time=${ACC_TIME:-0}
acc_notify=${ACC_NOTIFY:-0}
load_modules=${LOAD_MODULES:-1}

mkdir -p "$work_dir"

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

# Load the standard JUWELS Booster toolchain when requested.
if [ "$load_modules" = 1 ] && command -v ml >/dev/null 2>&1; then
  ml Stages/2026 GCCcore/14.3.0
  ml CMake/4.0.3
  ml ecBuild
  ml nvidia-compilers ParaStationMPI
fi

export LANG=C
export LC_ALL=C
export LD_LIBRARY_PATH="$repo_root/libs/build/lib:$repo_root/libs/build/lib64:${LD_LIBRARY_PATH:-}"

cd "$work_dir"

# Materialize a run-local control file with the chosen LUT base name.
active_ctl="$work_dir/${case_name}.ctl"
awk -v tblbase="$bench_tblbase" '{ if ($1 == "TBLBASE") print "TBLBASE = " tblbase; else print $0; }' "$ctl_template" > "$active_ctl"

# Build fresh atmospheric and observation inputs once for both runs.
prepare_inputs() {
  rm -rf data
  mkdir -p data
  "$src_dir/climatology" "$active_ctl" data/atm.tab
  "$src_dir/$geometry" "$active_ctl" data/obs.tab
}

# Build the CPU binary used for the scalar reference run.
build_cpu() {
  cd "$src_dir"
  make clean
  make -j MPI="$mpi_cpu" COMPILER="$compiler_cpu" GPU=0
  cd "$work_dir"
}

# Build the GPU binary used for the OpenACC comparison run.
build_gpu() {
  cd "$src_dir"
  make clean
  make -j MPI="$mpi_gpu" MPICC="$mpicc_gpu" COMPILER="$compiler_gpu" GPU=1 GPU_PIN="$gpu_pin"
  cd "$work_dir"
}

# Extract only the brightness-temperature lines with NaNs from a log file.
extract_bt_nan_lines() {
  local input_log=$1
  local output_txt=$2
  if grep -F 'Brightness temperature' "$input_log" | grep -F -- '-nan' > "$output_txt"; then
    :
  else
    : > "$output_txt"
  fi
}

# Create compact summary files for later inspection.
write_summary() {
  local diff_status=$1
  {
    printf 'case_name=%s
' "$case_name"
    printf 'geometry=%s
' "$geometry"
    printf 'bench_tblbase=%s
' "$bench_tblbase"
    printf 'compiler_cpu=%s
' "$compiler_cpu"
    printf 'compiler_gpu=%s
' "$compiler_gpu"
    printf 'gpu_batch_size=%s
' "$gpu_batch_size"
    printf 'acc_time=%s
' "$acc_time"
    printf 'acc_notify=%s
' "$acc_notify"
    printf 'diff_status=%s
' "$diff_status"
    printf 'cpu_bt_nan_count=%s
' "$(wc -l < bt_nan_cpu.txt)"
    printf 'gpu_bt_nan_count=%s
' "$(wc -l < bt_nan_gpu.txt)"
  } > "$run_dir/summary.txt"
}

printf 'case_name=%s
geometry=%s
ctl_template=%s
active_ctl=%s
bench_tblbase=%s
mpi_cpu=%s
mpi_gpu=%s
mpicc_gpu=%s
compiler_cpu=%s
compiler_gpu=%s
gpu_pin=%s
gpu_batch_size=%s
acc_time=%s
acc_notify=%s
load_modules=%s
'   "$case_name" "$geometry" "$ctl_template" "$active_ctl" "$bench_tblbase" "$mpi_cpu" "$mpi_gpu" "$mpicc_gpu" "$compiler_cpu" "$compiler_gpu" "$gpu_pin" "$gpu_batch_size" "$acc_time" "$acc_notify" "$load_modules" > "$run_dir/config.txt"

build_cpu
prepare_inputs
OMP_NUM_THREADS=1 "$src_dir/formod" "$active_ctl" data/obs.tab data/atm.tab "$work_dir/rad_cpu.tab" > "$work_dir/log.cpu" 2>&1
extract_bt_nan_lines "$work_dir/log.cpu" "$work_dir/bt_nan_cpu.txt"

build_gpu
if [ "$acc_time" = 1 ] || [ "$acc_notify" != 0 ]; then
  env_cmd=(env)
  if [ "$acc_time" = 1 ]; then
    env_cmd+=(NVCOMPILER_ACC_TIME=1)
  fi
  if [ "$acc_notify" != 0 ]; then
    env_cmd+=("NVCOMPILER_ACC_NOTIFY=$acc_notify")
  fi
  "${env_cmd[@]}" "$src_dir/formod" "$active_ctl" data/obs.tab data/atm.tab "$work_dir/rad_gpu.tab" TASK time BATCH_SIZE "$gpu_batch_size" > "$work_dir/log.gpu" 2>&1
else
  "$src_dir/formod" "$active_ctl" data/obs.tab data/atm.tab "$work_dir/rad_gpu.tab" TASK time BATCH_SIZE "$gpu_batch_size" > "$work_dir/log.gpu" 2>&1
fi
extract_bt_nan_lines "$work_dir/log.gpu" "$work_dir/bt_nan_gpu.txt"

if diff -q "$work_dir/rad_cpu.tab" "$work_dir/rad_gpu.tab" > "$work_dir/diff.txt" 2>&1; then
  diff_status=IDENTICAL
else
  diff_status=DIFFERENT
fi
write_summary "$diff_status"

cp "$work_dir/log.cpu" "$run_dir/"
cp "$work_dir/log.gpu" "$run_dir/"
cp "$work_dir/rad_cpu.tab" "$run_dir/"
cp "$work_dir/rad_gpu.tab" "$run_dir/"
cp "$work_dir/bt_nan_cpu.txt" "$run_dir/"
cp "$work_dir/bt_nan_gpu.txt" "$run_dir/"
cp "$work_dir/diff.txt" "$run_dir/"
cp "$active_ctl" "$run_dir/"
cp -a "$work_dir/data" "$run_dir/"

printf 'Run directory: %s
' "$run_dir"
printf 'Diff status: %s
' "$diff_status"
printf 'CPU BT NaN lines: %s
' "$(wc -l < "$work_dir/bt_nan_cpu.txt")"
printf 'GPU BT NaN lines: %s
' "$(wc -l < "$work_dir/bt_nan_gpu.txt")"
