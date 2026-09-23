#!/bin/bash
#SBATCH --account=slmet
#SBATCH --partition=booster
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --time=00:20:00
#SBATCH --exclusive
#SBATCH --disable-perfparanoid

set -euo pipefail
set -x                               
trap 'echo "FAILED at line $LINENO: $BASH_COMMAND" >&2' ERR 

script_source=${BASH_SOURCE[0]:-$0}
script_dir=$(cd "$(dirname "$script_source")" && pwd)
repo_root="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)}"

src_dir="$repo_root/src"
runs_root="$repo_root/runs"
run_id=${RUN_ID:-juwels_booster_likwid_${SLURM_JOB_ID:-manual}}
run_dir="$runs_root/$run_id"
work_dir="$run_dir/work"

max_iter=${MAX_ITER:-20}
fun_name=${FUN_NAME:-drb059}
rebuild=${REBUILD:-1}

# LIKWID & Perf Setup

if command -v ml >/dev/null 2>&1; then
  ml Stages/2026 GCC/14.3.0
  ml likwid/5.4.1
fi

likwid_threads=${LIKWID_THREADS:-"12 24 48"}
likwid_groups=${LIKWID_GROUPS:-"MEM_DP FLOPS_DP CACHE"} 

mkdir -p "$work_dir" 
mkdir -p "$run_dir"

echo "perf_event_paranoid: $(cat /proc/sys/kernel/perf_event_paranoid 2>/dev/null || echo unavailable)" > "$run_dir/perf_paranoid_status.txt"

cd "$work_dir"
export LANG=C
export LC_ALL=C

if ! command -v likwid-perfctr >/dev/null 2>&1; then
  echo "likwid-perfctr nicht gefunden!" >&2
  exit 1
fi

if ! command -v perf >/dev/null 2>&1; then
  perf_available=0
else
  perf_available=1
  perf --version > "$run_dir/perf_version.txt" 2>&1 || true
fi

likwid-perfctr -a > "$run_dir/likwid_available_groups.txt" 2>&1 || true

perf_topdown_available=0
if [ "$perf_available" -eq 1 ]; then
  if perf stat --help 2>&1 | grep -q -- '--topdown'; then
    perf_topdown_available=1
  fi
fi

collect_hardware_info() {
  if command -v lscpu >/dev/null 2>&1; then lscpu > "$run_dir/lscpu.txt"; fi
  if command -v numactl >/dev/null 2>&1; then numactl --hardware > "$run_dir/numactl_hardware.txt"; fi
}
collect_hardware_info

printf 'run_id=%s\nfun_name=%s\nmax_iter=%s\nlikwid_threads=%s\nlikwid_groups=%s\n' \
  "$run_id" "$fun_name" "$max_iter" "$likwid_threads" "$likwid_groups" > "$run_dir/config.txt"


build_cpu() {
  cd "$src_dir" || return 1
  make clean || return 1
  make -j
  cd "$work_dir" 2>/dev/null || true
  return 0
}

if [ "$rebuild" = 1 ]; then
  build_cpu
fi

valid_groups=()
for group in $likwid_groups; do
  if grep -q -w "$group" "$run_dir/likwid_available_groups.txt"; then
    valid_groups+=("$group")
  fi
done

mkdir -p "$work_dir/likwid"
cd "$work_dir/likwid"

export LIKWID_FILEPATH="$work_dir/likwid/likwid_marker_${SLURM_JOB_ID:-manual}"

if [ ${#valid_groups[@]} -gt 0 ]; then
  unset OMP_PLACES OMP_PROC_BIND
  for omp in $likwid_threads; do
    core_list="S0:0-$((omp - 1))"
  
    for group in "${valid_groups[@]}"; do
      log_txt="log.omp${omp}.${group}.txt"
      log_csv="log.omp${omp}.${group}.csv"
  
      echo "Running LIKWID group=$group OMP_NUM_THREADS=$omp ..."
  
      cp "$src_dir/benchmark" ./benchmark_exec

      OMP_NUM_THREADS=$omp likwid-perfctr -C "$core_list" -g "$group" -m \
        -o "$log_csv" \
        ./benchmark_exec "$max_iter" "$fun_name" \
        > "$log_txt" 2>&1
        
      rm -f ./benchmark_exec
    done
  done
  cp -a log.omp*.txt log.omp*.csv "$run_dir/" 2>/dev/null || true
fi

echo "Done. Results written to: $run_dir"
