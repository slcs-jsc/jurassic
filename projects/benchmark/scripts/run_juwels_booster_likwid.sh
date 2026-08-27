#!/bin/bash
#SBATCH --account=slmet
#SBATCH --partition=booster
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --time=00:45:00
#SBATCH --exclusive

# LIKWID reads hardware performance counters (cache misses, memory bandwith, FLOP rate etc.)

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
run_id=${RUN_ID:-juwels_booster_likwid_${SLURM_JOB_ID:-manual}}
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
cpu_batch_size=${CPU_BATCH_SIZE:-48}
compiler_cpu=${COMPILER_CPU:-gcc}
mpicc=${MPICC:-mpicc}
mpi=${MPI:-0}
rebuild=${REBUILD:-1}

# LIKWID Setup
likwid_threads=${LIKWID_THREADS:-"1 12 48"}
likwid_groups=${LIKWID_GROUPS:-"MEM_DP FLOPS_DP CACHES"} # performance groups to collect

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

cd "$work_dir"
export LANG=C
export LC_ALL=C

# Load the compiler and plotting stack expected on JUWELS Booster.
if command -v ml >/dev/null 2>&1; then
  ml Stages/2026 GCC/14.3.0
  module load likwid/5.4.1
  ml CMake/4.0.3
  ml ecBuild
fi

if ! command -v likwid-perfctr >/dev/null 2>&1; then
  echo "likwid-perfctr not found on PATH after 'ml likwid'." >&2
  echo "Check module availability with: ml spider likwid" >&2
  exit 1
fi

likwid-perfctr -a > "$run_dir/likwid_available_groups.txt" 2>&1 || true

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

printf 'case_name=%s\ngeometry=%s\nctl_template=%s\nactive_ctl=%s\nbench_tblbase=%s\ncpu_batch_size=%s\ncompiler_cpu=%s\nmpicc=%s\nmpi=%s\nrebuild=%s\nlikwid_threads=%s\nlikwid_groups=%s\n' \
  "$case_name" "$geometry" "$ctl_template" "$active_ctl" "$bench_tblbase" "$cpu_batch_size" "$compiler_cpu" "$mpicc" "$mpi" "$rebuild" "$likwid_threads" "$likwid_groups" > "$run_dir/config.txt"


# Rebuild a CPU-only binary when the run requests it.
build_cpu() {
  cd "$src_dir"
  make clean
  make -j MPI="$mpi" MPICC="$mpicc" COMPILER="$compiler_cpu" GPU=0
  cd "$work_dir"
}

# Create atmospheric and observation inputs for the selected geometry.
prepare_inputs() {
  rm -rf data
  mkdir -p data
  "$src_dir/climatology" "$active_ctl" data/atm.tab
  "$src_dir/$geometry" "$active_ctl" data/obs.tab
}

if [ "$rebuild" = 1 ]; then
  build_cpu
fi

# Perform Validation before Profiling
validation_status="$run_dir/validation_status.txt"
run_profiling=1
if [ "${SKIP_VALIDATION:-0}" != "1" ]; then
  ( cd "$repo_root/projects/validation" && \
    VALIDATION_TBLBASE="$bench_tblbase" scripts/run_validation.py \
    > "$run_dir/validation.log" 2>&1 )
  validation_rc=$?
  echo "exit_code=$validation_rc" > "$validation_status"
  
  latest_validation_run=$(ls -td "$repo_root/projects/validation"/runs/validation_*/ 2>/dev/null | head -n1)
  if [ -n "$latest_validation_run" ]; then
    cp -a "${latest_validation_run}summary.tsv" "$run_dir/validation_summary.tsv" 2>/dev/null || true
  fi
  if [ "$validation_rc" -ne 0 ]; then
    echo "Validation FAILED (exit $validation_rc) -- skipping LIKWID profiling for this candidate." >&2
    run_profiling=0
  fi
else
  echo "exit_code=skipped" > "$validation_status"
fi

mkdir -p likwid
cd likwid
prepare_inputs

if [ "$run_profiling" = 1 ]; then
  # Validate LIKWID groups against what this node actually supports
  valid_groups=()
  for group in $likwid_groups; do
    if grep -q -w "$group" "$run_dir/likwid_available_groups.txt"; then
      valid_groups+=("$group")
    else
      echo "WARNING: LIKWID group '$group' not available on this node -- skipping. See likwid_available_groups.txt for valid options." >&2
    fi
  done

  if [ ${#valid_groups[@]} -eq 0 ]; then
    echo "No requested LIKWID groups are available on this node -- skipping LIKWID sweep." >&2
    echo "Check $run_dir/likwid_available_groups.txt for valid group names and rerun with LIKWID_GROUPS set accordingly." >&2
    run_profiling=0
  fi
fi

if [ "$run_profiling" = 1 ]; then
  # unset OMP_PLACES OMP_PROC_BIND ?
  for omp in $likwid_threads; do
    core_list="S0:0-$((omp - 1))"
    out_tab="/tmp/jurassic_bench_${run_id}_likwid_omp${omp}.tab"
  
    for group in "${valid_groups[@]}"; do
      log_txt="log.omp${omp}.${group}.txt"
      log_csv="log.omp${omp}.${group}.csv"
  
      echo "Running LIKWID group=$group OMP_NUM_THREADS=$omp ..."
  
      OMP_NUM_THREADS=$omp likwid-perfctr -C "$core_list" -g "$group" \
        -o "$log_csv" \
        "$src_dir/formod" "$active_ctl" data/obs.tab data/atm.tab "$out_tab" \
        TASK time BATCH_SIZE "$cpu_batch_size" \
        > "$log_txt" 2>&1
  
      printf 'OMP_NUM_THREADS=%s\nLIKWID_GROUP=%s\nCPU_BATCH_SIZE=%s\nCORE_LIST=%s\n' \
        "$omp" "$group" "$cpu_batch_size" "$core_list" >> "$log_txt"
    done
  done
else 
  echo "Skipped LIKWID sweep -- validation failed for this candidate, or no requested LIKWID groups were available on this node." > skipped_profiling.txt
fi

cp -a data "$run_dir/data.likwid"
cp -a log.omp*.txt log.omp*.csv "$run_dir/" 2>/dev/null || true
cp -a "$active_ctl" "$run_dir/"
 
echo "LIKWID run directory: $run_dir"
echo "Raw per-config output: $run_dir/log.omp<N>.<GROUP>.txt (human-readable) and .csv (machine-readable)"
echo "Available groups on this node: $run_dir/likwid_available_groups.txt"