#!/bin/bash
#SBATCH --account=slmet
#SBATCH --partition=dc-cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=02:00:00
#SBATCH --exclusive
#SBATCH --disable-perfparanoid
#SBATCH --job-name=e4_cache

# Cache behaviour benchmark for JURASSIC (CPU-version), LIKWID group CACHES.
# Target: JURECA DC -> AMD EPYC 7742 (2x 64 cores, 2.25 GHz)
#
# Analysis Axes:
#   1. Working-set sweep
#   2. Thread-count sweep
#
# Output: out/<category>.t<threads>.CACHES.b<batch>.rep<rep>.csv

set -euo pipefail

JR_EXPERIMENT=e4_cache
RUN_ID=${RUN_ID:-e4_cache_${SLURM_JOB_ID:-manual}}

if [ -n "${SLURM_SUBMIT_DIR:-}" ] && [ -f "$SLURM_SUBMIT_DIR/base.sh" ]; then
  jr_scripts_dir="$SLURM_SUBMIT_DIR"
else
  jr_scripts_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi

export JR_SCRIPTS_DIR_OVERRIDE="$jr_scripts_dir"
source "$jr_scripts_dir/base.sh"

bench_init

reps=${REPS:-5}
groups=${LIKWID_GROUPS:-"CACHES"}

threads_fixed=${THREADS:-${JR_PHYS_PER_SOCKET:-64}}
thread_list=${THREAD_LIST:-"1 2 4 8 16 32 64"}
batch_fixed=${BATCH_SIZE:-1024}
batch_list=${BATCH_LIST:-"16 64 256 1024 4096 16384 65536"}

bench_analyse_topology
bench_build_forward "${VARIANT:-base}"
bench_prepare_inputs
bench_check_groups "$groups"

if bench_cache_workingsets; then
  {
    echo "l1_ws_kb_per_thread=$JR_L1_WS_KB"
    echo "l2_ws_kb_per_thread=$JR_L2_WS_KB"
    echo "l3_ws_kb_per_thread=$JR_L3_WS_KB"
  } | tee -a "$JR_RUN_DIR/topology.txt"
else
  echo "WARNING: could not read cache sizes from sysfs -> skipping cache-size annotation." >&2
fi

for rep in $(seq 1 "$reps"); do

  # 1. Working-set sweep (fixed physical cores, growing batch size)
  cores_fixed=$(cpus_phys "$threads_fixed")
  for b in $batch_list; do
    current_batch=$b
    [ "$current_batch" -lt "$threads_fixed" ] && current_batch=$threads_fixed

    for group in "${JR_GROUPS[@]}"; do
      bench_run_forward "cache_workingset" "$threads_fixed" "$group" "$current_batch" "$rep" "$cores_fixed" "-m"
    done
  done

  # 2. Thread-count sweep (fixed batch size, growing physical-core count)
  for t in $thread_list; do
    current_batch=$batch_fixed
    [ "$current_batch" -lt "$t" ] && current_batch=$t
    cores=$(cpus_phys "$t")

    for group in "${JR_GROUPS[@]}"; do
      bench_run_forward "cache_threads" "$t" "$group" "$current_batch" "$rep" "$cores" "-m"
    done
  done
done

bench_finish
