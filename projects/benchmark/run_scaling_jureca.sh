#!/bin/bash
#SBATCH --account=slmet
#SBATCH --partition=dc-cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=04:00:00
#SBATCH --exclusive
#SBATCH --job-name=e1_scaling

# Benchmarking Script for JURASSIC (CPU-version)
# Target: JURECA DC -> AMD EPYC 7742 (2× 64 cores, 2.25 GHz)    
#
# Scaling Axes:
#   1. Intra-Socket: Scale across physical cores on Socket 0
#   2. SMT: Separate run using all logical threads on Socket 0
#   3. Inter-Socket: Compare at a fixed core count (one socket's worth of cores):
#      - Compact: All threads executed on Socket 0
#      - Spread: Threads split evenly across Sockets
#
# Scaling Technique:
#   - Strong Scaling (SCALING_MODE="strong"): Total problem size stays constant. Runtime time directly equals "Time-to-Solution".
#   - Weak Scaling (SCALING_MODE="weak"): Workload per core stays constant. Total problem size scales linearly with the thread count.
#
# Metrics:
#   - Measure Runtime -> compute Speedup & Efficiency
# Output: out/<category>_<strong/weak>.t<threads>.b<batch>.rep<rep>.txt (runtime only, no LIKWID)

get_batch_size() {
  local threads=$1
  local batch
  if [ "$SCALING_MODE" = "strong" ]; then
    batch="$BATCH_SIZE"
  elif [ "$SCALING_MODE" = "weak" ]; then
    batch=$(( BATCH_SIZE * threads ))
  else
    echo "Error: Unknown SCALING_MODE '$SCALING_MODE'" >&2
    exit 1
  fi
  
  if [ "$batch" -lt "$threads" ]; then
    batch="$threads"
  fi

  echo "$batch"
}

set -euo pipefail

JR_EXPERIMENT=e2_scaling
RUN_ID=${RUN_ID:-e2_scaling_${SLURM_JOB_ID:-manual}}

if [ -n "${SLURM_SUBMIT_DIR:-}" ] && [ -f "$SLURM_SUBMIT_DIR/base.sh" ]; then
  jr_scripts_dir="$SLURM_SUBMIT_DIR"
else
  jr_scripts_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi

# Strong scaling: batch elements per thread, Weak scaling: size for single-thread
SCALING_MODE=${SCALING_MODE:-"strong"}
BATCH_SIZE=${BATCH_SIZE:-64}    

export JR_SCRIPTS_DIR_OVERRIDE="$jr_scripts_dir"
JR_USE_LIKWID=0
source "$jr_scripts_dir/base.sh"
 
bench_init

target_threads=${JR_PHYS_PER_SOCKET:-64}
reps=${REPS:-5}
thread_list=${THREAD_LIST:-"1 2 4 8 16 32 64"}
smt_threads=${SMT_THREADS:-128}
 
bench_analyse_topology
bench_build_forward "${VARIANT:-base}"
bench_prepare_inputs
 
if [ -z "${THREAD_LIST:-}" ]; then
  thread_list=""
  t=1
  while [ "$t" -le "$JR_PHYS_PER_SOCKET" ]; do
    thread_list="$thread_list $t"
    t=$(( t * 2 ))
  done
  [ "$t" -ne $(( JR_PHYS_PER_SOCKET * 2 )) ] && thread_list="$thread_list $JR_PHYS_PER_SOCKET"
fi
smt_threads=${SMT_THREADS:-$(( JR_PHYS_PER_SOCKET * JR_SMT ))}

for rep in $(seq 1 "$reps"); do

  # Intra-Socket
  for t in $thread_list; do
    bench_run_time "intra_socket_${SCALING_MODE}" "$t" "$(get_batch_size "$t")" "$rep" "$(cpus_phys "$t")"
  done

  # SMT: all logical CPUs of socket 0
  bench_run_time "smt_socket_${SCALING_MODE}" "$smt_threads" "$(get_batch_size "$smt_threads")" "$rep" "$(cpus_smt "$smt_threads")"

  # Inter-Socket, same thread count: all on socket 0 vs. spread evenly across sockets
  current_batch=$(get_batch_size "$target_threads")
  bench_run_time "inter_compact_${SCALING_MODE}" "$target_threads" "$current_batch" "$rep" "$(cpus_phys "$target_threads")"

  if [ $(( target_threads / JR_N_SOCKETS )) -gt 0 ]; then
    bench_run_time "inter_spread_${SCALING_MODE}" "$target_threads" "$current_batch" "$rep" "$(cpus_spread "$target_threads")"
  else
    echo "WARNING: target_threads ($target_threads) < sockets ($JR_N_SOCKETS), skipping spread run." >&2
  fi
done
 
bench_finish
