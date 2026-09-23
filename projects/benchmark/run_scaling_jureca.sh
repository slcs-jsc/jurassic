#!/bin/bash
#SBATCH --account=slmet
#SBATCH --partition=dc-cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=01:00:00
#SBATCH --exclusive
#SBATCH --disable-perfparanoid
#SBATCH --job-name=e2_scaling

# Benchmarking Script for JURASSIC (CPU-version)
# Target: JURECA DC -> AMD EPYC 7742 (2× 64 cores, 2.25 GHz)    
#
# Scaling Axes:
#   1. Intra-Socket: Scale across physical cores on Socket 0
#   2. SMT: Separate run using all 128 logical threads on Socket 0
#   3. Inter-Socket: Compare at a fixed core count (e.g., 64 threads):
#      - Compact: All threads executed on Socket 0
#      - Spread: Threads split evenly across Sockets
#
# Scaling Technique:
#   - Strong Scaling (SCALING_MODE="strong"): Total problem size stays constant. Runtime time directly equals "Time-to-Solution".
#   - Weak Scaling (SCALING_MODE="weak"): Workload per core stays constant. Total problem size scales linearly with the thread count.
#
# Metrics:
#   - Measure Runtime -> compute Speedup & Efficiency
#   - Analyze memory volume/bandwidth (MEM_DP)
#
# Thread pinning is handled via LIKWID hardware expressions (E:S0:...)
# Memory affinity is dynamically set via Membind (-m) or Interleave (-i)
#
# Output: out/<category>_<strong/weak>.t<threads>.<GROUP>.b<batch>.rep<rep>.csv

get_batch_size() {
  local threads=$1
  if [ "$SCALING_MODE" = "strong" ]; then
    echo "$BATCH_SIZE"
  elif [ "$SCALING_MODE" = "weak" ]; then
    echo $(( BATCH_SIZE * threads ))
  else
    echo "Error: Unknown SCALING_MODE '$SCALING_MODE'" >&2
    exit 1
  fi
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
batch=${BATCH_SIZE:-64}    

export JR_SCRIPTS_DIR_OVERRIDE="$jr_scripts_dir"
source "$jr_scripts_dir/base.sh"
 
bench_init

target_threads=${JR_PHYS_PER_SOCKET:-64}
reps=${REPS:-5}
thread_list=${THREAD_LIST:-"1 2 4 8 16 32 64"}
smt_threads=${SMT_THREADS:-128}
groups=${LIKWID_GROUPS:-"MEM_DP FLOPS_DP"}
 
bench_analyse_topology
bench_build_forward "${VARIANT:-base}"
bench_prepare_inputs
bench_check_groups "$groups"
 
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
spread_threads=${SPREAD_THREADS:-"$(( JR_PHYS_PER_SOCKET + 4 )) $(( JR_PHYS_PER_SOCKET * 2 ))"}

for rep in $(seq 1 "$reps"); do

  # Intra-Socket 
  numa_flag="-m" 
  for t in $thread_list; do
    current_batch=$(get_batch_size "$t")
    likwid_cores="E:S0:${t}"

    for group in "${JR_GROUPS[@]}"; do
      bench_run_forward "intra_socket_${SCALING_MODE}" "$t" "$group" "$current_batch" "$rep" "$likwid_cores" "$numa_flag"
    done
  done
 
  # SMT (Hyperthreading)
  current_batch=$(get_batch_size "$smt_threads")
  likwid_cores="E:S0:${smt_threads}"

  for group in "${JR_GROUPS[@]}"; do
    bench_run_forward smt_socket_${SCALING_MODE} "$smt_threads" "$group" "$current_batch" "$rep" "$likwid_cores" "-m"
  done

  # Inter-Socket
  current_batch=$(get_batch_size "$target_threads")

  # 1. All Threads on Socket 0
  cores_compact="E:S0:${target_threads}"
  for group in "${JR_GROUPS[@]}"; do
    bench_run_forward "inter_compact_${SCALING_MODE}" "$target_threads" "$group" "$current_batch" "$rep" "$cores_compact" "-i"
  done

  # 2. Spread threads evenly across sockets
  threads_per_sock=$(( target_threads / JR_N_SOCKETS ))

  if [ "$threads_per_sock" -gt 0 ]; then
    cores_spread=""
    for (( s=0; s<JR_N_SOCKETS; s++ )); do
      if [ -z "$cores_spread" ]; then
        cores_spread="E:S${s}:${threads_per_sock}"
      else
        cores_spread="${cores_spread}@E:S${s}:${threads_per_sock}"
      fi
    done

    for group in "${JR_GROUPS[@]}"; do
        bench_run_forward "inter_spread_${SCALING_MODE}" "$target_threads" "$group" "$current_batch" "$rep" "$cores_spread"
    done
  else
    echo "WARNING: target_threads ($target_threads) ist kleiner als die Anzahl der Sockets ($JR_N_SOCKETS). Full Spread wird übersprungen." >&2
  fi
done
 
bench_finish