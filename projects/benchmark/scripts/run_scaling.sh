#!/bin/bash
#SBATCH --account=slmet
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --time=02:00:00
#SBATCH --exclusive
#SBATCH --disable-perfparanoid
#SBATCH --job-name=e2_scaling

# Physical Core Scaling with (1, 2, 4, 8, 12, 24)
# Record SMT (48 logical CPUs on socket 0) separately 

# Strong Scaling: BATCH_SIZE is fixed across all thread counts, so every configuration solves the same total problem and 
# wall-clock time is directly comparable (time-to-solution). 

# Goal: Measure how runtime, calls and MEM_DP volumes compare scale with threadcount

# Output: out/{phys,smt}.t<N>.<GROUP>.b<N>.rep<N>.csv

# TODO: consider NUMA topology, report efficiency speedup/p

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
batch=${BATCH_SIZE:-48}    

source "$jr_scripts_dir/base.sh"
 
bench_init

reps=${REPS:-3}
thread_list=${THREAD_LIST:-"1 2 4 8 12 24"}
smt_threads=${SMT_THREADS:-48}
groups=${LIKWID_GROUPS:-"MEM_DP FLOPS_DP"}
 
bench_analyse_topology
bench_build_forward "${VARIANT:-base}"
bench_validate
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
  target_threads=$JR_PHYS_PER_SOCKET
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
