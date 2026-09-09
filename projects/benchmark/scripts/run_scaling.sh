#!/bin/bash
#SBATCH --account=slmet
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --time=06:00:00
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

set -euo pipefail

JR_EXPERIMENT=e2_scaling
RUN_ID=${RUN_ID:-e2_scaling_${SLURM_JOB_ID:-manual}}

if [ -n "${SLURM_SUBMIT_DIR:-}" ] && [ -f "$SLURM_SUBMIT_DIR/base.sh" ]; then
  jr_scripts_dir="$SLURM_SUBMIT_DIR"
else
  jr_scripts_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi
source "$jr_scripts_dir/base.sh"
 
reps=${REPS:-3}
thread_list=${THREAD_LIST:-"1 2 4 8 12 24"}
smt_threads=${SMT_THREADS:-48}
batch=${BATCH_SIZE:-48}    # batch elements per thread
groups=${LIKWID_GROUPS:-"MEM_DP FLOPS_DP"}
 
bench_init
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
spread_threads=${SPREAD_THREADS:-$JR_N_PHYS}

for rep in $(seq 1 "$reps"); do

  # one socket, physical cores  
  for t in $thread_list; do
    for group in "${JR_GROUPS[@]}"; do
      bench_run_forward phys "$t" "$group" "$batch" "$rep" "$(cpus_phys "$t")"
    done
  done
 
  # one socket, SMT 
  for group in "${JR_GROUPS[@]}"; do
    bench_run_forward smt "$smt_threads" "$group" "$batch" "$rep" "$(cpus_smt "$smt_threads")"
  done

  # physical cores spread across sockets
  for t in $thread_list $spread_threads; do
    [ "$t" -le "$JR_N_PHYS" ] || continue
    [ "$t" -le "$JR_PHYS_PER_SOCKET" ] && continue
    for group in "${JR_GROUPS[@]}"; do
      bench_run_forward spread "$t" "$group" "$batch" "$rep" "$(cpus_spread "$t")"
    done
  done
done
 
bench_finish
