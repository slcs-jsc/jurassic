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
# Scale BATCH_SIZE with thread count so the work per thread stays constant

# Goal: Measure how runtime, calls and MEM_DP volumes compare scale with threadcount

# Output: out/{phys,smt}.t<N>.<GROUP>.b<N>.rep<N>.csv

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
per_thread=${WORK_PER_THREAD:-10}     # batch elements per thread
groups=${LIKWID_GROUPS:-"MEM_DP FLOPS_DP"}
 
bench_init
bench_build "${VARIANT:-base}"
bench_validate
bench_prepare_inputs
bench_check_groups "$groups"
 
# physical cores
for rep in $(seq 1 "$reps"); do
  for t in $thread_list; do
    batch=$(( t * per_thread ))
    for group in "${JR_GROUPS[@]}"; do
      bench_run phys "$t" "$group" "$batch" "$rep" "$(cores_phys "$t")"
    done
  done
done
 
# SMT data point
for rep in $(seq 1 "$reps"); do
  batch=$(( smt_threads * per_thread ))
  for group in "${JR_GROUPS[@]}"; do
    bench_run smt "$smt_threads" "$group" "$batch" "$rep" "$(cores_smt "$smt_threads")"
  done
done
 
bench_finish

