#!/bin/bash
#SBATCH --account=slmet
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --time=02:00:00
#SBATCH --exclusive
#SBATCH --disable-perfparanoid
#SBATCH --job-name=e3_tma

# Top-down Microarchitecture Analysis (Retiring / Bad Speculation / Frontend Bound / Backend Bound)
# Reference: https://www.intel.com/content/www/us/en/docs/vtune-profiler/cookbook/2024-0/top-down-microarchitecture-analysis-method.html
# Measured with two different thread counts:
#   1 thread -> analyse intrinsic behavior of code
# 24 threads -> analyse what happens once all cores access the same memory path i.e. problems that arise with shared-resources

# Output: out/tma_<variant>.t<N>.<GROUP>.b<N>.rep<N>.csv

set -euo pipefail
 
JR_EXPERIMENT=e3_tma
RUN_ID=${RUN_ID:-e3_tma_${SLURM_JOB_ID:-manual}}

if [ -n "${SLURM_SUBMIT_DIR:-}" ] && [ -f "$SLURM_SUBMIT_DIR/base.sh" ]; then
  jr_scripts_dir="$SLURM_SUBMIT_DIR"
else
  jr_scripts_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi
source "$jr_scripts_dir/base.sh"

reps=${REPS:-3}
thread_list=${THREAD_LIST:-"1 24"}
per_thread=${WORK_PER_THREAD:-10}
groups=${LIKWID_GROUPS:-"TMA CACHE L3"}
variants=${VARIANTS:-"base nomemset"} 
  
bench_init
bench_check_groups "$groups"
 
for variant in $variants; do
  bench_build "$variant"
  bench_validate
  if [ "$JR_VALIDATION_OK" -ne 1 ]; then
    echo "Variant '$variant' failed validation -- skipping." >&2
    continue
  fi
  bench_prepare_inputs
 
  for rep in $(seq 1 "$reps"); do
    for t in $thread_list; do
      batch=$(( t * per_thread ))
      for group in "${JR_GROUPS[@]}"; do
        bench_run "tma_${variant}" "$t" "$group" "$batch" "$rep" "$(cores_phys "$t")"
      done
    done
  done
done
 
bench_finish
