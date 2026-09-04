#!/bin/bash
#SBATCH --account=slmet
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --time=01:30:00
#SBATCH --exclusive
#SBATCH --disable-perfparanoid
#SBATCH --job-name=e4_compare_ab

# Within a single Slurm allocation, run variant A, then run variant B on the exact same physical node
# Filter out noise by running on a single thread

# Output: out/{base,nomemset}.t<N>.<GROUP>.b<N>.rep<N>.csv

set -euo pipefail
 
JR_EXPERIMENT=e4_compare_ab
RUN_ID=${RUN_ID:-e1_compare_ab_${SLURM_JOB_ID:-manual}}

if [ -n "${SLURM_SUBMIT_DIR:-}" ] && [ -f "$SLURM_SUBMIT_DIR/base.sh" ]; then
  jr_scripts_dir="$SLURM_SUBMIT_DIR"
else
  jr_scripts_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi
source "$jr_scripts_dir/base.sh"
 
reps=${REPS:-5}
threads=${THREADS:-1}
batch=${BATCH_SIZE:-240}
groups=${LIKWID_GROUPS:-"MEM_DP CACHE"}
 
bench_init
bench_check_groups "$groups"
 
for variant in base nomemset; do
  bench_build "$variant"
 
  # Correctness gate per variant: a faster wrong answer is not an optimization.
  bench_validate
  cp -a "$JR_RUN_DIR/validation_status.txt" \
        "$JR_RUN_DIR/validation_status.${variant}.txt" 2>/dev/null || true
  cp -a "$JR_RUN_DIR/validation.log" \
        "$JR_RUN_DIR/validation.${variant}.log" 2>/dev/null || true
  if [ "$JR_VALIDATION_OK" -ne 1 ]; then
    echo "Variant '$variant' failed validation -> aborting A/B comparison." >&2
    exit 1
  fi
 
  bench_prepare_inputs
 
  for rep in $(seq 1 "$reps"); do
    for group in "${JR_GROUPS[@]}"; do
      bench_run "$variant" "$threads" "$group" "$batch" "$rep"
    done
  done
done
 
bench_finish
