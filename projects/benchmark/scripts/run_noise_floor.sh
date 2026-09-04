#!/bin/bash
#SBATCH --account=slmet
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --time=01:00:00
#SBATCH --exclusive
#SBATCH --disable-perfparanoid
#SBATCH --job-name=e0_noise

# Runs identical binary + config for N times 
# Goal: determine measurement noise 

# Output: out/noise.t1.<GROUP>.b<N>.rep<1..N>.csv


set -euo pipefail
set -x                               
trap 'echo "FAILED at line $LINENO: $BASH_COMMAND" >&2' ERR 

JR_EXPERIMENT=e0_noise
RUN_ID=${RUN_ID:-e0_noise_${SLURM_JOB_ID:-manual}}

if [ -n "${SLURM_SUBMIT_DIR:-}" ] && [ -f "$SLURM_SUBMIT_DIR/base.sh" ]; then
  jr_scripts_dir="$SLURM_SUBMIT_DIR"
else
  jr_scripts_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi
source "$jr_scripts_dir/base.sh"

reps=${REPS:-12}
threads=${THREADS:-1}
batch=${BATCH_SIZE:-240}
groups=${LIKWID_GROUPS:-"MEM_DP"}

bench_init
bench_build "${VARIANT:-base}"
bench_validate
bench_prepare_inputs
bench_check_groups "$groups"

for rep in $(seq 1 "$reps"); do
  for group in "${JR_GROUPS[@]}"; do
    bench_run noise "$threads" "$group" "$batch" "$rep"
  done
done

bench_finish