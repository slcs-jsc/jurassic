#!/bin/bash
#SBATCH --account=slmet
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --time=02:00:00
#SBATCH --exclusive
#SBATCH --disable-perfparanoid
#SBATCH --job-name=e1_roofline

# Performs LIKWID profiling (FLOPS_DP, MEM_DP) to measure operational intensity 
# i.e. the ratio of computational work to data movement: 
#     operational intensity = DP [MFLOP/s] / (Memory bandwidth [MBytes/s])
#                           = FLOPs per byte moved
# Goal: analyse whether algorithm is memory- or compute-bound. Cache-residency crossover
# Reference: https://github.com/RRZE-HPC/likwid/wiki/Tutorial:-Empirical-Roofline-Model

# Output: out/size_<case>.t1.<GROUP>.b<N>.rep<N>.csv

set -euo pipefail
 
JR_EXPERIMENT=e1_roofline
RUN_ID=${RUN_ID:-e1_roofline_${SLURM_JOB_ID:-manual}}

if [ -n "${SLURM_SUBMIT_DIR:-}" ] && [ -f "$SLURM_SUBMIT_DIR/base.sh" ]; then
  jr_scripts_dir="$SLURM_SUBMIT_DIR"
else
  jr_scripts_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi
source "$jr_scripts_dir/base.sh"
 
reps=${REPS:-3}
threads=${THREADS:-1}
batch=${BATCH_SIZE:-240}
groups=${LIKWID_GROUPS:-"FLOPS_DP MEM_DP CACHE"}
 
# Cases to sweep. Defaults to the three geometries: zenith, nadir, limb
# Can be overwritten with a list of case names from configs/baseline_cases.tsv 
# that differ in channel count, gas count or lookup-table size.
case_list=${CASE_LIST:-"zenith_baseline nadir_baseline limb_baseline"}
 
first=1
for case_name in $case_list; do
  export CASE_NAME="$case_name"
  export RUN_ID="${RUN_ID%%_case_*}"
 
  bench_init
  if [ "$first" -eq 1 ]; then
    bench_build "${VARIANT:-base}"
    bench_validate
    bench_check_groups "$groups"
    first=0
  fi
  bench_prepare_inputs
 
  for rep in $(seq 1 "$reps"); do
    for group in "${JR_GROUPS[@]}"; do
      bench_run "size_${case_name}" "$threads" "$group" "$batch" "$rep"
    done
  done
done
 
bench_finish
