#!/bin/bash
#SBATCH --account=slmet
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --exclusive
#SBATCH --disable-perfparanoid
#SBATCH --job-name=e5_configsearch
#
# MPI/OpenMP configuration search
# Benchmark combinations of MPI RANKS x OPENMP THREADS PER RANK for fixed number of cores(=48)
# Output: out/retrieval.r<ranks>.t<threads>.<GROUP>.rep<N>.csv
# TODO: consider NUMA topology, report efficiency speedup/p

set -euo pipefail

JR_EXPERIMENT=e5_configsearch
RUN_ID=${RUN_ID:-e5_configsearch_${SLURM_JOB_ID:-manual}}

if [ -n "${SLURM_SUBMIT_DIR:-}" ] && [ -f "$SLURM_SUBMIT_DIR/base.sh" ]; then
  jr_scripts_dir="$SLURM_SUBMIT_DIR"
else
  jr_scripts_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi
source "$jr_scripts_dir/base.sh"

reps=${REPS:-3}
groups=${LIKWID_GROUPS:-"MEM_DP TMA"}
cores=${TOTAL_CORES:-48}
JR_RET_DIRLIST=${RET_DIRLIST:-"$JR_REPO_ROOT/projects/validation/dirlist.txt"}

# ranks * threads = cores: the divisor pairs of $cores
# yields configs: (1,48), (2,24), (4,12), (6,8), (8,6), (12,4), (16,3), (24,2), (48,1)
pairs=""
for r in $(seq 1 "$cores"); do
  if [ $(( cores % r )) -eq 0 ]; then
    pairs="$pairs $r:$(( cores / r ))"
  fi
done
pairs=${RANK_THREAD_PAIRS:-$pairs}

bench_init
bench_build_retrieval "${VARIANT:-base}"
bench_validate
bench_check_groups "$groups"

for rep in $(seq 1 "$reps"); do
  for pair in $pairs; do
    ranks=${pair%:*}
    threads=${pair#*:}
    for group in "${JR_GROUPS[@]}"; do
      bench_run_retrieval retrieval "$ranks" "$threads" "$group" "$rep"
    done
  done
done

bench_finish
cat <<'EOF'

Next step: for each (ranks,threads) pair, take the median wall-clock time
(PRINT_TIMERS TIMER_RET_KERNEL_INIT for the whole run, or wrap the retrieval
call in the shell with `time` if that timer isn't granular enough) and find
the minimum across the swept pairs -- that is the "best config" from E2's
question 3, now for the real production parallelization mechanism.
EOF