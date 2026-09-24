#!/bin/bash
#SBATCH --account=slmet
#SBATCH --partition=dc-cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=03:00:00
#SBATCH --exclusive
#SBATCH --disable-perfparanoid
#SBATCH --job-name=cpu_baseline

# CPU baseline thread sweep for JURASSIC (TASK time) on a full node.
#
# Scaling modes (SCALING_MODE):
#   strong: BATCH_SIZE = total number of forward models (same at every thread count)
#   weak:   BATCH_SIZE = forward models per thread (total = BATCH_SIZE * threads)
#
# Threads fill socket 0 first, then socket 1 (compact)
#
# Thread list, in order of precedence:
#   THREAD_LIST="1 2 4 ..."            explicit list
#   SWEEP_MATRIX=<script> [SWEEP_FAMILY=<family>]
#                                      omp_threads of the cpu rows of the sweep
#                                      matrix for the selected geometry. Ray set,
#                                      band, channel count and gas set of the
#                                      matrix rows are NOT applied here.
#   default                            powers of two up to all physical cores
#
# Output: out/cpu_<mode>.t<threads>.<GROUP>.b<batch>.rep<rep>.{csv,txt}
#         summary.cpu.tsv (aggregated over reps), summary.cpu.per_rep.tsv,
#         summary.cpu.txt, plot_cpu_scaling.png

set -euo pipefail

JR_EXPERIMENT=cpu_baseline
RUN_ID=${RUN_ID:-cpu_baseline_${SLURM_JOB_ID:-manual}}

if [ -n "${SLURM_SUBMIT_DIR:-}" ] && [ -f "$SLURM_SUBMIT_DIR/base.sh" ]; then
  jr_scripts_dir="$SLURM_SUBMIT_DIR"
else
  jr_scripts_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi

SCALING_MODE=${SCALING_MODE:-strong}
case "$SCALING_MODE" in
  strong) BATCH_SIZE=${BATCH_SIZE:-1024} ;;
  weak)   BATCH_SIZE=${BATCH_SIZE:-8} ;;
  *) echo "Unknown SCALING_MODE '$SCALING_MODE' (strong|weak)" >&2; exit 1 ;;
esac
reps=${REPS:-3}
group=${LIKWID_GROUP:-FLOPS_DP}
variant=${VARIANT:-base}

# Fixed sample count per run
export JURASSIC_MAX_ITER=${SAMPLES:-10}
export JURASSIC_TIME_BUDGET=${JURASSIC_TIME_BUDGET:-1e9}

export JR_SCRIPTS_DIR_OVERRIDE="$jr_scripts_dir"
source "$jr_scripts_dir/base.sh"

bench_init                       
bench_build_forward "$variant"
bench_prepare_inputs
bench_check_groups "$group"

maybe_plot() {
  local summary_tsv=$1
  local output_png=$2
  local title=$3

  if python3 -c "import matplotlib" >/dev/null 2>&1; then
    python3 "$JR_SCRIPT_DIR/plot_benchmark_results.py" "$summary_tsv" \
      --output "$output_png" --title "$title"
  else
    echo "WARNING: matplotlib not available; skipping plot generation for $output_png" >&2
  fi
}

# thread list
if [ -n "${THREAD_LIST:-}" ]; then
  thread_list=$THREAD_LIST
elif [ -n "${SWEEP_MATRIX:-}" ]; then
  # Columns of the sweep matrix: family geometry ray_set nr band nd gas_set
  # target omp_threads batch_size
  thread_list=$(python3 "$SWEEP_MATRIX" "${SWEEP_FAMILY:-cpu_strong_scaling}" \
    | awk -F'\t' -v geo="$JR_GEOMETRY" \
        'NR > 1 && $2 == geo && $8 == "cpu" && $9 != "" { print $9 }' \
    | sort -nu | tr '\n' ' ')
  if [ -z "$thread_list" ]; then
    echo "Sweep matrix has no cpu rows for family '${SWEEP_FAMILY:-cpu_strong_scaling}'," \
         "geometry '$JR_GEOMETRY'" >&2
    exit 1
  fi
else
  thread_list=""
  t=1
  while [ "$t" -lt "$JR_N_PHYS" ]; do
    thread_list="$thread_list $t"
    t=$(( t * 2 ))
  done
  thread_list="$thread_list $JR_N_PHYS"
fi

max_t=0
for t in $thread_list; do
  [ "$t" -gt "$JR_N_PHYS" ] && {
    echo "Thread count $t exceeds the $JR_N_PHYS physical cores" >&2; exit 1; }
  [ "$t" -gt "$max_t" ] && max_t=$t
done

if [ "$SCALING_MODE" = strong ] && [ "$BATCH_SIZE" -lt $(( 2 * max_t )) ]; then
  echo "WARNING: strong scaling with BATCH_SIZE=$BATCH_SIZE gives < 2 elements" \
       "per thread at $max_t threads; use a multiple of $max_t (>= 4x)." >&2
fi
if [ "$SCALING_MODE" = weak ] && [ "$BATCH_SIZE" -lt 2 ]; then
  echo "BATCH_SIZE per thread must be >= 2 in weak mode" \
       "(batch 1 takes the scalar path in formod.c)" >&2
  exit 1
fi

get_batch_size() {
  local threads=$1
  if [ "$SCALING_MODE" = strong ]; then
    echo "$BATCH_SIZE"
  else
    echo $(( BATCH_SIZE * threads ))
  fi
}

# Compact placement: fill socket 0, then socket 1, ...
cores_expr() {
  local remaining=$1 p=$JR_PHYS_PER_SOCKET expr="" n s
  for (( s=0; s<JR_N_SOCKETS && remaining>0; s++ )); do
    n=$(( remaining < p ? remaining : p ))
    expr="${expr:+$expr@}E:S${s}:${n}"
    remaining=$(( remaining - n ))
  done
  echo "$expr"
}

{
  echo "scaling_mode=$SCALING_MODE"
  echo "batch_size_param=$BATCH_SIZE"
  echo "thread_list=$thread_list"
  echo "reps=$reps"
  echo "likwid_group=$group"
  echo "variant=$variant"
  echo "jurassic_max_iter=$JURASSIC_MAX_ITER"
  echo "jurassic_time_budget=$JURASSIC_TIME_BUDGET"
} >> "$JR_RUN_DIR/config.txt"

# sweep
for rep in $(seq 1 "$reps"); do
  for t in $thread_list; do
    batch=$(get_batch_size "$t")
    bench_run_forward "cpu_${SCALING_MODE}" "$t" "$group" "$batch" "$rep" \
                      "$(cores_expr "$t")" ""
  done
done

# summary
for f in "$JR_WORK_DIR"/out/cpu_"${SCALING_MODE}".t*.txt; do
  [ -e "$f" ] || continue
  n=$(sed -n 's/^threads=//p' "$f" | tail -n1)
  if [ -n "$n" ]; then
    echo "OMP_NUM_THREADS=$n" >> "$f"
  fi
done

( cd "$JR_WORK_DIR/out" \
  && python3 "$JR_SCRIPT_DIR/summarize_time_logs.py" omp \
       "cpu_${SCALING_MODE}.t*.txt" \
       --tsv-out "$JR_RUN_DIR/summary.cpu.tsv" \
       --per-rep-tsv-out "$JR_RUN_DIR/summary.cpu.per_rep.tsv" ) \
  | tee "$JR_RUN_DIR/summary.cpu.txt"

maybe_plot "$JR_RUN_DIR/summary.cpu.tsv" "$JR_RUN_DIR/plot_cpu_scaling.png" \
           "JURASSIC CPU baseline (${SCALING_MODE} scaling)"

bench_finish