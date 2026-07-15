#!/bin/bash
set -euo pipefail

script_dir=$(cd "$(dirname "$0")" && pwd)
repo_root=$(cd "$script_dir/../../.." && pwd)
src_dir="$repo_root/src"
runs_root="$repo_root/projects/benchmark/runs"
run_id=${RUN_ID:-local_$(date +%Y%m%d_%H%M%S)}
run_dir="$runs_root/$run_id"
work_dir="$run_dir/work"

case_name=${CASE_NAME:-${GEOMETRY:-zenith}_baseline}
baseline_cases="$repo_root/projects/benchmark/configs/baseline_cases.tsv"
case_row=$(awk -F'	' -v key="$case_name" 'NR > 1 && $1 == key { print; exit }' "$baseline_cases")
if [ -z "$case_row" ]; then
  echo "Unknown benchmark baseline case: $case_name" >&2
  exit 1
fi
geometry=$(printf '%s
' "$case_row" | awk -F'	' '{print $2}')
ctl_rel=$(printf '%s
' "$case_row" | awk -F'	' '{print $3}')
ctl_template=${CTLFILE:-$repo_root/$ctl_rel}
bench_tblbase=${BENCH_TBLBASE:-$HOME/wrk/jurassic/tab/tria_1cm/nc_1e-6/tria}
threads=${THREADS:-"1 2 4 8 12"}
cpu_batch_size=${CPU_BATCH_SIZE:-64}
compiler=${COMPILER:-gcc}
mpi=${MPI:-0}
rebuild=${REBUILD:-1}

mkdir -p "$work_dir"

if [ ! -f "$ctl_template" ]; then
  echo "Control file not found: $ctl_template" >&2
  exit 1
fi
bench_tbl_dir=$(dirname "$bench_tblbase")
if [ ! -d "$bench_tbl_dir" ]; then
  echo "Benchmark LUT directory not found: $bench_tbl_dir" >&2
  exit 1
fi

case "$geometry" in
  zenith|nadir|limb)
    ;;
  *)
    echo "Unsupported geometry: $geometry" >&2
    exit 1
    ;;
esac

maybe_plot() {
  local summary_tsv=$1
  local output_png=$2
  local title=$3

  if python3 -c "import matplotlib" >/dev/null 2>&1; then
    python3 "$script_dir/plot_benchmark_results.py" "$summary_tsv" --output "$output_png" --title "$title"
  else
    echo "WARNING: matplotlib not available; skipping plot generation for $output_png" >&2
  fi
}

cd "$work_dir"
export LANG=C
export LC_ALL=C

active_ctl="$work_dir/${case_name}.ctl"
awk -v tblbase="$bench_tblbase" '{ if ($1 == "TBLBASE") print "TBLBASE = " tblbase; else print $0; }' "$ctl_template" > "$active_ctl"

printf 'case_name=%s
geometry=%s
ctl_template=%s
active_ctl=%s
bench_tblbase=%s
threads=%s
cpu_batch_size=%s
compiler=%s
mpi=%s
rebuild=%s
'   "$case_name" "$geometry" "$ctl_template" "$active_ctl" "$bench_tblbase" "$threads" "$cpu_batch_size" "$compiler" "$mpi" "$rebuild" > "$run_dir/config.txt"

if [ "$rebuild" = 1 ]; then
  cd "$src_dir"
  make clean
  make -j MPI="$mpi" COMPILER="$compiler" GPU=0
  cd "$work_dir"
fi

rm -rf data
mkdir -p data

"$src_dir/climatology" "$active_ctl" data/atm.tab
"$src_dir/$geometry" "$active_ctl" data/obs.tab

for omp in $threads; do
  log="log.omp${omp}"
  out="/tmp/jurassic_bench_${run_id}_omp${omp}.tab"
  OMP_NUM_THREADS=$omp "$src_dir/formod" "$active_ctl" data/obs.tab data/atm.tab "$out" TASK time BATCH_SIZE "$cpu_batch_size" > "$log" 2>&1
  printf 'OMP_NUM_THREADS=%s
CPU_BATCH_SIZE=%s
' "$omp" "$cpu_batch_size" >> "$log"
done

python3 "$script_dir/summarize_time_logs.py" omp 'log.omp*' --tsv-out "$run_dir/summary.tsv" | tee "$run_dir/summary.md"
maybe_plot "$run_dir/summary.tsv" "$run_dir/plot_cpu_scaling.png" "JURASSIC local CPU baseline"
cp -a data "$run_dir/"
cp -a log.omp* "$run_dir/"
cp -a "$active_ctl" "$run_dir/"

echo "Run directory: $run_dir"
