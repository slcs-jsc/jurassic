#!/bin/bash
#SBATCH --account=slmet
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --time=01:00:00
#SBATCH --exclusive
#SBATCH --disable-perfparanoid

# LIKWID profile of the retrieval (optimal_estimation) for the hermes
# optimization loop -- the retrieval counterpart of run_hermes_profile.sh.
#
# The case is built so every run does the same amount of work:
#   - measurements are simulated from the climatology ("truth"),
#   - the a priori is the truth with temperature shifted by RET_DT_APR K, so
#     Levenberg-Marquardt actually iterates,
#   - CONV_DMIN=0 disables early convergence, so exactly RET_CONV_ITMAX outer
#     iterations (and a fixed number of kernel recomputations) run,
#   - RET_CASES identical copies of the case are retrieved in one process;
#     retrieval prints per-case time statistics as a RUNTIME: line.
# The number of rejected LM steps can still change if a code change alters
# the numerics; the ret_lm_formod call count in the CSV shows when it does.
#
# Output per thread count/group: log.omp<N>.<GROUP>.{csv,txt}, marker regions
# retrieval, kernel_jacobian, formod, ret_lm_formod, ret_linalg, ret_err_ana.

set -euo pipefail
set -x
trap 'echo "FAILED at line $LINENO: $BASH_COMMAND" >&2' ERR

script_source=${BASH_SOURCE[0]:-$0}
script_dir=$(cd "$(dirname "$script_source")" && pwd)
repo_root=$(cd "$script_dir/../../.." && pwd)
if [ ! -f "$repo_root/projects/benchmark/configs/baseline_cases.tsv" ] && [ -n "${SLURM_SUBMIT_DIR:-}" ]; then
  submit_repo_root=$(cd "$SLURM_SUBMIT_DIR/../../.." && pwd)
  if [ -f "$submit_repo_root/projects/benchmark/configs/baseline_cases.tsv" ]; then
    repo_root=$submit_repo_root
    script_dir="$repo_root/projects/benchmark/scripts"
  fi
fi
src_dir="$repo_root/src"
runs_root="$repo_root/projects/benchmark/runs"
run_id=${RUN_ID:-hermes_retrieval_profile_${SLURM_JOB_ID:-manual}}
run_dir="$runs_root/$run_id"
work_dir="$run_dir/work"

case_name=${CASE_NAME:-limb_baseline}
baseline_cases="$repo_root/projects/benchmark/configs/baseline_cases.tsv"
case_row=$(awk -F'\t' -v key="$case_name" 'NR > 1 && $1 == key { print; exit }' "$baseline_cases")
if [ -z "$case_row" ]; then
  echo "Unknown benchmark baseline case: $case_name" >&2
  exit 1
fi
geometry=$(printf '%s\n' "$case_row" | awk -F'\t' '{print $2}')
ctl_rel=$(printf '%s\n' "$case_row" | awk -F'\t' '{print $3}')
ctl_template=${CTLFILE:-$repo_root/$ctl_rel}
bench_tblbase=${BENCH_TBLBASE:-/p/data1/slmet/model_data/jurassic/tab/tria_1cm/nc_1e-6/tria}
compiler_cpu=${COMPILER_CPU:-gcc}
rebuild=${REBUILD:-1}

ret_cases=${RET_CASES:-3}
ret_conv_itmax=${RET_CONV_ITMAX:-3}
ret_kernel_recomp=${RET_KERNEL_RECOMP:-3}
ret_err_ana=${RET_ERR_ANA:-1}
ret_dt_apr=${RET_DT_APR:-3}

likwid_threads=${LIKWID_THREADS:-"24"}
likwid_groups=${LIKWID_GROUPS:-"MEM_DP FLOPS_DP"}
likwid_socket=${LIKWID_SOCKET:-0}

mkdir -p "$work_dir"
echo "perf_event_paranoid: $(cat /proc/sys/kernel/perf_event_paranoid 2>/dev/null || echo unavailable)" > "$run_dir/perf_paranoid_status.txt"

if [ ! -f "$ctl_template" ]; then
  echo "Control file not found: $ctl_template" >&2
  exit 1
fi
if [ ! -d "$(dirname "$bench_tblbase")" ]; then
  echo "Benchmark LUT directory not found: $(dirname "$bench_tblbase")" >&2
  exit 1
fi
case "$geometry" in
  zenith|nadir|limb) ;;
  *) echo "Unsupported geometry: $geometry" >&2; exit 1 ;;
esac

cd "$work_dir"
export LANG=C
export LC_ALL=C

if command -v ml >/dev/null 2>&1; then
  ml Stages/2026 GCC/14.3.0
  ml likwid/5.4.1
  ml CMake/4.0.3
  ml ecBuild
  ml SciPy-bundle/2025.07
  ml netcdf4-python/1.7.2
fi
if ! command -v likwid-perfctr >/dev/null 2>&1; then
  echo "likwid-perfctr not found on PATH after 'ml likwid'." >&2
  exit 1
fi
likwid-perfctr -a > "$run_dir/likwid_available_groups.txt" 2>&1 || true

export LD_LIBRARY_PATH="$repo_root/libs/build/lib:$repo_root/libs/build/lib64:${LD_LIBRARY_PATH:-}"

active_ctl="$work_dir/${case_name}_retrieval.ctl"
awk -v tblbase="$bench_tblbase" '{ if ($1 == "TBLBASE") print "TBLBASE = " tblbase; else print $0; }' "$ctl_template" > "$active_ctl"

# Fixed-work retrieval settings, passed as KEY VALUE overrides. The case ctl
# files only define the retrieval targets (RET*_ZMIN/ZMAX); all a priori and
# measurement errors default to 0, which optimal_estimation() rejects.
ret_args=(CONV_ITMAX "$ret_conv_itmax" KERNEL_RECOMP "$ret_kernel_recomp"
          CONV_DMIN 0 ERR_ANA "$ret_err_ana" WRITE_MATRIX 0
          ERR_PRESS 10 ERR_PRESS_CZ 5 ERR_PRESS_CH 200
          ERR_TEMP 5 ERR_TEMP_CZ 5 ERR_TEMP_CH 200
          "ERR_Q[*]" 50 "ERR_Q_CZ[*]" 5 "ERR_Q_CH[*]" 200
          "ERR_K[*]" 1e-3 "ERR_K_CZ[*]" 5 "ERR_K_CH[*]" 200
          "ERR_NOISE[*]" 1e-5 "ERR_FORMOD[*]" 1)

lscpu > "$run_dir/lscpu.txt" 2>/dev/null || true
numactl --hardware > "$run_dir/numactl_hardware.txt" 2>/dev/null || true
if git -C "$repo_root" rev-parse --is-inside-work-tree >/dev/null 2>&1; then
  {
    echo "commit=$(git -C "$repo_root" rev-parse HEAD)"
    echo "branch=$(git -C "$repo_root" rev-parse --abbrev-ref HEAD)"
    if [ -z "$(git -C "$repo_root" status --porcelain)" ]; then echo "dirty=0"; else echo "dirty=1"; fi
  } > "$run_dir/git_info.txt"
fi

printf 'case_name=%s\ngeometry=%s\nctl_template=%s\nactive_ctl=%s\nbench_tblbase=%s\nret_cases=%s\nret_args=%s\nret_dt_apr=%s\nlikwid_threads=%s\nlikwid_groups=%s\n' \
  "$case_name" "$geometry" "$ctl_template" "$active_ctl" "$bench_tblbase" \
  "$ret_cases" "${ret_args[*]}" "$ret_dt_apr" "$likwid_threads" "$likwid_groups" \
  > "$run_dir/config.txt"

# libs/build (GSL, netCDF, HDF5, ...) is excluded from the optimization loop's
# rsync, so a fresh remote work dir has none. Seed it from LIBS_BUILD_DIR (an
# existing libs/build, e.g. another worktree's); rsync leaves the excluded copy
# alone afterwards. Building the bundled libs (with their test suites) is the
# slow last resort.
ensure_libs() {
  local libs_dir="$repo_root/libs/build"
  [ -f "$libs_dir/include/gsl/gsl_math.h" ] && return 0
  if [ -n "${LIBS_BUILD_DIR:-}" ] && [ -f "$LIBS_BUILD_DIR/include/gsl/gsl_math.h" ]; then
    echo "Seeding $libs_dir from LIBS_BUILD_DIR=$LIBS_BUILD_DIR"
    mkdir -p "$libs_dir"
    cp -a "$LIBS_BUILD_DIR/." "$libs_dir/"
  else
    echo "No compiled libs in $libs_dir and LIBS_BUILD_DIR='${LIBS_BUILD_DIR:-}' has none either -- building bundled libs (slow)." >&2
    ( cd "$repo_root/libs" && bash build.sh > "$run_dir/libs_build.log" 2>&1 )
  fi
}

ensure_libs

if [ "$rebuild" = 1 ]; then
  ( cd "$src_dir" && make clean && make -j MPI=0 COMPILER="$compiler_cpu" GPU=0 LIKWID=1 )
fi

# Forward-model validation, as in run_hermes_profile.sh: the retrieval is only
# as correct as the forward model and kernel it is built on.
validation_status="$run_dir/validation_status.txt"
run_profiling=1
skip_reason=""
if [ "${SKIP_VALIDATION:-0}" != "1" ]; then
  trap - ERR
  set +e
  ( cd "$repo_root/projects/validation" && \
    VALIDATION_TBLBASE="$bench_tblbase" scripts/run_validation.py \
      > "$run_dir/validation.log" 2>&1 )
  validation_rc=$?
  set -e
  trap 'echo "FAILED at line $LINENO: $BASH_COMMAND" >&2' ERR
  echo "exit_code=$validation_rc" > "$validation_status"
  if [ "$validation_rc" -ne 0 ]; then
    run_profiling=0
    skip_reason="forward-model validation failed (exit $validation_rc)"
  fi
else
  echo "exit_code=skipped" > "$validation_status"
fi

# Inputs: truth -> simulated measurements; perturbed a priori; RET_CASES copies.
prepare_inputs() {
  rm -rf data
  mkdir -p data
  "$src_dir/climatology" "$active_ctl" data/atm_true.tab
  "$src_dir/$geometry" "$active_ctl" data/obs.tab
  "$src_dir/formod" "$active_ctl" data/obs.tab data/atm_true.tab data/obs_meas.tab
  awk -v dt="$ret_dt_apr" '/^#/ || NF == 0 { print; next } { $6 = $6 + dt; print }' \
    data/atm_true.tab > data/atm_apr.tab
  : > data/dirlist.txt
  for ic in $(seq 1 "$ret_cases"); do
    mkdir -p "data/case$ic"
    cp data/atm_apr.tab data/obs_meas.tab "data/case$ic/"
    echo "data/case$ic" >> data/dirlist.txt
  done
}

mkdir -p likwid
cd likwid
prepare_inputs

if [ "$run_profiling" = 1 ]; then
  valid_groups=()
  for group in $likwid_groups; do
    if grep -q -w "$group" "$run_dir/likwid_available_groups.txt"; then
      valid_groups+=("$group")
    else
      echo "WARNING: LIKWID group '$group' not available on this node -- skipping." >&2
    fi
  done
  if [ ${#valid_groups[@]} -eq 0 ]; then
    run_profiling=0
    skip_reason="none of the requested LIKWID groups ($likwid_groups) is available on this node"
  fi
fi

if [ "$run_profiling" = 1 ]; then
  unset OMP_PLACES OMP_PROC_BIND
  for omp in $likwid_threads; do
    core_list="S${likwid_socket}:0-$((omp - 1))"
    for group in "${valid_groups[@]}"; do
      log_txt="log.omp${omp}.${group}.txt"
      log_csv="log.omp${omp}.${group}.csv"
      echo "Running retrieval: LIKWID group=$group OMP_NUM_THREADS=$omp ..."
      trap - ERR
      set +e
      OMP_NUM_THREADS=$omp likwid-perfctr -C "$core_list" -g "$group" -m -o "$log_csv" \
        "$src_dir/retrieval" "$active_ctl" data/dirlist.txt "${ret_args[@]}" \
        > "$log_txt" 2>&1
      rc=$?
      set -e
      trap 'echo "FAILED at line $LINENO: $BASH_COMMAND" >&2' ERR
      printf 'OMP_NUM_THREADS=%s\nLIKWID_GROUP=%s\nCORE_LIST=%s\nexit_code=%s\n' \
        "$omp" "$group" "$core_list" "$rc" >> "$log_txt"
      if [ "$rc" -ne 0 ]; then
        skip_reason="retrieval exited with code $rc (OMP_NUM_THREADS=$omp, group $group); see $log_txt"
      elif grep -qiE "chi\^2/m= *-?(nan|inf)" "$log_txt"; then
        skip_reason="retrieval produced a non-finite chi^2 (OMP_NUM_THREADS=$omp, group $group)"
      fi
    done
  done
fi

if [ -n "$skip_reason" ]; then
  echo "Skipped/invalid retrieval profiling: $skip_reason" > "$run_dir/skipped_profiling.txt"
fi

cp -a data "$run_dir/data.likwid"
cp -a log.omp*.txt log.omp*.csv "$run_dir/" 2>/dev/null || true
cp -a "$active_ctl" "$run_dir/"

echo "LIKWID run directory: $run_dir"

if [ "${validation_rc:-0}" -ne 0 ]; then
  echo "Exiting non-zero: validation failed (see $validation_status)." >&2
  exit 3
fi
