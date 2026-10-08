#!/bin/bash
#SBATCH --account=slmet
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --time=01:45:00
#SBATCH --exclusive
#SBATCH --disable-perfparanoid
#SBATCH --job-name=hermes_compare

# Profile the unmodified baseline and the optimized code on the same node with
# run_hermes_profile.sh, then plot both reports in one diagram.
#
# Output: runs/<RUN_ID>_{baseline,optimized}/   one hermes profile report each
#         runs/<RUN_ID>_compare/                plot (png + pdf) and summary.tsv
#
# Environment (all optional; the rest is passed through to run_hermes_profile.sh):
#   BASELINE_SRC   source tree of the baseline  (default: projects/benchmark/baseline/src)
#   OPTIMIZED_SRC  source tree of the changes   (default: src)
#   RUN_ID         prefix of the run directories (default: hermes_compare_<job id>)
#   CASE_NAME, LIKWID_THREADS, LIKWID_GROUPS, BENCH_TBLBASE, ...  see run_hermes_profile.sh

set -euo pipefail
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

baseline_src=${BASELINE_SRC:-$repo_root/projects/benchmark/baseline/src}
optimized_src=${OPTIMIZED_SRC:-$repo_root/src}
runs_root="$repo_root/projects/benchmark/runs"
run_id=${RUN_ID:-hermes_compare_${SLURM_JOB_ID:-manual}}

for d in "$baseline_src" "$optimized_src"; do
  if [ ! -f "$d/formod.c" ]; then
    echo "Not a JURASSIC source tree: $d" >&2
    exit 1
  fi
done

# Baseline first, then the changes, on the same node and with identical settings.
# Validation gates only the optimized code: the baseline is the reference, and
# run_validation.py would otherwise check the optimized src/formod either way.
for variant in baseline optimized; do
  if [ "$variant" = baseline ]; then
    src=$baseline_src
    skip_validation=1
  else
    src=$optimized_src
    skip_validation=${SKIP_VALIDATION:-0}
  fi

  echo "=== Profiling $variant ($src) ==="
  RUN_ID="${run_id}_${variant}" SRC_DIR="$src" BUILD_VERSION="$variant" \
    SKIP_VALIDATION="$skip_validation" \
    bash "$script_dir/run_hermes_profile.sh"
done

# Plot both reports in one diagram.
if command -v ml >/dev/null 2>&1; then
  ml Stages/2026 GCC/14.3.0 SciPy-bundle/2025.07 || true
fi
python3 "$repo_root/projects/benchmark/utils/eval_hermes_compare.py" \
  "$runs_root/${run_id}_baseline" "$runs_root/${run_id}_optimized" \
  --out "$runs_root/${run_id}_compare"

echo "Comparison: $runs_root/${run_id}_compare"
