#!/bin/bash
#SBATCH --account=slmet
#SBATCH --partition=booster
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --gpus-per-task=1
#SBATCH --disable-dcgm
#SBATCH --disable-perfparanoid
#SBATCH --time=00:45:00
#SBATCH --exclusive

set -euo pipefail

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
run_id=${RUN_ID:-juwels_ncu_test_${SLURM_JOB_ID:-manual}}
run_dir="$runs_root/$run_id"
work_dir="$run_dir/work"
out="$run_dir/formod.tab"

case_name="zenith_baseline"
baseline_cases="$repo_root/projects/benchmark/configs/baseline_cases.tsv"

case_row=$(awk -F'	' -v key="$case_name" 'NR > 1 && $1 == key { print; exit }' "$baseline_cases")
if [ -z "$case_row" ]; then
  echo "Unknown benchmark baseline case: $case_name" >&2
  exit 1
fi

geometry=$(printf '%s\n' "$case_row" | awk -F'	' '{print $2}')
ctl_rel=$(printf '%s\n' "$case_row" | awk -F'	' '{print $3}')
ctl_template=${CTLFILE:-$repo_root/$ctl_rel}

bench_tblbase=${BENCH_TBLBASE:-/p/data1/slmet/model_data/jurassic/tab/tria_1cm/nc_1e-6/tria}
slurm_cpus_per_task=${SLURM_CPUS_PER_TASK:-12}

compiler_gpu=${COMPILER_GPU:-nvc}
mpicc=${MPICC:-mpicc}
mpi=${MPI:-0}
gpu_pin=${GPU_PIN:-1}
gpu_target=${GPU_TARGET:-gpu}
info=${INFO:-0}

# Einzelne Test-Batchgröße für Nsight Compute festlegen
nvidia_profile_batch=${PROFILE_BATCH:-256}
nvidia_profile_output="$run_dir/ncu"

mkdir -p "$work_dir"
mkdir -p "$nvidia_profile_output"

if [ ! -f "$ctl_template" ]; then
  echo "Control file not found: $ctl_template" >&2
  exit 1
fi

bench_tbl_dir=$(dirname "$bench_tblbase")
if [ ! -d "$bench_tbl_dir" ]; then
  echo "Benchmark LUT directory not found: $bench_tbl_dir" >&2
  exit 1
fi

cd "$work_dir"
export LANG=C
export LC_ALL=C
export OMP_PLACES=${OMP_PLACES:-cores}
export OMP_PROC_BIND=${OMP_PROC_BIND:-close}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-$slurm_cpus_per_task}

if command -v ml >/dev/null 2>&1; then
  ml Stages/2026 GCCcore/14.3.0
  ml CMake/4.0.3
  ml ecBuild
  ml nvidia-compilers ParaStationMPI
  ml Nsight-Compute/2025.3.1
  ml Nsight-Systems/2025.5.1
fi

echo "=== Nsight Compute ==="
command -v ncu
ncu --version

export LD_LIBRARY_PATH="$repo_root/libs/build/lib:$repo_root/libs/build/lib64:${LD_LIBRARY_PATH:-}"

active_ctl="$work_dir/${case_name}.ctl"
awk -v tblbase="$bench_tblbase" '{ if ($1 == "TBLBASE") print "TBLBASE = " tblbase; else print $0; }' "$ctl_template" > "$active_ctl"

cd "$src_dir"
make clean
make -j MPI="$mpi" MPICC="$mpicc" COMPILER="$compiler_gpu" GPU=1 GPU_TARGET="$gpu_target" GPU_PIN="$gpu_pin" INFO="$info"
cd "$work_dir"

rm -rf data
mkdir -p data
"$src_dir/climatology" "$active_ctl" data/atm.tab
"$src_dir/$geometry" "$active_ctl" data/obs.tab

ncu_output="$nvidia_profile_output/formod_batch${nvidia_profile_batch}"
ncu_log="$run_dir/ncu.log"
out="/tmp/jurassic_ncu_${run_id}_b${nvidia_profile_batch}.tab"

srun -n1 -N1 ncu \
  --target-processes all \
  --metrics sm__throughput.avg.pct_of_peak_sustained_elapsed,dram__throughput.avg.pct_of_peak_sustained_elapsed,sm__warps_active.avg.pct_of_peak_sustained_active \
  --kernel-name-base function \
  --force-overwrite \
  --export "$ncu_output" \
  "$src_dir/formod" "$active_ctl" data/obs.tab data/atm.tab "$out" TASK time BATCH_SIZE "$nvidia_profile_batch" \
  2>&1 | tee "$ncu_log"

echo "Done."
echo "Nsight Systems report: ${nsys_output}.nsys-rep"
