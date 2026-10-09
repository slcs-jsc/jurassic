#!/bin/bash
#SBATCH --account=slmet
#SBATCH --partition=dc-cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=10:00:00
#SBATCH --array=0-2
#SBATCH --exclusive
#SBATCH --job-name=e4_scaling_axes

# OpenMP scaling over geometry, channels and gas sets (see scaling_axes.sh),
# one array task per case. Measure the cost per setting first:
#   MODES=cost REPS=1 sbatch --time=02:00:00 run_scaling_axes_jureca.sh
#   sbatch run_scaling_axes_jureca.sh
#   MODES=t1check AXES=geometry REPS=1 sbatch --time=02:00:00 run_scaling_axes_jureca.sh
#   MODES=batches AXES=geometry sbatch --time=02:00:00 run_scaling_axes_jureca.sh
#   python3 eval_scaling_axes.py runs/scaling_axes_<array job id>

if [ -n "${SLURM_SUBMIT_DIR:-}" ] && [ -f "$SLURM_SUBMIT_DIR/base.sh" ]; then
  jr_scripts_dir="$SLURM_SUBMIT_DIR"
else
  jr_scripts_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi

source "$jr_scripts_dir/scaling_axes.sh"
