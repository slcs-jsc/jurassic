#!/bin/bash
#SBATCH --account=slmet
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --time=09:00:00
#SBATCH --array=0-2
#SBATCH --exclusive
#SBATCH --job-name=e4_scaling_axes

# OpenMP scaling over geometry, channels and gas sets (see scaling_axes.sh),
# one array task per case; zenith takes longest.
#   sbatch run_scaling_axes_juwels.sh
#   MODES=t1check AXES=geometry sbatch --time=01:00:00 run_scaling_axes_juwels.sh
#   MODES=batches AXES=geometry sbatch --time=01:00:00 run_scaling_axes_juwels.sh
#   python3 eval_scaling_axes.py runs/scaling_axes_<array job id>

if [ -n "${SLURM_SUBMIT_DIR:-}" ] && [ -f "$SLURM_SUBMIT_DIR/base.sh" ]; then
  jr_scripts_dir="$SLURM_SUBMIT_DIR"
else
  jr_scripts_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi

source "$jr_scripts_dir/scaling_axes.sh"
