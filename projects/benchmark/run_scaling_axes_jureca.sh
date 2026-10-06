#!/bin/bash
#SBATCH --account=slmet
#SBATCH --partition=dc-cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=16:00:00
#SBATCH --array=0-2
#SBATCH --exclusive
#SBATCH --job-name=e4_scaling_axes

# OpenMP scaling over geometry, channels and gas sets (see scaling_axes.sh).
# One array task per case (zenith, nadir, limb), strong scaling over threads for
# every setting: about 2-14 h for the longest task (zenith).
# Cheaper: CURVE_SETTINGS="nd128 priority_full" sbatch --time=06:00:00 run_scaling_axes_jureca.sh
# With weak scaling too (MODES="strong weak") use --time=24:00:00.
#   sbatch run_scaling_axes_jureca.sh
#   python3 eval_scaling_axes.py runs/scaling_axes_<array job id>

if [ -n "${SLURM_SUBMIT_DIR:-}" ] && [ -f "$SLURM_SUBMIT_DIR/base.sh" ]; then
  jr_scripts_dir="$SLURM_SUBMIT_DIR"
else
  jr_scripts_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi

source "$jr_scripts_dir/scaling_axes.sh"
