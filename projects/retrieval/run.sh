#!/usr/bin/env bash
set -euo pipefail

# Run from the example directory with the bundled libraries and one CPU thread.
cd "$(dirname "$0")"
export LD_LIBRARY_PATH="../../libs/build/lib:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export LANG=C LC_ALL=C

# Check executable/table prerequisites and preserve any existing result.
for tool in climatology limb formod retrieval; do
    test -x "../../src/$tool" || { echo "Build src/$tool first" >&2; exit 1; }
done
test -f ../../tab/tria_1cm/nc_1e-6/tria_HNO3.nc || {
    echo "Missing TRIA tables; see README.md" >&2; exit 1;
}
if test -e work; then
    echo "work/ exists; move it before rerunning" >&2
    exit 1
fi
mkdir -p work/case
printf '%s\n' 'work/case' > work/dirlist.txt
start=$SECONDS

# Generate the climatological a priori.
../../src/climatology retrieval.ctl work/case/atm_apr.tab Z0 1 Z1 71 DZ 2 > work/climatology.log 2>&1

# Apply deterministic Gaussian anomalies to make the synthetic truth.
python3 perturb.py work/case/atm_apr.tab work/case/atm_truth.tab

# Generate the 16-ray limb geometry.
../../src/limb retrieval.ctl work/case/obs.tab OBSZ 780 Z0 15 Z1 45 DZ 2 > work/limb.log 2>&1

# Calculate noise-free synthetic observations from truth with EGA.
../../src/formod retrieval.ctl work/case/obs.tab work/case/atm_truth.tab work/case/obs_meas.tab > work/truth_formod.log 2>&1

# Calculate initial radiances from the unperturbed a priori.
../../src/formod retrieval.ctl work/case/obs.tab work/case/atm_apr.tab work/case/obs_initial.tab > work/prior_formod.log 2>&1

# Run optimal estimation from the a priori state.
../../src/retrieval retrieval.ctl work/dirlist.txt > work/retrieval.log 2>&1

# Record elapsed time through retrieval, then analyze and plot its outputs.
printf '%s\n' "$((SECONDS - start))" > work/total_runtime_s.txt
python3 analyze.py
