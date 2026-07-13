#! /bin/bash

# Set environment...
export LD_LIBRARY_PATH=../../libs/build/lib:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}
export LANG=C
export LC_ALL=C
set -euo pipefail

# Setup...
jurassic=${JURASSIC:-../../src}
ctl=${CTLFILE:-zenith.ctl}

# Create directory...
rm -rf data && mkdir -p data

# Create atmospheric and observation data...
$jurassic/climatology $ctl data/atm_zenith.tab
$jurassic/zenith $ctl data/obs_zenith.tab

# Call forward model and kernel...
$jurassic/formod $ctl data/obs_zenith.tab data/atm_zenith.tab data/rad_zenith.tab TASK time
$jurassic/kernel $ctl data/obs_zenith.tab data/atm_zenith.tab data/kernel_zenith.tab

# Compare files...
echo -e "\nCompare results..."
error=0
for f in $(ls data.ref/*.tab) ; do
    diff -q -s data/"$(basename "$f")" "$f" || error=1
done
exit $error
