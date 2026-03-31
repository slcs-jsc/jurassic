#! /bin/bash

# Set environment...
export LD_LIBRARY_PATH=../../libs/build/lib:${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=4
export LANG=C
export LC_ALL=C
set -euo pipefail

# Setup...
jurassic=../../src

# Create directory...
rm -rf data && mkdir -p data
echo "data" > data/dirlist.txt

# Create atmospheric data file...
$jurassic/climatology ret.ctl data/atm_apr.tab

# Create observation geomtry...
$jurassic/limb ret.ctl data/obs.tab

# Call forward model...
$jurassic/formod ret.ctl data/obs.tab data/atm_apr.tab data/obs_meas.tab

# Retrieval...
$jurassic/retrieval ret.ctl data/dirlist.txt

# Compare files...
echo -e "\nCompare results..."
error=0
for f in data.ref/*.tab ; do
    diff -q -s data/"$(basename "$f")" "$f" || error=1
done
exit $error
