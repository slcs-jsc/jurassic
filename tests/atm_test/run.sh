#! /bin/bash

# Set environment...
export LD_LIBRARY_PATH=../../libs/build/lib:${LD_LIBRARY_PATH}
export OMP_NUM_THREADS=4
export LANG=C
export LC_ALL=C

set -euo pipefail

# Setup...
jurassic=../../src

# Create directory...
rm -rf data && mkdir -p data

# Create atmospheric data file...
$jurassic/climatology atm.ctl data/atm.tab

# Test interpolation...
$jurassic/climatology atm.ctl data/atm_3km.tab DZ 3
$jurassic/interpolate atm.ctl data/atm_3km.tab data/atm.tab data/atm_intpol.tab

# Test hydrostatic build up...
$jurassic/hydrostatic atm.ctl data/atm.tab data/atm_hyd.tab HYDZ 0.0

# Test binary file-I/O...
$jurassic/atmfmt atm.ctl data/atm.tab 1 data/atm_asc2bin.bin 2
$jurassic/atmfmt atm.ctl data/atm_asc2bin.bin 2 data/atm_bin2asc.tab 1

# Test netCDF file-I/O...
$jurassic/atmfmt atm.ctl data/atm.tab 1 data/atm_asc2nc.nc 3
$jurassic/atmfmt atm.ctl data/atm_asc2nc.nc 3 data/atm_nc2asc.tab 1

# Compare files...
echo -e "\nCompare results..."
error=0
for f in $(ls data.ref/*.tab) ; do
    diff -q -s data/"$(basename "$f")" "$f" || error=1
done
exit $error
