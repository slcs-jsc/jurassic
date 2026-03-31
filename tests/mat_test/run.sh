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

# Create atmospheric and observation data...
$jurassic/climatology mat.ctl data/atm.tab
$jurassic/nadir mat.ctl data/obs.tab

# Create source matrix in y/x/r ordering...
$jurassic/kernel mat.ctl data/obs.tab data/atm.tab data/matrix.tab

# Test binary file-I/O...
$jurassic/matfmt mat.ctl data/atm.tab data/obs.tab y x r \
  data/matrix.tab 1 data/matrix_asc2bin.bin 2
$jurassic/matfmt mat.ctl data/atm.tab data/obs.tab y x r \
  data/matrix_asc2bin.bin 2 data/matrix_bin2asc.tab 1

# Test netCDF file-I/O...
$jurassic/matfmt mat.ctl data/atm.tab data/obs.tab y x r \
  data/matrix.tab 1 data/matrix_asc2nc.nc 3
$jurassic/matfmt mat.ctl data/atm.tab data/obs.tab y x r \
  data/matrix_asc2nc.nc 3 data/matrix_nc2asc.tab 1

# Compare files...
echo -e "\nCompare results..."
diff -q -s data/matrix_bin2asc.tab data/matrix.tab
diff -q -s data/matrix_nc2asc.tab data/matrix.tab
