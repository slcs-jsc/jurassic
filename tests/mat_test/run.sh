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
$jurassic/nadir mat.ctl data/obs_alt.tab OBSZ 650

# Create source matrix in y/x/r ordering...
$jurassic/kernel mat.ctl data/obs.tab data/atm.tab data/matrix.tab
$jurassic/kernel mat.ctl data/obs_alt.tab data/atm.tab data/matrix_alt.tab

# Verify batched kernel calculations through the command-line interface...
mkdir -p data/batch0 data/batch1 data/batch2
cp data/atm.tab data/batch0/atm.tab
cp data/atm.tab data/batch1/atm.tab
cp data/atm.tab data/batch2/atm.tab
cp data/obs.tab data/batch0/obs.tab
cp data/obs_alt.tab data/batch1/obs.tab
cp data/obs.tab data/batch2/obs.tab
printf "%s\n" data/batch0 data/batch1 data/batch2 > data/dirlist
$jurassic/kernel mat.ctl obs.tab atm.tab matrix.tab DIRLIST data/dirlist KERNEL_BATCH 2
diff -q -s data/batch0/matrix.tab data/matrix.tab
diff -q -s data/batch1/matrix.tab data/matrix_alt.tab
diff -q -s data/batch2/matrix.tab data/matrix.tab

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

# Test netCDF file-I/O with appending...
$jurassic/matfmt mat.ctl data/atm.tab data/obs.tab y x r \
  data/matrix.tab 1 data/matrix_append.nc 3
$jurassic/matfmt mat.ctl data/atm.tab data/obs_alt.tab y x r \
  data/matrix_alt.tab 1 data/matrix_append.nc 3 PROF_OUT 1
$jurassic/matfmt mat.ctl data/atm.tab data/obs.tab y x r \
  data/matrix_append.nc 3 data/matrix_append_p0.tab 1
$jurassic/matfmt mat.ctl data/atm.tab data/obs_alt.tab y x r \
  data/matrix_append.nc 3 data/matrix_append_p1.tab 1 PROF_IN 1

# Compare files...
echo -e "\nCompare results..."
diff -q -s data/matrix_bin2asc.tab data/matrix.tab
diff -q -s data/matrix_nc2asc.tab data/matrix.tab
diff -q -s data/matrix_append_p0.tab data/matrix.tab
diff -q -s data/matrix_append_p1.tab data/matrix_alt.tab
