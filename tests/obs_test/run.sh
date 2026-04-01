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
$jurassic/climatology obs.ctl data/atm.tab

# Create observation geomtry...
$jurassic/zenith obs.ctl data/obs.tab
$jurassic/zenith obs.ctl data/obs_alt.tab OBSZ 9.0

# Test raytracer...
$jurassic/raytrace obs.ctl data/obs.tab data/atm.tab data/raytrace.tab LOSBASE data/los

# Convert to spectrum...
$jurassic/obs2spec obs.ctl data/obs.tab data/spec.tab

# Test binary file-I/O...
$jurassic/obsfmt obs.ctl data/obs.tab 1 data/obs_asc2bin.bin 2
$jurassic/obsfmt obs.ctl data/obs_asc2bin.bin 2 data/obs_bin2asc.tab 1

# Test netCDF file-I/O...
$jurassic/obsfmt obs.ctl data/obs.tab 1 data/obs_asc2nc.nc 3
$jurassic/obsfmt obs.ctl data/obs_asc2nc.nc 3 data/obs_nc2asc.tab 1

# Test netCDF file-I/O with appending...
$jurassic/obsfmt obs.ctl data/obs.tab 1 data/obs_append.nc 3
$jurassic/obsfmt obs.ctl data/obs_alt.tab 1 data/obs_append.nc 3 PROF_OUT 1
$jurassic/obsfmt obs.ctl data/obs_append.nc 3 data/obs_append_p0.tab 1
$jurassic/obsfmt obs.ctl data/obs_append.nc 3 data/obs_append_p1.tab 1 PROF_IN 1




# Compare files...
echo -e "\nCompare results..."
error=0
for f in $(ls data.ref/*.tab) ; do
    diff -q -s data/"$(basename "$f")" "$f" || error=1
done
diff -q -s data/obs_append_p0.tab data/obs.tab || error=1
diff -q -s data/obs_append_p1.tab data/obs_alt.tab || error=1
exit $error
