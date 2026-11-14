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

# Test netCDF file-I/O with appending...
$jurassic/atmfmt atm.ctl data/atm.tab 1 data/atm_append.nc 3
$jurassic/atmfmt atm.ctl data/atm_hyd.tab 1 data/atm_append.nc 3 PROF_OUT 1
$jurassic/atmfmt atm.ctl data/atm_append.nc 3 data/atm_append_p0.tab 1
$jurassic/atmfmt atm.ctl data/atm_append.nc 3 data/atm_append_p1.tab 1 PROF_IN 1

# Check ASCII profile capacity, including a comment after the final point...
np_limit=$(awk '$1 == "#define" && $2 == "NP" { print $3; exit }' "$jurassic/jurassic.h")
awk -v n="$np_limit" '!/^#/ && NF { for (i = 0; i < n; i++) print; exit }' \
    data.ref/atm.tab > data/atm_boundary_max.tab
echo '# trailing comment' >> data/atm_boundary_max.tab
$jurassic/atmfmt atm.ctl data/atm_boundary_max.tab 1 data/atm_boundary_out.tab 1 \
    > data/atm_boundary_max.log 2>&1
awk -v n="$np_limit" '!/^#/ && NF { count++ } END { exit(count != n) }' \
    data/atm_boundary_out.tab
cp data/atm_boundary_max.tab data/atm_boundary_over.tab
awk '!/^#/ && NF { print; exit }' data.ref/atm.tab >> data/atm_boundary_over.tab
if $jurassic/atmfmt atm.ctl data/atm_boundary_over.tab 1 data/atm_boundary_rejected.tab 1 \
    > data/atm_boundary_over.log 2>&1; then
    echo 'Oversized atmospheric profile was accepted' >&2
    exit 1
fi
grep -q 'Too many data points!' data/atm_boundary_over.log

# Compare files...
echo -e "\nCompare results..."
error=0
for f in $(ls data.ref/*.tab) ; do
    diff -q -s data/"$(basename "$f")" "$f" || error=1
done
diff -q -s data/atm_append_p0.tab data/atm.tab || error=1
diff -q -s data/atm_append_p1.tab data/atm_hyd.tab || error=1
exit $error
