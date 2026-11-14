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

# Create observation geometry...
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

# Check ASCII ray capacity, including a comment after the final ray...
nr_limit=$(awk '$1 == "#define" && $2 == "NR" { print $3; exit }' "$jurassic/jurassic.h")
awk -v n="$nr_limit" '!/^#/ && NF { for (i = 0; i < n; i++) print; exit }' \
    data.ref/obs.tab > data/obs_boundary_max.tab
echo '# trailing comment' >> data/obs_boundary_max.tab
$jurassic/obsfmt obs.ctl data/obs_boundary_max.tab 1 data/obs_boundary_out.tab 1 \
    > data/obs_boundary_max.log 2>&1
awk -v n="$nr_limit" '!/^#/ && NF { count++ } END { exit(count != n) }' \
    data/obs_boundary_out.tab
cp data/obs_boundary_max.tab data/obs_boundary_over.tab
awk '!/^#/ && NF { print; exit }' data.ref/obs.tab >> data/obs_boundary_over.tab
if $jurassic/obsfmt obs.ctl data/obs_boundary_over.tab 1 data/obs_boundary_rejected.tab 1 \
    > data/obs_boundary_over.log 2>&1; then
    echo 'Oversized observation file was accepted' >&2
    exit 1
fi
grep -q 'Too many rays!' data/obs_boundary_over.log
if grep -Eq 'AddressSanitizer|runtime error:' data/obs_boundary_over.log; then
    echo 'Memory error in oversized observation file' >&2
    exit 1
fi



# Compare files...
echo -e "\nCompare results..."
error=0
for f in $(ls data.ref/*.tab) ; do
    diff -q -s data/"$(basename "$f")" "$f" || error=1
done
diff -q -s data/obs_append_p0.tab data/obs.tab || error=1
diff -q -s data/obs_append_p1.tab data/obs_alt.tab || error=1
exit $error
