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

# Exercise ray entry and exit at the atmospheric boundaries.
printf '0 0 0 0 1000 280 0\n0 100 0 0 1 220 0\n' > data/ray_boundary_atm.tab
for case in inside bottom top above_enter; do
    case "$case" in
        inside) obsz=50; vpz=100 ;;
        bottom) obsz=0; vpz=10 ;;
        top) obsz=100; vpz=150 ;;
        above_enter) obsz=150; vpz=0 ;;
    esac
    printf '0 %s 0 0 %s 0 0 0 0 0\n' "$obsz" "$vpz" > "data/ray_${case}_obs.tab"
    "$jurassic/raytrace" - "data/ray_${case}_obs.tab" data/ray_boundary_atm.tab \
        "data/ray_${case}.tab" LOSBASE "data/los_${case}" REFRAC 0 \
        RAYDS 10 RAYDZ 0 > "data/ray_${case}.log" 2>&1
    awk '!/^#/ && NF { found=1 } END { exit(!found) }' "data/los_${case}.0.tab"
done

# Rays starting below the surface or missing the atmosphere must fail.
for case in below above_miss; do
    if [ "$case" = below ]; then obsz=-1; vpz=10; else obsz=150; vpz=200; fi
    printf '0 %s 0 0 %s 0 0 0 0 0\n' "$obsz" "$vpz" > "data/ray_${case}_obs.tab"
    if "$jurassic/raytrace" - "data/ray_${case}_obs.tab" data/ray_boundary_atm.tab \
        "data/ray_${case}.tab" LOSBASE "data/los_${case}" REFRAC 0 \
        RAYDS 10 RAYDZ 0 > "data/ray_${case}.log" 2>&1; then
        echo "Invalid ray geometry was accepted: $case" >&2
        exit 1
    fi
    grep -q 'Ray tracing failed!' "data/ray_${case}.log"
done

# Compare files...
echo -e "\nCompare results..."
error=0
for f in $(ls data.ref/*.tab) ; do
    diff -q -s data/"$(basename "$f")" "$f" || error=1
done
diff -q -s data/obs_append_p0.tab data/obs.tab || error=1
diff -q -s data/obs_append_p1.tab data/obs_alt.tab || error=1
exit $error
