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

# Create filter functions...
nu=667.5000
$jurassic/filter - data/boxcar_$nu.filt FILTER_TYPE 0 FILTER_CENTER $nu
$jurassic/filter - data/triangle_$nu.filt FILTER_TYPE 1 FILTER_CENTER $nu
$jurassic/filter - data/gaussian_$nu.filt FILTER_TYPE 2 FILTER_CENTER $nu
$jurassic/filter - data/sinc_$nu.filt FILTER_TYPE 3 FILTER_CENTER $nu
$jurassic/filter - data/norton_beer_$nu.filt FILTER_TYPE 4 FILTER_CENTER $nu

# Create table files...
for filter in boxcar triangle gaussian sinc norton_beer ; do
    $jurassic/tblgen 1000 285 opt_01000.asc data/${filter}_$nu.filt > data/${filter}_${nu}_CO2.tab
done

# Convert an existing two-channel table set to netCDF and back...
filter=boxcar
ctl="ND 2 NU[0] 792.0000 NU[1] 832.0000 NG 1 EMITTER[0] CO2"
$jurassic/tblfmt - ../data/${filter} 1 data/${filter}_multi_asc2bin 2 $ctl
$jurassic/tblfmt - ../data/${filter} 1 data/${filter}_multi_asc2nc 3 $ctl
$jurassic/tblfmt - data/${filter}_multi_asc2bin 2 data/${filter}_multi_bin2asc 1 $ctl
$jurassic/tblfmt - data/${filter}_multi_asc2nc 3 data/${filter}_multi_nc2asc 1 $ctl

# Compress table...
ctl="ND 1 NU[0] $nu NG 1 EMITTER[0] CO2"
$jurassic/tblrdc - data/${filter} 1 data/${filter}_rdc 1 $ctl

# Compare files...
echo -e "\nCompare results..."
error=0
for f in data.ref/boxcar_667.5000.filt \
         data.ref/boxcar_rdc_667.5000.filt \
         data.ref/gaussian_667.5000.filt \
         data.ref/norton_beer_667.5000.filt \
         data.ref/sinc_667.5000.filt \
         data.ref/triangle_667.5000.filt \
         data.ref/boxcar_667.5000_CO2.tab \
         data.ref/boxcar_rdc_667.5000_CO2.tab \
         data.ref/gaussian_667.5000_CO2.tab \
         data.ref/norton_beer_667.5000_CO2.tab \
         data.ref/sinc_667.5000_CO2.tab \
         data.ref/triangle_667.5000_CO2.tab ; do
    diff -q -s data/"$(basename "$f")" "$f" || error=1
done
for f in data/boxcar_multi_bin2asc_792.0000.filt \
         data/boxcar_multi_bin2asc_792.0000_CO2.tab \
         data/boxcar_multi_bin2asc_832.0000.filt \
         data/boxcar_multi_bin2asc_832.0000_CO2.tab ; do
    base="$(basename "$f")"
    diff -q -s "data/boxcar_multi_nc2asc_${base#boxcar_multi_bin2asc_}" "$f" || error=1
done
exit $error
