#! /bin/bash

# Set environment...
export LD_LIBRARY_PATH=../../libs/build/lib:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=4
export LANG=C
export LC_ALL=C

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

# Convert table files from ASCII to binary and netCDF...
filter=boxcar
ctl="ND 1 NU[0] $nu NG 1 EMITTER[0] CO2"
$jurassic/tblfmt - data/${filter} 1 data/${filter}_asc2bin 2 $ctl
$jurassic/tblfmt - data/${filter} 1 data/${filter}_asc2nc  3 $ctl
$jurassic/tblfmt - data/${filter}_asc2bin 2 data/${filter}_bin2asc 1 $ctl
$jurassic/tblfmt - data/${filter}_asc2nc  3 data/${filter}_nc2asc  1 $ctl

# Compare files...
echo -e "\nCompare results..."
error=0
for f in data.ref/*.filt data.ref/*.tab ; do
    diff -q -s data/"$(basename "$f")" "$f" || error=1
done
exit $error
