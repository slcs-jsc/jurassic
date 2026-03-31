#! /bin/bash

# Set environment...
export LD_LIBRARY_PATH=../../libs/build/lib:${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=4
set -euo pipefail

# Setup...
trac=../../src

# Create directories...
rm -rf data && mkdir -p data

echo "checking time..."
(for year in 1900 1980 2000 2020 2100 ; do
    for mon in 1 7 12 ; do
	for day in 1 15 31 ; do
	    for hour in 0 12 24 ; do
		t0=$($trac/time2jsec $year $mon $day $hour 0 0 0)
		echo "$($trac/jsec2time "$t0") = $t0"
		d0=$($trac/day2doy $year $mon $day)
		echo "$($trac/doy2day $d0) = $d0"
	    done
	done
    done
done) > data/time.tab

echo "checking planck and brightness..."
$trac/planck 180 320 10 500 3000 100 > data/planck.tab
$trac/brightness 1e-8 0.1 0.01 500 3000 100 > data/brightness.tab

# Compare files...
echo -e "\nCompare results..."
error=0
for f in data.ref/*.tab ; do
    diff -q -s data/"$(basename "$f")" "$f" || error=1
done
exit $error
