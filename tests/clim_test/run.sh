#! /bin/bash

# Set environment...
export LD_LIBRARY_PATH=../../libs/build/lib:${LD_LIBRARY_PATH}
export OMP_NUM_THREADS=4
export LANG=C
export LC_ALL=C
set -euo pipefail

# Setup...
jurassic=../../src
clim=../../data/jumacs_emac_refd1_1985-2014_5deg_1km.nc

# Create directory...
rm -rf data && mkdir -p data

# Reference profiles using the JuMACS climatology...
$jurassic/climatology jumacs.ctl data/jan45.tab \
    T0 0 T1 0 DT 1 Z0 0 Z1 90 DZ 10 LAT 45 CLIM_FILE $clim
$jurassic/climatology jumacs.ctl data/jul45.tab \
    T0 15724800 T1 15724800 DT 1 Z0 0 Z1 90 DZ 10 LAT 45 CLIM_FILE $clim
$jurassic/climatology jumacs.ctl data/jan_south.tab \
    T0 0 T1 0 DT 1 Z0 0 Z1 90 DZ 10 LAT -45 CLIM_FILE $clim
$jurassic/climatology jumacs.ctl data/feb29.tab \
    T0 383788800 T1 383788800 DT 1 Z0 0 Z1 75 DZ 15 LAT 45 CLIM_FILE $clim

# Built-in climatology for the special and unmapped species...
$jurassic/climatology jumacs.ctl data/legacy.tab \
    T0 0 T1 0 DT 1 Z0 0 Z1 90 DZ 10 LAT 45 CLIM_FILE -

# Polar clamp: LAT 89 must behave like the grid edge LAT 87.5...
$jurassic/climatology jumacs.ctl data/polar89.tab \
    T0 0 T1 0 DT 1 Z0 0 Z1 90 DZ 10 LAT 89 CLIM_FILE $clim
$jurassic/climatology jumacs.ctl data/polar875.tab \
    T0 0 T1 0 DT 1 Z0 0 Z1 90 DZ 10 LAT 87.5 CLIM_FILE $clim

# New-Year and first-of-March invariances across leap and common years...
$jurassic/climatology jumacs.ctl data/ny2010.tab \
    T0 315619200 T1 315619200 DT 1 Z0 0 Z1 75 DZ 15 LAT 45 CLIM_FILE $clim
$jurassic/climatology jumacs.ctl data/ny2012.tab \
    T0 378691200 T1 378691200 DT 1 Z0 0 Z1 75 DZ 15 LAT 45 CLIM_FILE $clim
$jurassic/climatology jumacs.ctl data/mar2010.tab \
    T0 320716800 T1 320716800 DT 1 Z0 0 Z1 75 DZ 15 LAT 45 CLIM_FILE $clim
$jurassic/climatology jumacs.ctl data/mar2012.tab \
    T0 383875200 T1 383875200 DT 1 Z0 0 Z1 75 DZ 15 LAT 45 CLIM_FILE $clim

# Error handling: nonexistent and incompatible climatology files...
if $jurassic/climatology jumacs.ctl data/error.tab \
     T0 0 LAT 45 CLIM_FILE /nonexistent_jumacs.nc 2>/dev/null ; then
  echo "Error handling failed for nonexistent climatology file!"
  exit 1
fi
$jurassic/atmfmt jumacs.ctl data/jan45.tab 1 data/wrong.nc 3 > /dev/null
if $jurassic/climatology jumacs.ctl data/error.tab \
     T0 0 LAT 45 CLIM_FILE data/wrong.nc 2>/dev/null ; then
  echo "Error handling failed for incompatible climatology file!"
  exit 1
fi

# Compare reference files...
echo -e "\nCompare results..."
error=0
for f in $(ls data.ref/*.tab) ; do
    diff -q -s data/"$(basename "$f")" "$f" || error=1
done

# Polar clamp (all columns except time and latitude)...
diff -q -s <(awk '{print $2,$3,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14}' \
    data/polar89.tab) \
    <(awk '{print $2,$3,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14}' \
    data/polar875.tab) || error=1

# Invariances across calendar years (excluding time and CO2, which trends)...
diff -q -s <(awk '{print $2,$3,$4,$5,$6,$7,$9,$10,$11,$12,$13,$14}' \
    data/ny2010.tab) \
    <(awk '{print $2,$3,$4,$5,$6,$7,$9,$10,$11,$12,$13,$14}' \
    data/ny2012.tab) || error=1
diff -q -s <(awk '{print $2,$3,$4,$5,$6,$7,$9,$10,$11,$12,$13,$14}' \
    data/mar2010.tab) \
    <(awk '{print $2,$3,$4,$5,$6,$7,$9,$10,$11,$12,$13,$14}' \
    data/mar2012.tab) || error=1

# Special/unmapped species (CO2, O2, N2, SF6) must equal built-in values...
diff -q -s <(awk '{print $2,$8,$11,$12,$13}' data/jan45.tab) \
    <(awk '{print $2,$8,$11,$12,$13}' data/legacy.tab) || error=1

# Altitudes above the JuMACS top (z > 76 km) must equal built-in values...
diff -q -s <(awk '$2>76' data/jan45.tab) \
    <(awk '$2>76' data/legacy.tab) || error=1

exit $error
