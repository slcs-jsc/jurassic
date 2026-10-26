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

# December branch (past the December center) and negative time...
$jurassic/climatology jumacs.ctl data/dec.tab \
    T0 30758400 T1 30758400 DT 1 Z0 0 Z1 90 DZ 10 LAT 45 CLIM_FILE $clim
$jurassic/climatology jumacs.ctl data/neg.tab \
    T0 -600 T1 -600 DT 1 Z0 0 Z1 90 DZ 10 LAT 45 CLIM_FILE $clim

# Same phases as the December and negative-time tests one year later...
$jurassic/climatology jumacs.ctl data/dec2.tab \
    T0 62316000 T1 62316000 DT 1 Z0 0 Z1 90 DZ 10 LAT 45 CLIM_FILE $clim
$jurassic/climatology jumacs.ctl data/neg2.tab \
    T0 31557000 T1 31557000 DT 1 Z0 0 Z1 90 DZ 10 LAT 45 CLIM_FILE $clim

# Built-in climatology for the special and unmapped species...
$jurassic/climatology jumacs.ctl data/legacy.tab \
    T0 0 T1 0 DT 1 Z0 0 Z1 90 DZ 10 LAT 45 CLIM_FILE -

# Polar clamp: LAT 89 must behave like the grid edge LAT 87.5...
$jurassic/climatology jumacs.ctl data/polar89.tab \
    T0 0 T1 0 DT 1 Z0 0 Z1 90 DZ 10 LAT 89 CLIM_FILE $clim
$jurassic/climatology jumacs.ctl data/polar875.tab \
    T0 0 T1 0 DT 1 Z0 0 Z1 90 DZ 10 LAT 87.5 CLIM_FILE $clim

# Same position in the climatological year across common and leap years...
$jurassic/climatology jumacs.ctl data/ny2010.tab \
    T0 315619200 T1 315619200 DT 1 Z0 0 Z1 75 DZ 15 LAT 45 CLIM_FILE $clim
$jurassic/climatology jumacs.ctl data/ny2012.tab \
    T0 378734400 T1 378734400 DT 1 Z0 0 Z1 75 DZ 15 LAT 45 CLIM_FILE $clim
$jurassic/climatology jumacs.ctl data/mar2010.tab \
    T0 320716800 T1 320716800 DT 1 Z0 0 Z1 75 DZ 15 LAT 45 CLIM_FILE $clim
$jurassic/climatology jumacs.ctl data/mar2012.tab \
    T0 383832000 T1 383832000 DT 1 Z0 0 Z1 75 DZ 15 LAT 45 CLIM_FILE $clim

# Same climatological phase in the very first year...
$jurassic/climatology jumacs.ctl data/phase.tab \
    T0 43200 T1 43200 DT 1 Z0 0 Z1 75 DZ 15 LAT 45 CLIM_FILE $clim

# Nontrivial species name aliases (F11, F12, F22, H2CO, O3)...
$jurassic/climatology alias.ctl data/alias.tab \
    T0 0 T1 0 DT 1 Z0 0 Z1 90 DZ 10 LAT 45 CLIM_FILE $clim
$jurassic/climatology alias.ctl data/alias_legacy.tab \
    T0 0 T1 0 DT 1 Z0 0 Z1 90 DZ 10 LAT 45 CLIM_FILE -

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

# December and negative-time invariances across one climatological year...
diff -q -s <(awk '{print $2,$3,$4,$5,$6,$7,$9,$10,$11,$12,$13,$14}' \
    data/dec.tab) \
    <(awk '{print $2,$3,$4,$5,$6,$7,$9,$10,$11,$12,$13,$14}' \
    data/dec2.tab) || error=1
diff -q -s <(awk '{print $2,$3,$4,$5,$6,$7,$9,$10,$11,$12,$13,$14}' \
    data/neg.tab) \
    <(awk '{print $2,$3,$4,$5,$6,$7,$9,$10,$11,$12,$13,$14}' \
    data/neg2.tab) || error=1

# Same climatological phase in the first year (pre-January cyclic branch)...
diff -q -s <(awk '{print $2,$3,$4,$5,$6,$7,$9,$10,$11,$12,$13,$14}' \
    data/phase.tab) \
    <(awk '{print $2,$3,$4,$5,$6,$7,$9,$10,$11,$12,$13,$14}' \
    data/ny2010.tab) || error=1

# Special/unmapped species (CO2, O2, N2, SF6) must equal built-in values...
diff -q -s <(awk '{print $2,$8,$11,$12,$13}' data/jan45.tab) \
    <(awk '{print $2,$8,$11,$12,$13}' data/legacy.tab) || error=1

# Altitudes above the JuMACS top (z > 76 km) must equal built-in values...
diff -q -s <(awk '$2>76' data/jan45.tab) \
    <(awk '$2>76' data/legacy.tab) || error=1

# Mapped alias species must differ from their built-in values...
if diff -q <(awk '!/^#/{print $6,$7,$8,$9}' data/alias.tab) \
   <(awk '!/^#/{print $6,$7,$8,$9}' data/alias_legacy.tab) > /dev/null ; then
  echo "Alias species columns unexpectedly equal built-in values!"
  error=1
fi

exit $error
