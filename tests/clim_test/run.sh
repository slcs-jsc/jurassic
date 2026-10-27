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

# Smooth transition to the built-in climatology above the JuMACS top...
$jurassic/climatology jumacs.ctl data/trans_jul.tab \
    T0 15724800 T1 15724800 DT 1 Z0 70 Z1 100 DZ 1 LAT 45 CLIM_FILE $clim
$jurassic/climatology jumacs.ctl data/legacy_jul.tab \
    T0 15724800 T1 15724800 DT 1 Z0 70 Z1 100 DZ 1 LAT 45 CLIM_FILE -

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

# At the transition inlet (70 km), values must equal pure JuMACS behavior...
diff -q -s <(awk '$2==70' data/trans_jul.tab) \
    <(awk '$2==70' data.ref/jul45.tab) || error=1

# Altitudes above the transition region (z >= 96 km) must equal built-in
# values...
diff -q -s <(awk '$2>=96' data/trans_jul.tab) \
    <(awk '$2>=96' data/legacy_jul.tab) || error=1

# Within the transition region (76 < z < 96 km), values must differ from
# pure built-in values...
if diff -q <(awk '$2>76 && $2<96' data/trans_jul.tab) \
   <(awk '$2>76 && $2<96' data/legacy_jul.tab) > /dev/null ; then
  echo "Transition band unexpectedly equals built-in values!"
  error=1
fi

# The hard switch at the JuMACS top must be removed: the temperature step
# across 76 to 77 km must be much smaller than towards the built-in
# climatology...
if ! awk -v t76="$(awk '$2 == 76 {print $6}' data/trans_jul.tab)" \
        -v t77="$(awk '$2 == 77 {print $6}' data/trans_jul.tab)" \
        -v l77="$(awk '$2 == 77 {print $6}' data/legacy_jul.tab)" \
    'BEGIN { exit !((t77 - t76)^2 < (l77 - t76)^2) }' ; then
  echo "Temperature jump at the JuMACS top not removed!"
  error=1
fi

# Transition values must be finite, and temperatures and mixing ratios
# must remain non-negative...
if ! awk '$2 >= 70 && $2 <= 100 && NF {
       for (i = 6; i <= 14; i++)
         if ($i !~ /^[+-]?([0-9]+\.?[0-9]*|\.?[0-9]+)([eE][+-]?[0-9]+)?$/ \
             || $i + 0 < 0)
           exit 1
     }' data/trans_jul.tab ; then
  echo "Non-finite or negative transition values!"
  error=1
fi

# Mapped alias species must differ from their built-in values...
if diff -q <(awk '!/^#/{print $6,$7,$8,$9}' data/alias.tab) \
   <(awk '!/^#/{print $6,$7,$8,$9}' data/alias_legacy.tab) > /dev/null ; then
  echo "Alias species columns unexpectedly equal built-in values!"
  error=1
fi

exit $error
