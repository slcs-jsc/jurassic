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

# Nontrivial species name aliases (F11, F12, F22, H2CO, N2O5, O3,
# HOCl, ClONO2)...
$jurassic/climatology alias.ctl data/alias.tab \
    T0 0 T1 0 DT 1 Z0 0 Z1 90 DZ 10 LAT 45 CLIM_FILE $clim
$jurassic/climatology alias.ctl data/alias_legacy.tab \
    T0 0 T1 0 DT 1 Z0 0 Z1 90 DZ 10 LAT 45 CLIM_FILE -

# Transition profile for the alias species, with built-in reference...
$jurassic/climatology alias.ctl data/alias_trans_jul.tab \
    T0 15724800 T1 15724800 DT 1 Z0 70 Z1 100 DZ 1 LAT 45 CLIM_FILE $clim
$jurassic/climatology alias.ctl data/alias_legacy_jul.tab \
    T0 15724800 T1 15724800 DT 1 Z0 70 Z1 100 DZ 1 LAT 45 CLIM_FILE -

# Non-negativity of the final climatology output for all available
# species across the four emitter groups (built-in climatology)...
$jurassic/climatology jumacs.ctl data/nn_jumacs.tab \
    T0 0 T1 0 DT 1 Z0 0 Z1 160 DZ 1 LAT 45 CLIM_FILE -
$jurassic/climatology alias.ctl data/nn_alias.tab \
    T0 0 T1 0 DT 1 Z0 0 Z1 160 DZ 1 LAT 45 CLIM_FILE -
$jurassic/climatology tracers1.ctl data/nn_tracers1.tab \
    T0 0 T1 0 DT 1 Z0 0 Z1 160 DZ 1 LAT 45 CLIM_FILE -
$jurassic/climatology tracers2.ctl data/nn_tracers2.tab \
    T0 0 T1 0 DT 1 Z0 0 Z1 160 DZ 1 LAT 45 CLIM_FILE -

# Same invariant for the JuMACS path (WACCM-X, valid up to 150 km),
# exercising both file values and built-in fallback species...
$jurassic/climatology jumacs.ctl data/nn_jumacs_waccmx.tab \
    T0 15724800 T1 15724800 DT 1 Z0 0 Z1 150 DZ 1 LAT 75 \
    CLIM_FILE ../../data/jumacs_waccmx_1985-2014_5deg_1km.nc
$jurassic/climatology tracers1.ctl data/nn_tracers1_waccmx.tab \
    T0 15724800 T1 15724800 DT 1 Z0 0 Z1 150 DZ 1 LAT 75 \
    CLIM_FILE ../../data/jumacs_waccmx_1985-2014_5deg_1km.nc

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
# values for all fields except OCS ($14), which instead keeps its JuMACS
# top-level value...
diff -q -s <(awk '$2>=96 {print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13}' \
    data/trans_jul.tab) \
    <(awk '$2>=96 {print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13}' \
    data/legacy_jul.tab) || error=1

# Constant extrapolation above the JuMACS top: F22 at 77, 90, 96, and
# 100 km must equal its value at 76 km and must not revert to the built-in
# climatology above 96 km, while temperature ($3) and O3 ($6) must recover
# the built-in values...
if ! awk '$2 == 76 {f = $9}
         $2 == 77 || $2 == 90 || $2 == 96 || $2 == 100 {
           if ($9 != f) exit 1
         }' data/alias_trans_jul.tab ; then
  echo "F22 not held constant above the JuMACS top!"
  error=1
fi
if diff -q <(awk '$2>=96 {print $9}' data/alias_trans_jul.tab) \
   <(awk '$2>=96 {print $9}' data/alias_legacy_jul.tab) > /dev/null ; then
  echo "F22 unexpectedly reverts to built-in climatology above 96 km!"
  error=1
fi
if ! awk '$2 == 76 {f = $11}
         $2 == 77 || $2 == 90 || $2 == 96 || $2 == 100 {
           if ($11 != f) exit 1
         }' data/alias_trans_jul.tab ; then
  echo "N2O5 not held constant above the JuMACS top!"
  error=1
fi
if diff -q <(awk '$2>=96 {print $11}' data/alias_trans_jul.tab) \
   <(awk '$2>=96 {print $11}' data/alias_legacy_jul.tab) > /dev/null ; then
  echo "N2O5 unexpectedly reverts to built-in climatology above 96 km!"
  error=1
fi
diff -q -s <(awk '$2>=96 {print $3,$6}' data/alias_trans_jul.tab) \
    <(awk '$2>=96 {print $3,$6}' data/alias_legacy_jul.tab) || error=1

# HOCl ($13) uses the normal transition again: it must converge exactly
# to the built-in climatology at and above 96 km, and must not be held
# constant above the JuMACS top...
if ! diff -q <(awk '$2>=96 {print $13}' data/alias_trans_jul.tab) \
   <(awk '$2>=96 {print $13}' data/alias_legacy_jul.tab) > /dev/null ; then
  echo "HOCl does not converge to built-in climatology above 96 km!"
  error=1
fi
if awk '$2 == 76 {f = $13}
        $2 == 96 { exit !(f == $13) }' data/alias_trans_jul.tab ; then
  echo "HOCl unexpectedly held constant above the JuMACS top!"
  error=1
fi

# ClONO2 ($14) uses the logarithmic blend: it must equal the JuMACS
# value at 76 km and the built-in climatology at and above 96 km;
# within the transition region it must follow the geometric
# interpolation between both, remain finite and non-negative, and must
# not overshoot the larger of the two endpoints (no bulge)...
if ! awk '$2 == 76 {f = $14}
         END { exit !(f > 0) }' data/alias_trans_jul.tab ; then
  echo "ClONO2 top value not positive!"
  error=1
fi
if ! diff -q <(awk '$2>=96 {print $14}' data/alias_trans_jul.tab) \
   <(awk '$2>=96 {print $14}' data/alias_legacy_jul.tab) > /dev/null ; then
  echo "ClONO2 does not converge to built-in climatology above 96 km!"
  error=1
fi
if ! awk 'NR == FNR {
           if ($2 == 76) { q = $14 }
           t[$2] = $14
           next
         }
         {
           z = $2
           if (z > 76 && z < 96) {
             x = (z - 76) / 20
             w = 0.5 * (1 + cos(atan2(0, -1) * x))
             if ($14 <= 0 || q <= 0)
               e = w * q + (1 - w) * $14
             else
               e = exp(w * log(q) + (1 - w) * log($14))
             d = t[z] - e
             if (d > 1e-3 * e || d < -1e-3 * e) exit 1
             m = (q > $14 ? q : $14)
             if (t[z] > m * (1 + 1e-3)) exit 1
           }
         }' data/alias_trans_jul.tab data/alias_legacy_jul.tab ; then
  echo "ClONO2 transition not a bulge-free log-space blend!"
  error=1
fi

# Alias transition values must be finite and non-negative...
if ! awk '$2 >= 70 && $2 <= 100 && NF {
       for (i = 6; i <= 14; i++)
         if ($i !~ /^[+-]?([0-9]+\.?[0-9]*|\.?[0-9]+)([eE][+-]?[0-9]+)?$/ \
             || $i + 0 < 0)
           exit 1
     }' data/alias_trans_jul.tab ; then
  echo "Non-finite or negative alias transition values!"
  error=1
fi

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

# The final climatology output must be finite and non-negative for every
# species at every altitude (mixing ratios are columns 7-14, all groups
# use NG = 8)...
for f in nn_jumacs nn_alias nn_tracers1 nn_tracers2 \
         nn_jumacs_waccmx nn_tracers1_waccmx ; do
  if ! awk 'NF > 13 {
         for (i = 7; i <= 14; i++)
           if ($i !~ /^[+-]?([0-9]+\.?[0-9]*|\.?[0-9]+)([eE][+-]?[0-9]+)?$/ \
               || $i + 0 < 0)
             exit 1
       }' data/$f.tab ; then
    echo "Non-finite or negative mixing ratios in final climatology ($f)!"
    error=1
  fi
done

exit $error
