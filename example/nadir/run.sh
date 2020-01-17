#! /bin/bash

# Setup...
jurassic=../../src

# Create atmospheric data file...
$jurassic/climatology nadir.ctl atm.tab

# Create observation geomtry...
$jurassic/nadir nadir.ctl obs.tab

# Call forward model...
$jurassic/formod nadir.ctl obs.tab atm.tab rad.tab TASK time

# Plot results...
gnuplot <<EOF
set term png enh truecolor font "Helvetica,28" size 1600,1200 crop lw 2
set out "plot.png"
set xla "latitude [deg]"
set yla "brightness temperature [K]"
set mxtics
set mytics
set key spac 1.5
set key outside
plot "rad.org" u 10:11 w lp pt 1 t "ref (667.8 cm^{-1})", \
     "rad.tab" u 10:11 w lp pt 2 t "test (667.8 cm^{-1})", \
     "rad.org" u 10:12 w lp pt 1 t "ref (668.5 cm^{-1})", \
     "rad.tab" u 10:12 w lp pt 2 t "test (668.5 cm^{-1})", \
     "rad.org" u 10:13 w lp pt 1 t "ref (669.8 cm^{-1})", \
     "rad.tab" u 10:13 w lp pt 2 t "test (669.8 cm^{-1})"
EOF

# Get differences...
echo -e "\nCheck for differences..."
diff -sq rad.tab rad.org
