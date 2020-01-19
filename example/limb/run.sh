#! /bin/bash

# Setup...
jurassic=../../src

# Create atmospheric data file...
$jurassic/climatology limb.ctl atm.tab

# Create observation geomtry...
$jurassic/limb limb.ctl obs.tab

# Call forward model...
$jurassic/formod limb.ctl obs.tab atm.tab rad.tab TASK time

# Plot results...
gnuplot <<EOF
set term png enh truecolor font "Helvetica,28" size 1600,1200 crop lw 2
set out "plot.png"
set xla "radiance [nW/(cm^2 sr cm^{-1})]"
set yla "tangent height [km]"
set mxtics
set mytics
set log x
set key spac 1.5
set key left bot
plot "rad.org" u (\$11*1e5):8 w lp pt 1 t "ref (792 cm^{-1})", \
     "rad.org" u (\$12*1e5):8 w lp pt 1 t "ref (832 cm^{-1})", \
     "rad.tab" u (\$11*1e5):8 w lp pt 2 t "test (792 cm^{-1})", \
     "rad.tab" u (\$12*1e5):8 w lp pt 2 t "test (832 cm^{-1})"
EOF

# Get differences...
echo -e "\nCheck for differences..."
diff -sq rad.tab rad.org
