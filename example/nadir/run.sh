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

set out "plot_rad.png"
set xla "Latitude [deg]"
set yla "Brightness temperature [K]"
set mxtics
set mytics
set key spac 1.5
set key outside
plot "rad.org" u 10:11 w lp pt 1 t "REF (667.8 cm^{-1})", \
     "rad.tab" u 10:11 w lp pt 2 t "TEST (667.8 cm^{-1})", \
     "rad.org" u 10:12 w lp pt 1 t "REF (668.5 cm^{-1})", \
     "rad.tab" u 10:12 w lp pt 2 t "TEST (668.5 cm^{-1})", \
     "rad.org" u 10:13 w lp pt 1 t "REF (669.8 cm^{-1})", \
     "rad.tab" u 10:13 w lp pt 2 t "TEST (669.8 cm^{-1})"

set out "plot_tau.png"
set yla "Transmittance"
plot "rad.org" u 10:14 w lp pt 1 t "REF (667.8 cm^{-1})", \
     "rad.tab" u 10:14 w lp pt 2 t "TEST (667.8 cm^{-1})", \
     "rad.org" u 10:15 w lp pt 1 t "REF (668.5 cm^{-1})", \
     "rad.tab" u 10:15 w lp pt 2 t "TEST (668.5 cm^{-1})", \
     "rad.org" u 10:16 w lp pt 1 t "REF (669.8 cm^{-1})", \
     "rad.tab" u 10:16 w lp pt 2 t "TEST (669.8 cm^{-1})"

set out "plot_rad_diff_abs.png"
set yla "Brightness temperature difference [K]"
set yra [*:*]
set key right top
plot "< paste rad.tab rad.org" u 10:(\$27-\$11) w lp pt 1 t "TEST - REF (667.8 cm^{-1})", \
     "< paste rad.tab rad.org" u 10:(\$28-\$12) w lp pt 1 t "TEST - REF (668.5 cm^{-1})", \
     "< paste rad.tab rad.org" u 10:(\$29-\$13) w lp pt 1 t "TEST - REF (669.8 cm^{-1})"

set out "plot_rad_diff_rel.png"
set yla "Brightness temperature difference [%]"
set yra [-5:5]
plot "< paste rad.tab rad.org" u 10:(100.*(\$27-\$11)/\$11) w lp pt 1 t "TEST - REF (667.8 cm^{-1})", \
     "< paste rad.tab rad.org" u 10:(100.*(\$28-\$12)/\$12) w lp pt 1 t "TEST - REF (668.5 cm^{-1})", \
     "< paste rad.tab rad.org" u 10:(100.*(\$29-\$13)/\$13) w lp pt 1 t "TEST - REF (669.8 cm^{-1})"

set out "plot_tau_diff_abs.png"
set yla "Transmittance difference"
set yra [*:*]
plot "< paste rad.tab rad.org" u 10:(\$30-\$14) w lp pt 1 t "TEST - REF (667.8 cm^{-1})", \
     "< paste rad.tab rad.org" u 10:(\$31-\$15) w lp pt 1 t "TEST - REF (668.5 cm^{-1})", \
     "< paste rad.tab rad.org" u 10:(\$32-\$16) w lp pt 1 t "TEST - REF (669.8 cm^{-1})"

set out "plot_tau_diff_rel.png"
set yla "Transmittance difference [%]"
set yra [-5:5]
plot "< paste rad.tab rad.org" u 10:(100.*(\$30-\$14)/\$14) w lp pt 1 t "TEST - REF (667.8 cm^{-1})", \
     "< paste rad.tab rad.org" u 10:(100.*(\$31-\$15)/\$15) w lp pt 1 t "TEST - REF (668.5 cm^{-1})", \
     "< paste rad.tab rad.org" u 10:(100.*(\$32-\$16)/\$16) w lp pt 1 t "TEST - REF (669.8 cm^{-1})"
EOF

# Get differences...
echo -e "\nCheck for differences..."
diff -sq rad.tab rad.org
