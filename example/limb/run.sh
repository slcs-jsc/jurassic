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

set out "plot_rad.png"
set xla "Radiance [nW / (cm^2 sr cm^{-1})]"
set yla "Tangent height [km]"
set mxtics
set mytics
set log x
set key spac 1.5
set key left bot
plot "rad.org" u (\$11*1e5):8 w lp pt 1 t "REF (792 cm^{-1})", \
     "rad.org" u (\$12*1e5):8 w lp pt 1 t "REF (832 cm^{-1})", \
     "rad.tab" u (\$11*1e5):8 w lp pt 2 t "TEST (792 cm^{-1})", \
     "rad.tab" u (\$12*1e5):8 w lp pt 2 t "TEST (832 cm^{-1})"

set out "plot_tau.png"
set xla "Transmittance"
unset log x
set key left top
plot "rad.org" u (\$13*1.):8 w lp pt 1 t "REF (792 cm^{-1})", \
     "rad.org" u (\$14*1.):8 w lp pt 1 t "REF (832 cm^{-1})", \
     "rad.tab" u (\$13*1.):8 w lp pt 2 t "TEST (792 cm^{-1})", \
     "rad.tab" u (\$14*1.):8 w lp pt 2 t "TEST (832 cm^{-1})"

set out "plot_rad_diff_abs.png"
set xla "Radiance difference [nW / (cm^2 sr cm^{-1})]"
set xra [*:*]
set key right top
plot "< paste rad.tab rad.org" u ((\$25-\$11)*1e5):8 w lp pt 1 t "TEST - REF (792 cm^{-1})", \
     "< paste rad.tab rad.org" u ((\$26-\$12)*1e5):8 w lp pt 1 t "TEST - REF (832 cm^{-1})"

set out "plot_rad_diff_rel.png"
set xla "Radiance difference [%]"
set xra [-5:5]
plot "< paste rad.tab rad.org" u (100.*(\$25-\$11)/\$11):8 w lp pt 1 t "TEST - REF (792 cm^{-1})", \
     "< paste rad.tab rad.org" u (100.*(\$26-\$12)/\$12):8 w lp pt 1 t "TEST - REF (832 cm^{-1})"

set out "plot_tau_diff_abs.png"
set xla "Transmittance difference"
set xra [*:*]
plot "< paste rad.tab rad.org" u (\$27-\$13):8 w lp pt 1 t "TEST - REF (792 cm^{-1})", \
     "< paste rad.tab rad.org" u (\$28-\$14):8 w lp pt 1 t "TEST - REF (832 cm^{-1})"

set out "plot_tau_diff_rel.png"
set xla "Transmittance difference [%]"
set xra [-5:5]
plot "< paste rad.tab rad.org" u (100.*(\$27-\$13)/\$13):8 w lp pt 1 t "TEST - REF (792 cm^{-1})", \
     "< paste rad.tab rad.org" u (100.*(\$28-\$14)/\$14):8 w lp pt 1 t "TEST - REF (832 cm^{-1})"
EOF

# Get differences...
echo -e "\nCheck for differences..."
diff -sq rad.tab rad.org
