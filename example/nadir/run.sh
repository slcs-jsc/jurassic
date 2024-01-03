#! /bin/bash

# Setup...
jurassic=../../src

# Create atmospheric data file...
$jurassic/climatology nadir.ctl atm.tab

# Create observation geomtry...
$jurassic/nadir nadir.ctl obs.tab

# Call forward model...
$jurassic/formod nadir.ctl obs.tab atm.tab rad.tab TASK time

# Compute kernel...
$jurassic/kernel nadir.ctl obs.tab atm.tab kernel.tab

# Plot radiance and transmittance...
gnuplot <<EOF
set term pngcairo enh truecolor font "Helvetica,28" size 1600,1200 crop lw 3
set size ratio 0.75
set title "JURASSIC | Nadir test case"

set out "plot_rad.png"
set xla "Latitude [deg]"
set yla "Brightness temperature [K]"
set mxtics
set mytics
set key spac 1.5
set key box
set key top center
set key maxrow 3
set grid
set yra [:260]
plot "rad.org" u 10:11 w lp pt 1 t "REF (667.8 cm^{-1})", \
     "rad.tab" u 10:11 w lp pt 2 t "TEST (667.8 cm^{-1})", \
     "rad.org" u 10:12 w lp pt 1 t "REF (668.5 cm^{-1})", \
     "rad.tab" u 10:12 w lp pt 2 t "TEST (668.5 cm^{-1})", \
     "rad.org" u 10:13 w lp pt 1 t "REF (669.8 cm^{-1})", \
     "rad.tab" u 10:13 w lp pt 2 t "TEST (669.8 cm^{-1})"

set out "plot_tau.png"
set yla "Transmittance"
set key bot
set yra [*:*]
plot "rad.org" u 10:14 w lp pt 1 t "REF (667.8 cm^{-1})", \
     "rad.tab" u 10:14 w lp pt 2 t "TEST (667.8 cm^{-1})", \
     "rad.org" u 10:15 w lp pt 1 t "REF (668.5 cm^{-1})", \
     "rad.tab" u 10:15 w lp pt 2 t "TEST (668.5 cm^{-1})", \
     "rad.org" u 10:16 w lp pt 1 t "REF (669.8 cm^{-1})", \
     "rad.tab" u 10:16 w lp pt 2 t "TEST (669.8 cm^{-1})"

set out "plot_rad_diff_abs.png"
set yla "Brightness temperature difference [K]"
set yra [*:*]
set key top
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

# Plot kernel functions...
for nu in 667.7820 668.5410 669.8110 ; do
    gnuplot <<EOF
set term pngcairo enh truecolor font "Helvetica,28" size 1600,1200 crop lw 3
set size ratio 0.75
set title "JURASSIC | Nadir test case"

set yla "Altitude [km]"
set cbla "View point latitude [deg]"
set cbla offset (0,1)
set cbra [0:]
set mxtics
set mytics
set grid
set pal def

set out "plot_kernel_pressure_${nu}.png"
set xla "Kernel function (pressure at $nu cm^{-1}) [K / hPa]"
plot "kernel.tab" u ((strcol(8) eq "PRESSURE" && \$2==$nu) ? \$13 : 1/0):10:6 w l lc pal z t ""

set out "plot_kernel_temperature_${nu}.png"
set xla "Kernel function (temperature at $nu cm^{-1}) [K / K]"
plot "kernel.tab" u ((strcol(8) eq "TEMPERATURE" && \$2==$nu) ? \$13 : 1/0):10:6 w l lc pal z t ""

set out "plot_kernel_CO2_${nu}.png"
set xla "Kernel function (CO_2 at $nu cm^{-1}) [K / ppmv]"
plot "kernel.tab" u ((strcol(8) eq "CO2" && \$2==$nu) ? 1e-6*\$13 : 1/0):10:6 w l lc pal z t ""

set out "plot_kernel_extinction_${nu}.png"
set xla "Kernel function (extinction at $nu cm^{-1}) [K / km^{-1}]"
plot "kernel.tab" u ((strcol(8) eq "EXTINCT_WINDOW_0" && \$2==$nu) ? \$13 : 1/0):10:6 w l lc pal z t ""
EOF
done

# Compare files...
echo -e "\nCompare results..."
error=0
diff -sq rad.tab rad.org || error=1
diff -sq kernel.tab kernel.org || error=1
exit $error
