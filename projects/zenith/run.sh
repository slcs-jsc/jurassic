#! /bin/bash

# Set environment...
export LD_LIBRARY_PATH=../../libs/build/lib:$LD_LIBRARY_PATH

# Setup...
jurassic=../../src

# Create atmospheric data file...
$jurassic/climatology zenith.ctl atm.tab

# Create observation geomtry...
$jurassic/zenith zenith.ctl obs.tab

# Call forward model...
$jurassic/formod zenith.ctl obs.tab atm.tab rad.tab TASK time

# Compute kernel...
$jurassic/kernel zenith.ctl obs.tab atm.tab kernel.tab

# Plot radiance and transmittance...
gnuplot <<EOF
set term pngcairo enh truecolor font "Helvetica,28" size 1600,1200 crop lw 3
set size ratio 0.75
set title "JURASSIC | Zenith test case"

set out "plot_rad.png"
set xla "View point latitude [deg]"
set yla "Brightness temperature [K]"
set mxtics
set mytics
set key spac 1.5
set key box
set key top center
set key maxrow 3
set grid
set yra [:310]
plot "rad.org" u 7:11 w lp pt 1 t "REF (792 cm^{-1})", \
     "rad.tab" u 7:11 w lp pt 2 t "TEST (792 cm^{-1})", \
     "rad.org" u 7:12 w lp pt 1 t "REF (832 cm^{-1})", \
     "rad.tab" u 7:12 w lp pt 2 t "TEST (832 cm^{-1})"

set out "plot_tau.png"
set yla "Transmittance"
set yra [0:1]
plot "rad.org" u 7:13 w lp pt 1 t "REF (792 cm^{-1})", \
     "rad.tab" u 7:13 w lp pt 2 t "TEST (792 cm^{-1})", \
     "rad.org" u 7:14 w lp pt 1 t "REF (832 cm^{-1})", \
     "rad.tab" u 7:14 w lp pt 2 t "TEST (832 cm^{-1})"

set out "plot_rad_diff_abs.png"
set yla "Brightness temperature difference [K]"
set yra [*:*]
plot "< paste rad.tab rad.org" u 7:(\$25-\$11) w lp pt 1 t "TEST - REF (792 cm^{-1})", \
     "< paste rad.tab rad.org" u 7:(\$26-\$12) w lp pt 1 t "TEST - REF (832 cm^{-1})"

set out "plot_rad_diff_rel.png"
set yla "Brightness temperature difference [%]"
plot "< paste rad.tab rad.org" u 7:(100.*(\$25-\$11)/\$11) w lp pt 1 t "TEST - REF (792 cm^{-1})", \
     "< paste rad.tab rad.org" u 7:(100.*(\$26-\$12)/\$12) w lp pt 1 t "TEST - REF (832 cm^{-1})"

set out "plot_tau_diff_abs.png"
set yla "Transmittance difference"
plot "< paste rad.tab rad.org" u 7:(\$27-\$13) w lp pt 1 t "TEST - REF (792 cm^{-1})", \
     "< paste rad.tab rad.org" u 7:(\$28-\$14) w lp pt 1 t "TEST - REF (832 cm^{-1})"

set out "plot_tau_diff_rel.png"
set yla "Transmittance difference [%]"
plot "< paste rad.tab rad.org" u 7:(100.*(\$27-\$13)/\$13) w lp pt 1 t "TEST - REF (792 cm^{-1})", \
     "< paste rad.tab rad.org" u 7:(100.*(\$28-\$14)/\$14) w lp pt 1 t "TEST - REF (832 cm^{-1})"
EOF

# Plot kernel functions...
for nu in 792.0000 832.0000 ; do
    gnuplot <<EOF
set term pngcairo enh truecolor font "Helvetica,28" size 1600,1200 crop lw 3
set size ratio 0.75
set title "JURASSIC | Zenith test case"

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
plot "< awk -v nu=$nu 'NF==0 || (\$8==\"PRESSURE\" && \$2==nu)' kernel.tab" u (1.0*\$6 >= 0 ? 1.0*\$13 : 1/0):(\$10):(1.0*\$6) w l lc pal z t ""

set out "plot_kernel_temperature_${nu}.png"
set xla "Kernel function (temperature at $nu cm^{-1}) [K / K]"
plot "< awk -v nu=$nu 'NF==0 || (\$8==\"TEMPERATURE\" && \$2==nu)' kernel.tab" u (1.0*\$6 >= 0 ? 1.0*\$13 : 1/0):(\$10):(1.0*\$6) w l lc pal z t ""
EOF
done

# Compare files...
echo -e "\nCompare results..."
error=0
diff -sq rad.tab rad.org || error=1
exit $error
