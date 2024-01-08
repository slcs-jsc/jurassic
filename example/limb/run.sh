#! /bin/bash

# Setup...
jurassic=../../src

# Create atmospheric data file...
$jurassic/climatology limb.ctl atm.tab

# Create observation geomtry...
$jurassic/limb limb.ctl obs.tab

# Call forward model...
$jurassic/formod limb.ctl obs.tab atm.tab rad.tab TASK time

# Compute kernel...
$jurassic/kernel limb.ctl obs.tab atm.tab kernel.tab

# Plot radiance and transmittance...
gnuplot <<EOF
set term pngcairo enh truecolor font "Helvetica,28" size 1600,1200 crop lw 3
set size ratio 0.75
set title "JURASSIC | Limb test case"

set out "plot_rad.png"
set xla "Radiance [nW / (cm^2 sr cm^{-1})]"
set yla "Tangent height [km]"
set mxtics
set mytics
set log x
set key spac 1.5
set key left bot
set key box
set grid

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

# Plot kernel functions...
for nu in 792 832 ; do
    gnuplot <<EOF
set term pngcairo enh truecolor font "Helvetica,28" size 1600,1200 crop lw 3
set size ratio 0.75
set title "JURASSIC | Limb test case"
set yla "Altitude [km]"
set cbla "View point altitude [km]"
set cbla offset (0,1)
set mxtics
set mytics
set grid
set pal def

set out "plot_kernel_pressure_${nu}.png"
set xla "Kernel function (pressure at $nu cm^{-1}) [nW / (cm^2 sr cm^{-1}) / hPa]"
plot "kernel.tab" u ((strcol(8) eq "PRESSURE" && \$2==$nu) ? 1e5*\$13 : 1/0):10:4 w lp pt 1 lc pal z t ""

set out "plot_kernel_temperature_${nu}.png"
set xla "Kernel function (temperature at $nu cm^{-1}) [nW / (cm^2 sr cm^{-1}) / K]"
plot "kernel.tab" u ((strcol(8) eq "TEMPERATURE" && \$2==$nu) ? 1e5*\$13 : 1/0):10:4 w lp pt 1 lc pal z t ""

set out "plot_kernel_CO2_${nu}.png"
set xla "Kernel function (CO_2 at $nu cm^{-1}) [nW / (cm^2 sr cm^{-1}) / ppmv]"
plot "kernel.tab" u ((strcol(8) eq "CO2" && \$2==$nu) ? 1e-6*1e5*\$13 : 1/0):10:4 w lp pt 1 lc pal z t ""

set out "plot_kernel_H2O_${nu}.png"
set xla "Kernel function (H_2O at $nu cm^{-1}) [nW / (cm^2 sr cm^{-1}) / ppmv]"
plot "kernel.tab" u ((strcol(8) eq "H2O" && \$2==$nu) ? 1e-6*1e5*\$13 : 1/0):10:4 w lp pt 1 lc pal z t ""

set out "plot_kernel_O3_${nu}.png"
set xla "Kernel function (O_3 at $nu cm^{-1}) [nW / (cm^2 sr cm^{-1}) / ppmv]"
plot "kernel.tab" u ((strcol(8) eq "O3" && \$2==$nu) ? 1e-6*1e5*\$13 : 1/0):10:4 w lp pt 1 lc pal z t ""

set out "plot_kernel_F11_${nu}.png"
set xla "Kernel function (CFC-11 at $nu cm^{-1}) [nW / (cm^2 sr cm^{-1}) / pptv]"
plot "kernel.tab" u ((strcol(8) eq "F11" && \$2==$nu) ? 1e-12*1e5*\$13 : 1/0):10:4 w lp pt 1 lc pal z t ""

set out "plot_kernel_CCl4_${nu}.png"
set xla "Kernel function (CCl_4 at $nu cm^{-1}) [nW / (cm^2 sr cm^{-1}) / pptv]"
plot "kernel.tab" u ((strcol(8) eq "CCl4" && \$2==$nu) ? 1e-12*1e5*\$13 : 1/0):10:4 w lp pt 1 lc pal z t ""

set out "plot_kernel_extinction_${nu}.png"
set xla "Kernel function (extinction at $nu cm^{-1}) [10^{6} nW / (cm^2 sr cm^{-1}) / km^{-1}]"
plot "kernel.tab" u ((strcol(8) eq "EXTINCT_WINDOW_0" && \$2==$nu) ? 1e-6*1e5*\$13 : 1/0):10:4 w lp pt 1 lc pal z t ""
EOF
done

# Compare files...
diff -sq rad.tab rad.org
