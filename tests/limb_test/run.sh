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

# Compare files...
echo -e "\nCompare results..."
error=0
diff -sq kernel.tab kernel.org
diff -sq rad.tab rad.org || error=1
exit $error
