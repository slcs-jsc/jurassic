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

# Compare files...
echo -e "\nCompare results..."
error=0
diff -sq kernel.tab kernel.org
diff -sq rad.tab rad.org || error=1
exit $error
