#! /bin/bash

# Create observation geomtry...
../src/nadir nadir.ctl obs.tab

# Create atmospheric data file...
../src/climatology nadir.ctl atm.tab

# Call forward model...
../src/formod nadir.ctl obs.tab atm.tab rad.tab TASK time

# Get differences...
echo -e "\nCheck for differences..."
diff -sq rad.tab rad.org

