#! /bin/bash

# Create atmospheric data file...
../src/climatology nadir.ctl atm.tab

# Create observation geomtry...
../src/limb nadir.ctl obs_limb.tab
../src/nadir nadir.ctl obs_nadir.tab

# Call forward model...
../src/formod nadir.ctl obs_limb.tab atm.tab rad_limb.tab TASK time
../src/formod nadir.ctl obs_nadir.tab atm.tab rad_nadir.tab TASK time

# Get differences...
echo -e "\nCheck for differences..."
diff -sq rad_limb.tab rad_limb.org
diff -sq rad_nadir.tab rad_nadir.org
