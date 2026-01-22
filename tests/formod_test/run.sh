#! /bin/bash

# Set environment...
export LD_LIBRARY_PATH=../../libs/build/lib:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=4
export LANG=C
export LC_ALL=C
set -euo pipefail

# Setup...
jurassic=../../src

# Create directory...
rm -rf data && mkdir -p data

# Loop over observation geometries...
for geo in limb nadir zenith ; do
    
    # Create atmospheric data file...
    $jurassic/climatology ${geo}.ctl data/atm_${geo}.tab
    
    # Create observation geomtry...
    $jurassic/${geo} ${geo}.ctl data/obs_${geo}.tab
    
    # Call forward model...
    $jurassic/formod ${geo}.ctl data/obs_${geo}.tab data/atm_${geo}.tab data/rad_${geo}.tab OBSREF data.ref/rad_${geo}.tab TASK time
    
    # Test CGA...
    $jurassic/formod ${geo}.ctl data/obs_${geo}.tab data/atm_${geo}.tab data/rad_${geo}_cga.tab OBSREF data.ref/rad_${geo}.tab FORMOD 0
    
    # Test FOV...
    if [ ${geo} = "limb" ] ; then
	$jurassic/formod ${geo}.ctl data/obs_${geo}.tab data/atm_${geo}.tab data/rad_${geo}_fov.tab OBSREF data.ref/rad_${geo}.tab FOV fov.tab
    fi
    
    # Compute kernel...
    $jurassic/kernel ${geo}.ctl data/obs_${geo}.tab data/atm_${geo}.tab data/kernel_${geo}.tab
    
done

# Compare files...
echo -e "\nCompare results..."
error=0
for f in $(ls data.ref/*.tab) ; do
    diff -q -s data/"$(basename "$f")" "$f" || error=1
done
exit $error
