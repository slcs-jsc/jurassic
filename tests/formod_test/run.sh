#! /bin/bash

# Set environment...
export LD_LIBRARY_PATH=../../libs/build/lib:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=4
export LANG=C
export LC_ALL=C
set -euo pipefail

# Setup...
jurassic=../../src
refdir=${JURASSIC_REF_DIR:-data.ref}

# Create directory...
rm -rf data && mkdir -p data

# Loop over observation geometries...
for geo in limb nadir zenith ; do
    
    # Create atmospheric data file...
    $jurassic/climatology ${geo}.ctl data/atm_${geo}.tab
    
    # Create observation geometry...
    $jurassic/${geo} ${geo}.ctl data/obs_${geo}.tab
    
    # Call forward model...
    $jurassic/formod ${geo}.ctl data/obs_${geo}.tab data/atm_${geo}.tab data/rad_${geo}.tab OBSREF "$refdir/rad_${geo}.tab" TASK time
    
    # Test CGA...
    $jurassic/formod ${geo}.ctl data/obs_${geo}.tab data/atm_${geo}.tab data/rad_${geo}_cga.tab OBSREF "$refdir/rad_${geo}.tab" FORMOD 0
    
    # Test FOV...
    if [ ${geo} = "limb" ] ; then
	$jurassic/formod ${geo}.ctl data/obs_${geo}.tab data/atm_${geo}.tab data/rad_${geo}_fov.tab OBSREF "$refdir/rad_${geo}.tab" FOV fov.tab
    fi
    
    # Compute kernel...
    $jurassic/kernel ${geo}.ctl data/obs_${geo}.tab data/atm_${geo}.tab data/kernel_${geo}.tab
    
done

# Check the field-of-view shape reader at and beyond its capacity...
nshape_limit=$(awk '$1 == "#define" && $2 == "NSHAPE" { print $3; exit }' "$jurassic/jurassic.h")
awk -v n="$nshape_limit" 'BEGIN { for (i = 0; i < n; i++) print i, 1 }' \
    > data/shape_boundary_max.tab
echo '# trailing comment' >> data/shape_boundary_max.tab
$jurassic/obsfmt limb.ctl data/obs_limb.tab 1 data/obs_shape_boundary.tab 1 \
    FOV data/shape_boundary_max.tab > data/shape_boundary_max.log 2>&1
grep -q "Number of data points: $nshape_limit" data/shape_boundary_max.log
cp data/shape_boundary_max.tab data/shape_boundary_over.tab
echo "$nshape_limit 1" >> data/shape_boundary_over.tab
if $jurassic/obsfmt limb.ctl data/obs_limb.tab 1 data/obs_shape_rejected.tab 1 \
    FOV data/shape_boundary_over.tab > data/shape_boundary_over.log 2>&1; then
    echo 'Oversized shape function was accepted' >&2
    exit 1
fi
grep -q 'Too many data points!' data/shape_boundary_over.log

# Compare files...
echo -e "\nCompare results..."
error=0
for f in "$refdir"/*.tab ; do
    diff -q -s data/"$(basename "$f")" "$f" || error=1
done
exit $error
