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

# Exercise isolated RFM runs with a small mock executable...
test_dir=$(pwd)
mkdir -p data/rfm_work/base/sibling
rfm_work=$test_dir/data/rfm_work
rfm_base=$rfm_work/base
mock_bin=$test_dir/mock_rfm.sh
formod_bin=$test_dir/../../src/formod
ctl_file=$test_dir/limb.ctl
obs_file=$test_dir/data/obs_limb.tab
atm_file=$test_dir/data/atm_limb.tab
tbl_base=$test_dir/../data/boxcar
export MOCK_RECORD=$rfm_work/invocations
printf 'original\n' > "$rfm_work/rfm.drv"
printf 'original\n' > "$rfm_work/rad_00000.asc"
printf 'original\n' > "$rfm_base/sibling/sentinel"

(
    cd "$rfm_work"
    "$formod_bin" "$ctl_file" "$obs_file" "$atm_file" out1.tab \
        FORMOD 2 RFMBIN "$mock_bin" RFMHIT MOCK RFMTMP "$rfm_base" \
        TBLBASE "$tbl_base" > run1.log 2>&1 &
    pid1=$!
    "$formod_bin" "$ctl_file" "$obs_file" "$atm_file" out2.tab \
        FORMOD 2 RFMBIN "$mock_bin" RFMHIT MOCK RFMTMP "$rfm_base" \
        TBLBASE "$tbl_base" > run2.log 2>&1 &
    pid2=$!
    wait "$pid1"
    wait "$pid2"
)
test "$(wc -l < "$MOCK_RECORD")" -eq 4
test "$(awk '{ print $1 }' "$MOCK_RECORD" | sort -u | wc -l)" -eq 2
test "$(find "$rfm_base" -mindepth 1 -maxdepth 1 -name 'jurassic-rfm-*' | wc -l)" -eq 0
for result in out1 out2; do
    awk '$1 !~ /^#/ && NF { if ($11 != 1 || $12 != 2) exit 1; found=1 }
         END { if (!found) exit 1 }' "$rfm_work/$result.tab"
done
test "$(cat "$rfm_work/rfm.drv")" = original
test "$(cat "$rfm_work/rad_00000.asc")" = original
test "$(cat "$rfm_base/sibling/sentinel")" = original

# The default base is the current directory; a bare executable uses PATH...
(
    cd "$rfm_work"
    PATH="$test_dir:$PATH" "$formod_bin" "$ctl_file" "$obs_file" "$atm_file" \
        default.tab FORMOD 2 RFMBIN mock_rfm.sh RFMHIT MOCK \
        TBLBASE "$tbl_base" > default.log 2>&1
)
test "$(find "$rfm_work" -mindepth 1 -maxdepth 1 -name 'jurassic-rfm-*' | wc -l)" -eq 0
test "$(cat "$rfm_work/rfm.drv")" = original

if (cd "$rfm_work" && "$formod_bin" "$ctl_file" "$obs_file" "$atm_file" \
        failed.tab FORMOD 2 RFMBIN "$mock_bin" RFMHIT FAIL RFMTMP "$rfm_base" \
        TBLBASE "$tbl_base" > failed.log 2>&1); then
    echo 'Failed mock RFM run was accepted' >&2
    exit 1
fi
failed_dir=$(find "$rfm_base" -mindepth 1 -maxdepth 1 -type d -name 'jurassic-rfm-*')
test -n "$failed_dir"
test -f "$failed_dir/rfm.drv"
test -f "$failed_dir/rfm.log"
grep -Fq "files retained in $failed_dir" "$rfm_work/failed.log"
rm -rf -- "$failed_dir"

if (cd "$rfm_work" && "$formod_bin" "$ctl_file" "$obs_file" "$atm_file" \
        stale.tab FORMOD 2 RFMBIN "$mock_bin" RFMHIT STALE RFMTMP "$rfm_base" \
        TBLBASE "$tbl_base" > stale.log 2>&1); then
    echo 'Stale RFM spectra were accepted' >&2
    exit 1
fi
grep -q 'Cannot find RFM data file!' "$rfm_work/stale.log"
stale_dir=$(find "$rfm_base" -mindepth 1 -maxdepth 1 -type d -name 'jurassic-rfm-*')
test -n "$stale_dir"
rm -rf -- "$stale_dir"
test -d "$rfm_base/sibling"

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
