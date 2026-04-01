#! /bin/bash

# Set environment...
export LD_LIBRARY_PATH=../../libs/build/lib:${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=4
export LANG=C
export LC_ALL=C
set -euo pipefail

# Setup...
jurassic=../../src

# Create directory...
rm -rf data && mkdir -p data/case0 data/case1 data/shared
printf "data/case0\ndata/case1\n" > data/dirlist.txt
printf "0\n1\n" > data/proflist.txt

# Create atmospheric data files...
$jurassic/climatology ret.ctl data/case0/atm_apr.tab
$jurassic/climatology ret.ctl data/case1/atm_apr.tab

# Create observation geomtry...
$jurassic/limb ret.ctl data/case0/obs.tab
$jurassic/limb ret.ctl data/case1/obs.tab OBSZ 760

# Call forward model...
$jurassic/formod ret.ctl data/case0/obs.tab data/case0/atm_apr.tab \
    data/case0/obs_meas.tab
$jurassic/formod ret.ctl data/case1/obs.tab data/case1/atm_apr.tab \
    data/case1/obs_meas.tab

# Retrieval in legacy directory mode...
$jurassic/retrieval ret.ctl data/dirlist.txt

# Create shared netCDF retrieval inputs...
$jurassic/atmfmt ret.ctl data/case0/atm_apr.tab 1 data/shared/atm_apr.nc 3 \
    PROF_OUT 0
$jurassic/atmfmt ret.ctl data/case1/atm_apr.tab 1 data/shared/atm_apr.nc 3 \
    PROF_OUT 1
$jurassic/obsfmt ret.ctl data/case0/obs_meas.tab 1 data/shared/obs_meas.nc 3 \
    PROF_OUT 0
$jurassic/obsfmt ret.ctl data/case1/obs_meas.tab 1 data/shared/obs_meas.nc 3 \
    PROF_OUT 1

# Retrieval in shared netCDF mode...
$jurassic/retrieval ret.ctl - PROFLIST data/proflist.txt \
    ATMFMT 3 OBSFMT 3 MATRIXFMT 3 \
    ATM_APR_FILE data/shared/atm_apr.nc \
    OBS_MEAS_FILE data/shared/obs_meas.nc \
    ATM_FINAL_FILE data/shared/atm_final.nc \
    OBS_FINAL_FILE data/shared/obs_final.nc \
    MATRIX_COV_APR_FILE data/shared/matrix_cov_apr.nc \
    MATRIX_KERNEL_FILE data/shared/matrix_kernel.nc \
    MATRIX_COV_RET_FILE data/shared/matrix_cov_ret.nc \
    MATRIX_CORR_FILE data/shared/matrix_corr.nc \
    MATRIX_GAIN_FILE data/shared/matrix_gain.nc \
    MATRIX_AVK_FILE data/shared/matrix_avk.nc \
    ATM_ERR_TOTAL_FILE data/shared/atm_err_total.nc \
    ATM_ERR_NOISE_FILE data/shared/atm_err_noise.nc \
    ATM_ERR_FORMOD_FILE data/shared/atm_err_formod.nc \
    ATM_CONT_FILE data/shared/atm_cont.nc \
    ATM_RES_FILE data/shared/atm_res.nc

# Extract shared outputs for comparison...
$jurassic/atmfmt ret.ctl data/shared/atm_final.nc 3 data/shared/atm_final_p0.tab 1 \
    PROF_IN 0
$jurassic/atmfmt ret.ctl data/shared/atm_final.nc 3 data/shared/atm_final_p1.tab 1 \
    PROF_IN 1
$jurassic/obsfmt ret.ctl data/shared/obs_final.nc 3 data/shared/obs_final_p0.tab 1 \
    PROF_IN 0
$jurassic/obsfmt ret.ctl data/shared/obs_final.nc 3 data/shared/obs_final_p1.tab 1 \
    PROF_IN 1
$jurassic/matfmt ret.ctl data/shared/atm_final_p0.tab data/shared/obs_final_p0.tab \
    y x r data/shared/matrix_kernel.nc 3 data/shared/matrix_kernel_p0.tab 1 \
    PROF_IN 0
$jurassic/matfmt ret.ctl data/shared/atm_final_p1.tab data/shared/obs_final_p1.tab \
    y x r data/shared/matrix_kernel.nc 3 data/shared/matrix_kernel_p1.tab 1 \
    PROF_IN 1

# Compare files...
echo -e "\nCompare results..."
error=0
for f in data.ref/*.tab ; do
    diff -q -s data/case0/"$(basename "$f")" "$f" || error=1
done
diff -q -s data/shared/atm_final_p0.tab data/case0/atm_final.tab || error=1
diff -q -s data/shared/atm_final_p1.tab data/case1/atm_final.tab || error=1
diff -q -s data/shared/obs_final_p0.tab data/case0/obs_final.tab || error=1
diff -q -s data/shared/obs_final_p1.tab data/case1/obs_final.tab || error=1
diff -q -s data/shared/matrix_kernel_p0.tab data/case0/matrix_kernel.tab \
    || error=1
diff -q -s data/shared/matrix_kernel_p1.tab data/case1/matrix_kernel.tab \
    || error=1
exit $error
