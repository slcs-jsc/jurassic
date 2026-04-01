#! /bin/bash

# Set environment...
export LD_LIBRARY_PATH=../../libs/build/lib:${LD_LIBRARY_PATH:-}
export LANG=C
export LC_ALL=C
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-2}
set -euo pipefail

# Setup...
jurassic=../../src
ctl=../ret_test/ret.ctl
NPROF=${NPROF:-4}
MPI_NTASKS=${MPI_NTASKS:-2}
MPI_CMD=${MPI_CMD:-mpirun -np ${MPI_NTASKS}}
read -r -a mpi_cmd <<< "$MPI_CMD"

# Create directory...
rm -rf data && mkdir -p data/shared
: > data/dirlist.txt
: > data/proflist.txt

# Create atmospheric and observation data files...
for ((ip = 0; ip < NPROF; ip++)); do
    printf -v tag "%03d" "$ip"
    casedir="data/case${tag}"
    mkdir -p "$casedir"
    printf "%s\n" "$casedir" >> data/dirlist.txt
    printf "%d\n" "$ip" >> data/proflist.txt
    obsz=$((780 - ip))

    $jurassic/climatology "$ctl" "$casedir/atm_apr.tab"
    $jurassic/limb "$ctl" "$casedir/obs.tab" OBSZ "$obsz"
    $jurassic/formod "$ctl" "$casedir/obs.tab" "$casedir/atm_apr.tab" \
        "$casedir/obs_meas.tab"
done

# Retrieval in legacy directory mode...
$jurassic/retrieval "$ctl" data/dirlist.txt

# Create shared netCDF retrieval inputs...
for ((ip = 0; ip < NPROF; ip++)); do
    printf -v tag "%03d" "$ip"
    casedir="data/case${tag}"
    $jurassic/atmfmt "$ctl" "$casedir/atm_apr.tab" 1 data/shared/atm_apr.nc 3 \
        PROF_OUT "$ip"
    $jurassic/obsfmt "$ctl" "$casedir/obs_meas.tab" 1 data/shared/obs_meas.nc 3 \
        PROF_OUT "$ip"
done

# Retrieval in shared netCDF mode...
"${mpi_cmd[@]}" "$jurassic/retrieval" "$ctl" - \
    SHARED_IO_PROFLIST data/proflist.txt \
    ATMFMT 3 OBSFMT 3 MATRIXFMT 3 \
    SHARED_IO_ATM_APR_FILE data/shared/atm_apr.nc \
    SHARED_IO_OBS_MEAS_FILE data/shared/obs_meas.nc \
    SHARED_IO_ATM_FINAL_FILE data/shared/atm_final.nc \
    SHARED_IO_OBS_FINAL_FILE data/shared/obs_final.nc \
    SHARED_IO_MATRIX_COV_APR_FILE data/shared/matrix_cov_apr.nc \
    SHARED_IO_MATRIX_KERNEL_FILE data/shared/matrix_kernel.nc \
    SHARED_IO_MATRIX_COV_RET_FILE data/shared/matrix_cov_ret.nc \
    SHARED_IO_MATRIX_CORR_FILE data/shared/matrix_corr.nc \
    SHARED_IO_MATRIX_GAIN_FILE data/shared/matrix_gain.nc \
    SHARED_IO_MATRIX_AVK_FILE data/shared/matrix_avk.nc \
    SHARED_IO_ATM_ERR_TOTAL_FILE data/shared/atm_err_total.nc \
    SHARED_IO_ATM_ERR_NOISE_FILE data/shared/atm_err_noise.nc \
    SHARED_IO_ATM_ERR_FORMOD_FILE data/shared/atm_err_formod.nc \
    SHARED_IO_ATM_CONT_FILE data/shared/atm_cont.nc \
    SHARED_IO_ATM_RES_FILE data/shared/atm_res.nc

# Extract shared outputs and compare...
for ((ip = 0; ip < NPROF; ip++)); do
    printf -v tag "%03d" "$ip"
    casedir="data/case${tag}"
    atm_out="data/shared/atm_final_p${tag}.tab"
    obs_out="data/shared/obs_final_p${tag}.tab"
    mat_out="data/shared/matrix_kernel_p${tag}.tab"

    $jurassic/atmfmt "$ctl" data/shared/atm_final.nc 3 "$atm_out" 1 \
        PROF_IN "$ip"
    $jurassic/obsfmt "$ctl" data/shared/obs_final.nc 3 "$obs_out" 1 \
        PROF_IN "$ip"
    $jurassic/matfmt "$ctl" "$atm_out" "$obs_out" y x r \
        data/shared/matrix_kernel.nc 3 "$mat_out" 1 PROF_IN "$ip"

    diff -q -s "$atm_out" "$casedir/atm_final.tab"
    diff -q -s "$obs_out" "$casedir/obs_final.tab"
    diff -q -s "$mat_out" "$casedir/matrix_kernel.tab"
done
