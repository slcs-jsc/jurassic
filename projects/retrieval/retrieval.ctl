# Synthetic EGA limb retrieval; paths are relative to projects/retrieval.

# Lookup tables and forward model: packed 1 cm^-1 TRIA emissivity tables.
TBLBASE = ../../tab/tria_1cm/nc_1e-6/tria
TBLFMT = 3
FORMOD = 1

# Radiatively active gases; O3 and HNO3 are state-vector emitters 2 and 3.
NG = 4
EMITTER[0] = CO2
EMITTER[1] = H2O
EMITTER[2] = O3
EMITTER[3] = HNO3

# Spectral channels [cm^-1]: CO2 region, O3 region, HNO3 region.
ND = 15
NU[0] = 791
NU[1] = 792
NU[2] = 793
NU[3] = 810
NU[4] = 811
NU[5] = 812
NU[6] = 778
NU[7] = 779
NU[8] = 780
NU[9] = 787
NU[10] = 788
NU[11] = 863
NU[12] = 864
NU[13] = 867
NU[14] = 869

# Retrieved state vector: temperature, O3, and HNO3 at 15--45 km;
# pressure and background CO2/H2O stay fixed.
RETT_ZMIN = 15
RETT_ZMAX = 45
RETQ_ZMIN[2] = 15
RETQ_ZMAX[2] = 45
RETQ_ZMIN[3] = 15
RETQ_ZMAX[3] = 45

# A priori 1-sigma errors: K for temperature, percent of VMR for gases.
# Vertical covariance decays exponentially over 5 km. The horizontal
# length is effectively irrelevant for this single-profile experiment.
ERR_TEMP = 5
ERR_TEMP_CZ = 5
ERR_TEMP_CH = 100000
ERR_Q[2] = 30
ERR_Q_CZ[2] = 5
ERR_Q_CH[2] = 100000
ERR_Q[3] = 30
ERR_Q_CZ[3] = 5
ERR_Q_CH[3] = 100000

# Assumed fit errors: absolute radiance noise [W/(m^2 sr cm^-1)] and
# forward-model error [% of observed radiance]. No noise is added to truth.
ERR_NOISE[*] = 0.0001
ERR_FORMOD[*] = 0.1

# Standard convergence controls (made explicit for reproducibility).
CONV_ITMAX = 30
CONV_DMIN = 0.1
KERNEL_RECOMP = 3

# Write finite-difference kernels, averaging kernels, and covariances.
ERR_ANA = 1
WRITE_MATRIX = 1
