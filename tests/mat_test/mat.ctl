# ======================================================================
# Forward model...
# ======================================================================

# Table directory...
TBLBASE = ../data/boxcar

# Emitters...
NG = 3
EMITTER[0] = CO2
EMITTER[1] = H2O
EMITTER[2] = O3

# Channels...
ND = 2
NU[0] = 792.0000
NU[1] = 832.0000

# Kernel...
RETT_ZMIN = 9
RETT_ZMAX = 11
RETQ_ZMIN[2] = 21
RETQ_ZMAX[2] = 25
RETK_ZMIN[0] = 3
RETK_ZMAX[0] = 6
RET_CLZ = 1
RET_CLDZ = 1
RET_CLK = 1
RET_SFT = 1
RET_SFEPS = 1

# Nadir geometry...
OBSZ = 700
LAT0 = -0.18
LAT1 = 0.18
DLAT = 0.18
