# ======================================================================
# Benchmark baseline forward model...
# ======================================================================

# Table directory...
TBLBASE = __BENCH_TBLBASE__
TBLFMT = 3

# Emitters...
NG = 7
EMITTER[0] = H2O
EMITTER[1] = CO2
EMITTER[2] = O3
EMITTER[3] = N2O
EMITTER[4] = CH4
EMITTER[5] = O2
EMITTER[6] = N2

# Channels...
ND = 32
NU[0] = 900.0000
NU[1] = 919.0000
NU[2] = 939.0000
NU[3] = 958.0000
NU[4] = 977.0000
NU[5] = 997.0000
NU[6] = 1016.0000
NU[7] = 1035.0000
NU[8] = 1055.0000
NU[9] = 1074.0000
NU[10] = 1094.0000
NU[11] = 1113.0000
NU[12] = 1132.0000
NU[13] = 1152.0000
NU[14] = 1171.0000
NU[15] = 1190.0000
NU[16] = 1210.0000
NU[17] = 1229.0000
NU[18] = 1248.0000
NU[19] = 1268.0000
NU[20] = 1287.0000
NU[21] = 1306.0000
NU[22] = 1326.0000
NU[23] = 1345.0000
NU[24] = 1365.0000
NU[25] = 1384.0000
NU[26] = 1403.0000
NU[27] = 1423.0000
NU[28] = 1442.0000
NU[29] = 1461.0000
NU[30] = 1481.0000
NU[31] = 1500.0000

# Limb geometry baseline (64 rays)...
OBSZ = 780
Z0 = 5
Z1 = 68
DZ = 1

# Kernel...
RETP_ZMIN = -100
RETP_ZMAX = 88
RETT_ZMIN = -100
RETT_ZMAX = 88
RETQ_ZMIN[0] = -100
RETQ_ZMAX[0] = 88
RETQ_ZMIN[1] = -100
RETQ_ZMAX[1] = 88
RETQ_ZMIN[2] = -100
RETQ_ZMAX[2] = 88
RETQ_ZMIN[3] = -100
RETQ_ZMAX[3] = 88
RETQ_ZMIN[4] = -100
RETQ_ZMAX[4] = 88
RETQ_ZMIN[5] = -100
RETQ_ZMAX[5] = 88
RETQ_ZMIN[6] = -100
RETQ_ZMAX[6] = 88
RETK_ZMIN[0] = -100
RETK_ZMAX[0] = 88

# Output...
WRITE_BBT = 1
