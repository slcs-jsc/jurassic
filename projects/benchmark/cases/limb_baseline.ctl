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
NU[1] = 919.3548
NU[2] = 938.7097
NU[3] = 958.0645
NU[4] = 977.4194
NU[5] = 996.7742
NU[6] = 1016.1290
NU[7] = 1035.4839
NU[8] = 1054.8387
NU[9] = 1074.1935
NU[10] = 1093.5484
NU[11] = 1112.9032
NU[12] = 1132.2581
NU[13] = 1151.6129
NU[14] = 1170.9677
NU[15] = 1190.3226
NU[16] = 1209.6774
NU[17] = 1229.0323
NU[18] = 1248.3871
NU[19] = 1267.7419
NU[20] = 1287.0968
NU[21] = 1306.4516
NU[22] = 1325.8065
NU[23] = 1345.1613
NU[24] = 1364.5161
NU[25] = 1383.8710
NU[26] = 1403.2258
NU[27] = 1422.5806
NU[28] = 1441.9355
NU[29] = 1461.2903
NU[30] = 1480.6452
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
