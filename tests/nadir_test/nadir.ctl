# ======================================================================
# Forward model...
# ======================================================================

# Table directory...
TBLBASE = ./airs

# Emitters...
NG = 1
EMITTER[0] = CO2

# Channels...
ND = 3
NU[0] = 667.7820
NU[1] = 668.5410
NU[2] = 669.8110

# Kernel...
RETP_ZMIN = -100
RETP_ZMAX = 88
RETT_ZMIN = -100
RETT_ZMAX = 88
RETQ_ZMIN[0] = -100
RETQ_ZMAX[0] = 88
RETK_ZMIN[0] = -100
RETK_ZMAX[0] = 88

# Output...
WRITE_BBT = 1
