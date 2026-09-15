# ======================================================================
# Forward model...
# ======================================================================

# Emissivity lookup-table directory...
TBLBASE = ../../../tests/data/boxcar

# Radiatively active gases...
NG = 5
EMITTER[0] = CO2
EMITTER[1] = H2O
EMITTER[2] = O3
EMITTER[3] = F11
EMITTER[4] = CCl4

# Spectral channel center wavenumbers [cm^-1]...
ND = 2
NU[0] = 792.0000
NU[1] = 832.0000

# Kernel altitude ranges [km] for pressure (RETP), temperature (RETT),
# gases (RETQ), and extinction (RETK)...
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
RETK_ZMIN[0] = -100
RETK_ZMAX[0] = 88

# Write radiance as brightness temperature...
WRITE_BBT = 1
