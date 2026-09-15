# ======================================================================
# Forward model...
# ======================================================================

# Emissivity lookup-table directory...
TBLBASE = ../../../tests/data/airs

# Radiatively active gases...
NG = 1
EMITTER[0] = CO2

# Spectral channel center wavenumbers [cm^-1]...
ND = 3
NU[0] = 667.7820
NU[1] = 668.5410
NU[2] = 669.8110

# Kernel altitude ranges [km] for pressure (RETP), temperature (RETT),
# gases (RETQ), and extinction (RETK)...
RETP_ZMIN = -100
RETP_ZMAX = 88
RETT_ZMIN = -100
RETT_ZMAX = 88
RETQ_ZMIN[0] = -100
RETQ_ZMAX[0] = 88
RETK_ZMIN[0] = -100
RETK_ZMAX[0] = 88

# Write radiance as brightness temperature...
WRITE_BBT = 1
