# ======================================================================
# Forward model...
# ======================================================================

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
ND = 5
NU[0] = 700.0000
NU[1] = 800.0000
NU[2] = 900.0000
NU[3] = 1000.0000
NU[4] = 1100.0000

# Extinction...
NW = 3
WINDOW[0] = 0
WINDOW[1] = 0
WINDOW[2] = 1
WINDOW[3] = 1
WINDOW[4] = 2

# Cloud layer...
NCL = 2
CLNU[0] = 500.0000
CLNU[1] = 2500.0000
CLZ = 6.0
CLDZ = 1.0
CLK[0] = 1.0e-2
CLK[1] = 0.5e-2

# Surface layer...
NSF = 4
SFNU[0] = 0.0000
SFNU[1] = 1000.000
SFNU[2] = 2000.000
SFNU[3] = 3000.000
SFT = 300.0
SFEPS[0] = 1.00
SFEPS[1] = 0.95
SFEPS[2] = 0.80
SFEPS[3] = 0.50

# Raytracer...
THETA0 = -180.0
THETA1 = 180.0
DTHETA = 10.0
OBSZ = 8.0
