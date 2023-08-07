import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
from DDtransport import DDTransport

# ---- User input ----
# Set the mechanism file
mech = "mech.yaml"

# Set the pressure range. Please DO ensure a larger span 
# even in a constant pressure case, becase the pressure might go beyong the range 
# durng the simulation

prange = (0.7*1e5,1.5*1e5) 

# Set the temperature range. Please DO ensure a larger span
Trange = (290,3000)

# Set the target accuracy for the tabulation. The program uses a 
# first-ordered approximations to estimate the accuracy (i.e.
#  spacing between the table / the value of the table). 
# But in fact, the data distribution between two points are almost linear
# so a second-ordered accuracy can almost be guaranteed.
# Generally, when you choose 0.3, the accuracy will not be worse than 0.3, 
# but it is usually much better than 0.01.
accuracy = 0.1


# Set true if a more verbose output is expected
debug = False

# Max refinement level: If the accuracy of the table is not desired, it will refine itself.
mrf = 5

# ---- End of user input ----


# ---- Main program ----
DD = DDTransport(mech, prange, Trange, accuracy, debug)
DD.writeDiffPolyNomials()