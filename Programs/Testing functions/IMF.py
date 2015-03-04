import numpy as np
import math
from scipy.constants import sigma
import matplotlib.pyplot as plt

rangemass = 100 #Amount of intervals I have for mass
maxmass = 50.
minmass = 5.
dlogM = (math.log10(maxmass)-math.log10(minmass))/rangemass
print("dlogM="+str(dlogM)+"\n")

bins = np.array(range(rangemass))
masslabel = bins.astype(float)
massfunction = np.array(range(rangemass))
massfunction = massfunction.astype(float)

def IMF(mass):
    e= 1.35 * (1/(minmass**-1.35-maxmass**-1.35))
    massfactor =  math.log(10)*e*mass**(-1.35)
    return (massfactor)

norm=0.
for mass in range(0,rangemass):
    masss = 10**(math.log10(minmass) + mass*dlogM)
    #masss = minmass*10**((float(mass))/rangemass)
    #print("M="+str(math.log10(masss))+"\n")
    masslabel[mass]=masss
    dp = IMF(masss)*dlogM
    norm += dp
    massfunction[mass]=dp

print("Norm="+str(norm)+"\n")
plt.plot(masslabel, massfunction)
plt.xlabel('mass')
plt.ylabel('dp')
plt.show()
