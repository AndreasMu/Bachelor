import numpy as np
import math
from scipy.constants import sigma
#import matplotlib.pyplot as plt

rangeage=100

rangemass=100
minmass=5.
maxmass=50.
dlogM = (math.log10(maxmass)-math.log10(minmass))/rangemass

norm=np.zeros(rangemass)
norm=norm.astype(float)

def Mainsequenceage(mass):
    a1 = 1593.89
    a2 = 2706.708
    a3 = 146.6143
    a4 = 0.0414196
    a5 = 0.3426349
    a6 = 19.49814
    a7 = 4.90383
    a8 = 0.05212154
    a9 = 1.312179
    a10 = 0.8073972
    k = max(0.95,min(0.95-0.03*0.30103,0.99))
    mu = max(0.5,1.0-0.01*max(a6/mass**a7,a8+a9/mass**a10))
    tbgb = (a1+a2*mass**4+a3*mass**5.5+mass**7)/(a4*mass**2+a5*mass**7)
    thook = mu*tbgb
    tms = max(thook,k*tbgb)
    return tms

maxage = Mainsequenceage(minmass)
dlogt=(math.log10(maxage)+1)/rangeage
print("dlogt="+str(dlogt)+"\n")

def Agerelation(mass,realage):
    #This functions purpose is to make sure I don't include stars, that are already dead.
    #If the age of the star is greater, than its mainsequence age, the function will return 0.
    tms = Mainsequenceage(mass)
    if(realage>tms):
        agefactor = 0
    else:
        agefactor = math.log(10)*realage/tms
    return agefactor

maxdiff=0.
agestuff=np.array(range(rangeage))
agestuff=agestuff.astype(float)
for mass in range(0, rangemass):
    masss = 10**(math.log10(minmass) + mass*dlogM)
    #dlogt=(math.log10(Mainsequenceage(masss))/rangeage)
    for age in range(0,rangeage):
        realage = 10**(-1+age*dlogt)
        agestuff[age]=realage
        dp = Agerelation(masss,realage)*dlogt
        norm[mass] += dp
    if abs(1-norm[mass])>maxdiff:
        maxdiff=1-norm[mass]
