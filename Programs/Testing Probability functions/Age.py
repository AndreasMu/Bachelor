import numpy as np
import math
from scipy.constants import sigma
import Lumiradius as LR
import matplotlib.pyplot as plt

rangeage=1000
rangemass=1000
minmass=5.
maxmass=50.

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

def Agerelation(mass,realage):
    #This functions purpose is to make sure I don't include stars, that are already dead.
    #If the age of the star is greater, than its mainsequence age, the function will return 0.
    tms = Mainsequenceage(mass)
    if(realage>tms):
        agefactor = 0
    else:
        agefactor = 1./tms
    return agefactor

