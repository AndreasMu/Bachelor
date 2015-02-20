import numpy as np
import math
from scipy.constants import sigma
import Lumiradius as LR
import matplotlib.pyplot as plt

rangemass = 750 #Amount of intervals I have for mass
maxmass = 50.
minmass = 5.

def IMF(mass):
    #The IMF is: IMF=e*mass**{-2.3} To get e i calculate the integral IMF from minmass to 100.
    #1==int_{minmass}^{100}(e*m^{-2.3}}=\frac{-1}{1.3}*e*100^{-1.3} + \frac{1}{1.3}*e*minmass^{-1.3}
    #because minmass^{-1.3}>>100^{-1.3} I neglect the first term. Solving for e gives:
    e= 1.35*minmass**1.35
    massfactor = e*mass**(-2.35)
    return massfactor
