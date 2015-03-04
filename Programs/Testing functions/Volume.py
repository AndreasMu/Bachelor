import numpy as np
import math
from scipy.constants import sigma
#import matplotlib.pyplot as plt

rangedistance = 1000
maxdistance = 5.
dr = float(maxdistance)/rangedistance

bins = np.array(range(rangedistance))
volumelabel = bins.astype(float)
volumeintegral = np.array(range(rangedistance))
volumeintegral = volumeintegral.astype(float)

def DpdV(distance):
    dpdv=1/((4./3)*math.pi*maxdistance**3)
    return dpdv


norm=0.
for distance in range(0,rangedistance):
    volumelabel[distance]=Volume(distance)
    dV=4*math.pi*dr**3*distance**2
    #dV = ((distance+1)**3-(distance)**3)*(4./3)*math.pi*dr**3
    dp=DpdV(distance)*dV
    norm+=dp


print("Norm="+str(norm)+"\n")
"""plt.plot(bins, volumeintegral)
plt.xlabel('bin')
plt.ylabel('integral of volume')
plt.show()"""
