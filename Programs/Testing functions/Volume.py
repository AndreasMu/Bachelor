import numpy as np
import math
from scipy.constants import sigma
#import matplotlib.pyplot as plt

rangedistance = 1000
maxdistance = 5.
dr = float(maxdistance)/rangedistance
dV = (4./3)*math.pi*dr**3

bins = np.array(range(rangedistance))
volumelabel = bins.astype(float)
volumeintegral = np.array(range(rangedistance))
volumeintegral = volumeintegral.astype(float)

def Volume(distance):
    distance+=1.
    volumenorm=1/((4./3)*math.pi*maxdistance**3)
    volume = volumenorm*(distance**3-(distance-1)**3)
    #volume = 4*math.pi*dr**3*distance**3
    #The Volume of a spherical shell is 4/3*pi*(a^3-b^3) a>b. So this is the
    #volume of the nth shell.
    return volume

norm=0.
for distance in range(0,rangedistance):
    volumelabel[distance]=Volume(distance)
    dp=Volume(distance)*dV
    norm+=dp


print("Norm="+str(norm)+"\n")
"""plt.plot(bins, volumeintegral)
plt.xlabel('bin')
plt.ylabel('integral of volume')
plt.show()"""
