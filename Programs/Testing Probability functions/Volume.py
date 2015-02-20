import numpy as np
import math
from scipy.constants import sigma
import Lumiradius as LR
import matplotlib.pyplot as plt

rangedistance = 1000
maxdistance = 5.
stepdistance = float(maxdistance)/rangedistance

bins = np.array(range(rangedistance))
volumelabel = bins.astype(float)
volumeintegral = np.array(range(rangedistance))
volumeintegral = volumeintegral.astype(float)

def Volume(distance):
    distance+=1.
    volumenorm=1/((4./3)*math.pi*maxdistance**3)
    volume = volumenorm*(distance**3-(distance-1)**3)*(4./3)*math.pi*stepdistance**3
    #volume = 4*math.pi*stepdistance**3*distance**3
    #The Volume of a spherical shell is 4/3*pi*(a^3-b^3) a>b. So this is the
    #volume of the nth shell.
    return volume

for distance in range(0,rangedistance):
    volumelabel[distance]=Volume(distance)
    if distance==0:
        volumeintegral[distance]=Volume(distance)
    else:
        volumeintegral[distance]=volumeintegral[distance-1]+Volume(distance)

plt.plot(bins, volumeintegral)
plt.xlabel('bin')
plt.ylabel('integral of volume')
plt.show()
