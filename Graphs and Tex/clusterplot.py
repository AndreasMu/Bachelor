import math
import numpy as np
import matplotlib.pyplot as plt


clusters=np.loadtxt('clusters.txt')
logage=clusters[:,1]
distance=clusters[:,0]

plt.plot(distance, logage, 'r^')
#plt.plot(agelabel, compdatanocut, 'k--', label='No magnitude cut')
#plt.errorbar(agelabel, obsdata, yerr=error,fmt='', label='Observational data')
#legend=plt.legend(loc='upper left', shadow=True)
#plt.ylim(ymin=0.1)
plt.xlabel('distance modulus [mag]')
plt.ylabel('log(Age/yr)')
plt.show()
