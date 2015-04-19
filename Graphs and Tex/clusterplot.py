import math
import numpy as np
import matplotlib.pyplot as plt


clusters=np.loadtxt('clusters.txt')
logage=clusters[:,1]
distance=clusters[:,0]

ar=np.arange(0,20.1,0.1)
al=np.arange(0,20.1,0.1)
al[:]=math.log10(104e6)

plt.plot(ar,al,'k',linewidth=2)

plt.plot(distance, logage, 'r^')
#plt.plot(agelabel, compdatanocut, 'k--', label='No magnitude cut')
#plt.errorbar(agelabel, obsdata, yerr=error,fmt='', label='Observational data')
#legend=plt.legend(loc='upper left', shadow=True)
#plt.ylim(ymin=0.1)
plt.xlabel('distance modulus [mag]')
plt.ylabel('log($t$/yr)')
plt.show()
