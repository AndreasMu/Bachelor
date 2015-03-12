import math
import numpy as np
import matplotlib.pyplot as plt


compdata=np.loadtxt('compdata.txt',skiprows=1)
compdata=compdata[:,1]
compdatanocut=np.loadtxt('compdatanocut.txt',skiprows=1)
agelabel=compdatanocut[:,0]
compdatanocut=compdatanocut[:,1]
obsdata=np.loadtxt('obsdata.txt')
error=obsdata[:,3]
obsdata=obsdata[:,1]

plt.plot(agelabel, compdata, 'r')
plt.plot(agelabel, compdatanocut, 'k--')
plt.errorbar(agelabel, obsdata, yerr=error,fmt='')
plt.xlabel('t/tms')
plt.ylabel('probability density')
plt.show()
