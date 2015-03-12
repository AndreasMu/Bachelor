import math
import numpy as np
import matplotlib.pyplot as plt


compdata=np.loadtxt('compdata.txt')
agelabel=compdata[:,0]
compdata=compdata[:,1]
obsdata=np.loadtxt('obsdata.txt')
error=obsdata[:,3]
obsdata=obsdata[:,1]

plt.plot(agelabel, compdata, 'r')
plt.errorbar(agelabel, obsdata, yerr=error,fmt='')
plt.xlabel('t/tms')
plt.ylabel('probability density')
plt.show()
