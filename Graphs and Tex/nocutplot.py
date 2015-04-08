import math
import numpy as np
import matplotlib.pyplot as plt


compdata=np.loadtxt('compdata.txt')
compdata=compdata[:,1]
compdatanocut=np.loadtxt('compdatanocut.txt')
agelabel=compdatanocut[:,0]
compdatanocut=compdatanocut[:,1]
obsdata=np.loadtxt('obsdata.txt')
error=obsdata[:,3]
obsdata=obsdata[:,1]

plt.plot(agelabel, compdata, 'r', label='Magnitude cut at V=9')
plt.plot(agelabel, compdatanocut, 'k--', label='No magnitude cut')
plt.errorbar(agelabel, obsdata, yerr=error,fmt='', label='Observational data')
legend=plt.legend(loc='upper left', shadow=True)
plt.ylim(ymin=0.1)
plt.xlabel('$\mathrm{t}/\mathrm{t}_{ms}$')
plt.ylabel('probability density')
plt.show()
