import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


compdata=np.loadtxt('compdata.txt')
agelabel=compdata[:,0]
compdata=compdata[:,1]
obsdata=np.loadtxt('obsdata.txt')
error=obsdata[:,3]
obsdata=obsdata[:,1]

#plt.plot(agelabel, compdata, 'r',label='Computed data')
plt.errorbar(agelabel, obsdata, yerr=error,fmt='',label='Observational data')
#legend = plt.legend(loc='upper left', shadow=True)
plt.ylim(ymin=0.1)
plt.xlabel('$\mathrm{t}/\mathrm{t}_{ms}$')
plt.ylabel('Probability Density')
plt.show()
