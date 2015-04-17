import numpy as np
import math
from scipy.constants import sigma
import Lumiradius as LR
import matplotlib.pyplot as plt

compdata=np.loadtxt('compdata.txt')
compdata=compdata[:,1]

compdata1=np.loadtxt('compdata1.txt')
agelabel=compdata1[:,0]
compdata1=compdata1[:,1]

compdata2=np.loadtxt('compdata2.txt')
compdata2=compdata2[:,1]

compdata3=np.loadtxt('compdata3.txt')
compdata3=compdata3[:,1]

compdata4=np.loadtxt('compdata4.txt')
compdata4=compdata4[:,1]

fig = plt.figure()
ax = fig.add_subplot(111)
fig.subplots_adjust(top=0.85)
ax.set_title('Graph 4')
#fig.suptitle('Graph 4', fontsize=14, fontweight='bold')
#plt.plot(agelabel, compdata,  'b',linewidth=2,label='100-1500-750')
plt.plot(agelabel, compdata1, 'b',linewidth=2,label='100-100-100')
#plt.plot(agelabel, compdata2, 'r',linewidth=2,label='10-100-100')
#plt.plot(agelabel, compdata3, 'r',linewidth=2,label='100-50-100')
plt.plot(agelabel, compdata4, 'r',linewidth=2,label='100-100-50')
plt.rc('font', size=20)
legend = plt.legend(loc='upper left', shadow=True)
plt.ylim(ymin=0.1,ymax=2.5)
plt.xlabel('$t/\mathrm{t}_{\mathrm{ms}}$')
plt.ylabel('Probability Density')
plt.show()
