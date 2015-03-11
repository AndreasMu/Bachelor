import math
import numpy as np
import matplotlib.pyplot as plt


fh=open('data.txt','r')
x=np.loadtxt(fh,usecols=(0,1))

y=np.array(range(20))
y=y.astype(float)
for i in range(20):
    y[i]=x[i,1]
fh.close()

fk=open('data.txt','r')
x2=np.loadtxt(fk,usecols=(0,3))

yerr=np.array(range(20))
yerr=yerr.astype(float)
for i in range(20):
    yerr[i]=x2[i,1]
    
agelabel = np.array(range(20))
agelabel = agelabel.astype(float)
for age in range(0,20):
    agelabel[age]=age/20.

plt.figure()    
plt.errorbar(agelabel,y,yerr=y2)
