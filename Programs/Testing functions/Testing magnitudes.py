import numpy as np
import math
from scipy.constants import sigma
import Lumiradius as LR
import matplotlib.pyplot as plt

Z = 0.02
#The metalicity should be 0.02

fh=open('table.txt','r')
Data=np.loadtxt(fh,skiprows=1,usecols=(2,14,5,6,9,12,13,15,20))
Massdata=Data[:,2]
magdata=Data[:,0]
Teffdata=Data[:,6]
logLdata=Data[:,1]
Agedata=Data[:,3]
Fracmsdata=Data[:,4]
Rdata=Data[:,5]
Distdata=Data[:,7]
BCdata=Data[:,8]
minmass=10
for i in range(len(Massdata)):
    if minmass > Massdata[i]:
        minmass=Massdata[i]

logRdata=np.array(range(len(Rdata)))
logRdata=logRdata.astype(float)
for i in range(len(Rdata)):
    logRdata[i]=math.log10(Rdata[i])
for i in range(len(Rdata)):
    Distdata[i]=Distdata[i]/1000

#This function is an excerpt from Lumiradius.py. Since Z=0.02 in our case the formulas
#get rather simple in comparison.
#tms is in units of Myr
def Mainsequenceage(realmass):
    a1 = 1593.89
    a2 = 2706.708
    a3 = 146.6143
    a4 = 0.0414196
    a5 = 0.3426349
    a6 = 19.49814
    a7 = 4.90383
    a8 = 0.05212154
    a9 = 1.312179
    a10 = 0.8073972
    k = max(0.95,min(0.95-0.03*0.30103,0.99))
    mu = max(0.5,1.0-0.01*max(a6/realmass**a7,a8+a9/realmass**a10))
    tbgb = (a1+a2*realmass**4+a3*realmass**5.5+realmass**7)/(a4*realmass**2+a5*realmass**7)
    thook = mu*tbgb
    tms = max(thook,k*tbgb)
    return tms

maxage = Mainsequenceage(minmass)


#Here I make two matrices, which will contain Luminosity and Radius of stars of
#a specific age and mass. In the later parts of this program I can then simply
#use these matrices instead of always having to call the function from Lumiradius.py
logLcomp = np.arange(len(Rdata),dtype=np.float)
logRcomp = np.arange(len(Rdata),dtype=np.float)
for i in range(len(Rdata)):
    logLcomp[i], logRcomp[i]= LR.Lumiradius(Massdata[i], Z, Fracmsdata[i])


#This is the start of calculating the visual magnitude of the star.
def Temperature(logL,logR):
    #I know: L=sigma*A*T**4 and A=4*pi*R**2
    #Therefore: T=sqrt(sqrt(L/(sigma*4*pi*R**2)))
    temperature=math.sqrt(math.sqrt((10**logL*3.846e26)/(sigma*4*math.pi*(10**logR*696342000)**2)))
    return temperature
    
def BC(T):
    if math.log10(T)>4:
        bc=4.1940953-0.00070441042*T+3.4516521e-8*T**2-9.5565244e-13*T**3+1.2790825e-17*T**4-6.4741275e-23*T**5
    elif math.log10(T)>3.7:
        bc=-29.325541+0.018052720*T-4.4823439e-6*T**2+5.5894085e-10*T**3-3.4753865e-14*T**4+8.5372998e-19*T**5
    else:
        bc=-210.13793+0.19596489*T-7.4465325e-5*T**2+1.4337726e-8*T**3-1.3955426e-12*T**4+5.4925758e-17*T**5
    return bc

def Reddening(distance):
    #Reddening is obtained from fig.9 from Amores_Lepine_2005 interpolating between known values
    d1=0.9
    d2=2.25
    d5=5.3
    if distance<1:
        red=0.9*distance
    elif distance<2:
        red=0.9+1.35*(distance-1)
    else:
        red=2.25+1.023*(distance-2)
    return red

def Magnitude(logL,logR,distance):
    red=Reddening(distance)
    T=Temperature(logL,logR)
    bc=BC(T)
    V=5*math.log10(1000*distance)-5+red+4.72-logL/0.4-bc
    return V

magdiff1=np.arange(len(Rdata),dtype=np.float)
magcomp1=np.arange(len(Rdata),dtype=np.float)
magdiff2=np.arange(len(Rdata),dtype=np.float)
magcomp2=np.arange(len(Rdata),dtype=np.float)
maxdiff1=0
maxdiff2=0
averagediff1=0.
averagediff2=0.
for i in range(len(logRdata)):
    magcomp2[i]=Magnitude(logLcomp[i],logRcomp[i],Distdata[i])
    magcomp1[i]=Magnitude(logLdata[i],logRdata[i],Distdata[i])
    magdiff2[i]=magcomp2[i]-magdata[i]
    magdiff1[i]=magcomp1[i]-magdata[i]
    if maxdiff2<magdiff2[i]:
        maxdiff2=magdiff2[i]
    if maxdiff1<magdiff1[i]:
        maxdiff1=magdiff1[i]
    averagediff2+=abs(magdiff2[i])
    averagediff1+=abs(magdiff1[i])
averagediff1=averagediff1/len(Rdata)
averagediff2=averagediff2/len(Rdata)

np.savetxt('magdiff2.txt',magdiff2,fmt='%5.3f')
np.savetxt('magdiff1.txt',magdiff1,fmt='%5.3f')

for i in range(len(Rdata)):
    Distdata[i]=math.log10(Distdata[i])
    
plt.plot(logLdata, magdiff1, 'r^')
#plt.plot(Distdata, magdiff2, 'r^')
#plt.xlim(xmin=-1.4999)
plt.xlabel('log(L/L$_\odot$)')
plt.ylabel('$V_{comp}-V_{obs}$')
plt.show()
