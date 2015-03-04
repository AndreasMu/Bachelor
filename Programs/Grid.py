import numpy as np
import math
from scipy.constants import sigma
import Lumiradius as LR
import matplotlib.pyplot as plt

Z = 0.02
#The metalicity should be 0.02
numberfactor = 100000.
#The numberfactor is what I use to regulate the total number of stars

rangedistance = 100 #The amount of intervals I have for distance
maxdistance = 5.#The maximum distance in kpc
dr = float(maxdistance)/rangedistance
rho=1

rangemass = 100 #Amount of intervals I have for mass
maxmass = 50.
minmass = 5.
dlogM = (math.log10(maxmass)-math.log10(minmass))/rangemass

rangeage = 100 #Amount of intervals I have for age


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
dlogt=(math.log10(maxage)+1)/rangeage

def Volume(distance):
    distance+=1.
    volumenorm=1/((4./3)*math.pi*maxdistance**3)
    volume = volumenorm*(distance**3-(distance-1)**3)
    #volume = 4*math.pi*dr**3*distance**3
    #The Volume of a spherical shell is 4/3*pi*(a^3-b^3) a>b. So this is the
    #volume of the nth shell.
    return volume

"""volumelabel = np.array(range(rangedistance))
volumelabel = volumelabel.astype(float)
for distance in range(0,rangedistance):
    volumelabel[distance]=Volume(distance)
    #The volumelabel will save the fractional volume of each shell."""

def IMF(realmass):
    e= 1.35 * (1/(minmass**-1.35-maxmass**-1.35))
    massfactor =  math.log(10)*e*realmass**(-1.35)
    return (massfactor)

def Agerelation(realmass,realage):
    tms = Mainsequenceage(realmass)
    if(realage>tms):
        agefactor = 0
    else:
        agefactor = math.log(10)*realage/tms
    return agefactor

#Here an agerelation matrix is filled so I can easily check, whether a star of
#mass x age y actually exists. 
agerelation = np.zeros(rangemass*rangeage)
agerelation = agerelation.astype(int)
agerelation = agerelation.reshape((rangemass,rangeage))
for mass in range(0,rangemass):
    realmass = 10**(math.log10(minmass) + mass*dlogM)
    for age in range(0,rangeage):
        realage = 10**(-1.+age*dlogt)
        agerelation[mass,age]=Agerelation(realmass,realage)

#Here I make two matrices, which will contain Luminosity and Radius of stars of
#a specific age and mass. In the later parts of this program I can then simply
#use these matrices instead of always having to call the function from Lumiradius.py
logK = np.array(range(rangemass*rangeage))
logL = logK.reshape((rangemass,rangeage))
logR = logK.reshape((rangemass,rangeage))
logL = logL.astype(float)
logR = logR.astype(float)
for mass in range(0,rangemass):
    realmass = 10**(math.log10(minmass) + mass*dlogM)
    for age in range(0,rangeage):
        realage = 10**(-1.+age*dlogt)
        fracage = realage/Mainsequenceage(realmass)
        if (agerelation[mass,age]!=0):
            logL[mass,age], logR[mass,age]= LR.Lumiradius(realmass, Z, fracage)
        else:
            logL[mass,age] = 0
            logR[mass,age] = 0


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
    #I know: MV=V-5*log10(distance)+5-Reddening
    #Mbol=MV+BC
    #L/L_\odot=0.4*(4.72-Mbol)
    #Therefore: V=5*log10(distance)-5+Reddening+4.72-\frac[L}{0.4}-BC
    red=Reddening(distance)
    T=Temperature(logL,logR)
    bc=BC(T)
    if distance==0:
        return 10
    V=5*math.log10(1000*distance)-5+red+4.72-logL/0.4-bc
    return V

graph=np.zeros(20)
graph=graph.astype(float)
masslabel = np.array(range(rangemass))
masslabel = masslabel.astype(float)    
for distance in range(0,rangedistance):
    realdistance=distance*maxdistance/rangedistance
    dV=Volume(distance)*(4./3)*math.pi*dr**3
    for mass in range(0,rangemass):
        realmass = 10**(math.log10(minmass) + mass*dlogM)
        masslabel[mass] = realmass
        mainsequenceage=Mainsequenceage(realmass)
        dM=IMF(realmass)*dlogM
        for age in range(0,rangeage):
            realage = 10**(-1.+age*dlogt)
            fracage = realage/mainsequenceage
            dt=Agerelation(realmass,realage)*dlogt
            if Magnitude(logL[mass,age],logR[mass,age],realdistance)<9 and agerelation[mass,age]==1:
                graphage=int(20*fracage)
                graph[graphage]+=dV*dt*dM

agelabel = np.array(range(20))
agelabel = agelabel.astype(float)
for age in range(0,20):
    agelabel[age]=age/20.

plt.plot(agelabel, graph)
plt.xlabel('t/tms')
plt.ylabel('#stars V<9')
plt.show()
