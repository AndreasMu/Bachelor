import numpy as np
import math
from scipy.constants import sigma
import Lumiradius as LR
import matplotlib.pyplot as plt

Z = 0.02
#The metalicity should be 0.02

rangedistance = 2 #The amount of intervals I have for distance
maxdistance = 3.#The maximum distance in kpc
dr = float(maxdistance)/rangedistance #size of distance bins

rangemass = 100 #Amount of intervals I have for mass
maxmass = 50. #Maximum and minimum mass
minmass = 5.
dlogM = (math.log10(maxmass)-math.log10(minmass))/rangemass #Size of logarithmic mass bins

rangeage = 100

#Amount of intervals I have for age


#This function is an excerpt from Lumiradius.py. Since Z=0.02 in our case the formulas
#get rather simple in comparison. If Z is changed in further editing, this formula needs to
#be revised.
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

maxage = Mainsequenceage(minmass) #The star with the minimal mass is also the one, that lives the longest
dlogt=(math.log10(maxage)+2)/rangeage #Size of logarithmic time bins

#Probability density in distance. Right now the stars are homogenously distributed in a 3kpc Sphere.
def DpdV(distance):
    dpdv=1/((4./3)*math.pi*maxdistance**3)
    return dpdv

"""volumelabel = np.array(range(rangedistance))
volumelabel = volumelabel.astype(float)
for distance in range(0,rangedistance):
    volumelabel[distance]=Volume(distance)
    #The volumelabel will save the fractional volume of each shell."""

#Probability density in Mass. This is essentially the IMF converted for a logarithmic scale
def DpdlogM(realmass):
    A= 1.35 * (1/(minmass**-1.35-maxmass**-1.35))
    dpdlogM =  math.log(10)*A*realmass**(-1.35)
    return (dpdlogM)

#Probability density in age. Stars, that exceed their lifetime will not be counted.
def Dpdlogt(realmass,realage):
    tms = Mainsequenceage(realmass)
    if(realage>tms):
        dpdlogt = 0
    else:
        dpdlogt = math.log(10)*realage/tms
    return dpdlogt

#Here I make two matrices, which will contain Luminosity and Radius of stars of
#a specific age and mass. In the later parts of this program I can then simply
#use these matrices instead of always having to call the function from Lumiradius.py
logL = np.arange(rangemass*rangeage,dtype=np.float)
logR = np.arange(rangemass*rangeage,dtype=np.float)
logL = logL.reshape((rangemass,rangeage))
logR = logR.reshape((rangemass,rangeage))
for mass in range(0,rangemass):
    realmass = 10**(math.log10(minmass) + mass*dlogM)
    for age in range(0,rangeage):
        realage = 10**(-2.+age*dlogt)
        fracage = realage/Mainsequenceage(realmass)
        if (fracage<=1):
            logL[mass,age], logR[mass,age]= LR.Lumiradius(realmass, Z, fracage)
        else:
            logL[mass,age] = 0
            logR[mass,age] = 0


#This is the start of calculating the visual magnitude of the star. All of these functions are taken from
#the program: "luminosity.pro" and implemented here.
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
    V=5*math.log10(1000*distance)-5+red+4.72-logL/0.4-bc
    return V

dtau=1./20. #size of bins in tau (fractional main sequence age)
graph=np.zeros(20) #These two arrays will be filled with juicy data.
graph=graph.astype(float)
nocut=np.zeros(20)
nocut=nocut.astype(float)
for distance in range(1,rangedistance+1): #If distance starts at 0, there will be a problem with the logarithm in the magnitude function.
    realdistance=distance*maxdistance/rangedistance
    dV=4*math.pi*dr**3*distance**2 #The approximated volume of the nth shell
    dpv=DpdV(distance)*dV #Probability to find a star in this volume segment
    for mass in range(0,rangemass):
        realmass = 10**(math.log10(minmass) + mass*dlogM)
        mainsequenceage=Mainsequenceage(realmass)
        dpm=DpdlogM(realmass)*dlogM #Probability to find a star in this mass segment
        for age in range(0,rangeage):
            realage = 10**(-2+age*dlogt)
            fracage = realage/mainsequenceage #Fractional mainsequenceage
            graphage=int(20*fracage) #This is the bin in the array, that corresponds to fracage
            dpt=Dpdlogt(realmass,realage)*dlogt #Probability to find a star in this age segment
            if fracage<=1: #For the nocut array I just don't use the magnitude cut
                nocut[graphage]+=dpv*dpt*dpm/dtau
            if Magnitude(logL[mass,age],logR[mass,age],realdistance)<9 and fracage<=1: #Exactly the same as nocut, but only with stars of mag <9
                graph[graphage]+=dpv*dpt*dpm/dtau
                

agelabel = np.arange(20,dtype=np.float)#This is tau 
for age in range(0,20):
    agelabel[age]=age/20.

norm=0. #These loops normalize graph and nocut.
normnocut=0.
for bins in range(0,20):
    norm+=graph[bins]*dtau
    normnocut+=nocut[bins]*dtau
for bins in range(0,20):
    graph[bins]=graph[bins]/norm
    nocut[bins]=nocut[bins]/normnocut

#Now I write the data I get into a file for further handling.
Dataforfile=np.column_stack((agelabel,graph))
np.savetxt('compdata.txt',Dataforfile,fmt='%5.3f')

Dataforfilenocut=np.column_stack((agelabel,nocut))
np.savetxt('compdatanocut.txt',Dataforfilenocut,fmt='%5.3f')

plt.plot(agelabel, graph)
#plt.errorbar(agelabel, data, yerr=error,fmt='')
plt.xlabel('t/tms')
plt.ylabel('probability density')
plt.show()
