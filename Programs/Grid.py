import numpy as np
import math
from scipy.constants import sigma
import Lumiradius as LR

Z = 0.02
#The metalicity should be 0.02
numberfactor = 100000.
#The numberfactor is what I use to regulate the total number of stars

rangedistance = 10 #The amount of intervals I have for distance
maxdistance = 5.#The maximum distance in kpc
stepdistance = float(maxdistance)/rangedistance

rangemass = 10 #Amount of intervals I have for mass
maxmass = 9.
minmass = 1.

rangeage = 10 #Amount of intervals I have for age

rangeall = rangedistance*rangemass*rangeage
#Total amount of different possible stars.

array = np.array(range(rangeall))
matrix = array.reshape((rangedistance,rangemass,rangeage))
#First dimension: distance
#Second dimension: mass
#Third dimension: age

distancelabel = np.array(range(rangedistance))
distancelabel = distancelabel.astype(float)
masslabel = np.array(range(rangemass))
masslabel = masslabel.astype(float)
agelabel = np.array(range(rangeage))
agelabel = agelabel.astype(float)
volumelabel = distancelabel
#Changing every array from an int array to a float array.

#This function is an excerpt from Lumiradius.py. Since Z=0.02 in our case the formulas
#get rather simple in comparison.
#tms is in units of Myr
def Mainsequenceage(mass):
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
    mu = max(0.5,1.0-0.01*max(a6/mass**a7,a8+a9/mass**a10))
    tbgb = (a1+a2*mass**4+a3*mass**5.5+mass**7)/(a4*mass**2+a5*mass**7)
    thook = mu*tbgb
    tms = max(thook,k*tbgb)
    return tms

maxage = Mainsequenceage(minmass)

def Volume(distance):
    distance+=1.
    volumenorm=1/((4./3)*math.pi*maxdistance**3)
    volume = volumenorm*(distance**3-(distance-1)**3)*(4./3)*math.pi*stepdistance**3
    #volume = 4*math.pi*stepdistance**3*distance**3
    #The Volume of a spherical shell is 4/3*pi*(a^3-b^3) a>b. So this is the
    #volume of the nth shell.
    return volume

for distance in range(0,rangedistance):
    volumelabel[distance]=Volume(distance)
    #The volumelabel will save the fractional volume of each shell.

def IMF(mass):
    #The IMF is: IMF=e*mass**{-2.3} To get e i calculate the integral IMF from minmass to 100.
    #1==int_{minmass}^{100}(e*m^{-2.3}}=\frac{-1}{1.3}*e*100^{-1.3} + \frac{1}{1.3}*e*minmass^{-1.3}
    #because minmass^{-1.3}>>100^{-1.3} I neglect the first term. Solving for e gives:
    e= 1.35*minmass**1.35
    massfactor = e*mass**(-2.35)
    return massfactor


def Agerelation(mass,age):
    #This functions purpose is to make sure I don't include stars, that are already dead.
    #If the age of the star is greater, than its mainsequence age, the function will return 0.
    tms = Mainsequenceage(mass)
    age = maxage*age
    if(age>tms):
        agefactor = 0
    else:
        agefactor = 1
    return agefactor

for distance in range(0,rangedistance):
    for mass in range(0,rangemass):
        masss = 10**((float(mass)*2)/rangemass) #This will make masss logarithmic
        masslabel[mass] = masss
        #print"%g"%masss
        for age in range(0,rangeage):
            ages = float(age)/(rangeage)
            agelabel[age] = ages
            agefactor = Agerelation(masss,ages)
            volume = Volume(distance)
            massfactor = IMF(masss)
            matrix[distance,mass,age] = volume * massfactor * agefactor * numberfactor

#Here I make two matrices, which will contain Luminosity and Radius of stars of
#a specific age and mass. In the later parts of this program I can then simply
#use these matrices instead of always having to call the function from Lumiradius.py
logK = np.array(range(rangemass*rangeage))
logL = logK.reshape((rangemass,rangeage))
logR = logK.reshape((rangemass,rangeage))
logL = logL.astype(float)
logR = logR.astype(float)

for mass in range(0,rangemass):
    masss = 10**((float(mass)*2)/rangemass)
    for age in range(0,rangeage):
        ages = float(age)/(rangeage)
        if (matrix[rangedistance-1,mass,age]!=0):
            logL[mass,age], logR[mass,age]= LR.Lumiradius(masss, Z, ages)
        else:
            logL[mass,age] = 0
            logR[mass,age] = 0


magarray = np.array(range(rangeall))
magnitude = magarray.reshape((rangedistance,rangemass,rangeage))
magnitude = magnitude.astype(float)

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
    #Therefore: V=5*log10(distance)-5+Reddening+\frac{4.72*0.4-L}{0.4}-BC
    red=Reddening(distance)
    T=Temperature(logL,logR)
    bc=BC(T)
    if distance==0:
        distance=0.01*maxdistance/rangedistance
    V=5*math.log10(distance)-5+red+(4.72-logL/0.4)-bc
    return V
    
for distance in range(0,rangedistance):
    distances = distance*maxdistance/rangedistance
    for mass in range(0,rangemass):
        masss = 10**((float(mass)*2)/rangemass)
        for age in range(0,rangeage):
            ages = float(age)/(rangeage)
            if (Agerelation(masss,ages)!=0):
                magnitude[distance,mass,age]=Magnitude(logL[mass,age],logR[mass,age],distances)
            else:
                magnitude[distance,mass,age]=0
                
