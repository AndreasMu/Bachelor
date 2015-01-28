import numpy as np
import math
import Lumiradius as LR

Z = 0.02
#The metalicity should be 0.02
numberfactor = 1.
#The numberfactor is what I use to regulate the total number of stars

rangedistance = 10 #The amount of intervals I have for distance
maxdistance = 40.#The maximum distance in kpc
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


def Agerelation(age,mass):
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
            agefactor = Agerelation(ages,masss)
            volume = Volume(distance)
            massfactor = IMF(masss)
            matrix[distance,mass,age] = volume*massfactor*agefactor*numberfactor

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

matrix2=matrix

def Magnitude(logL,distance):
    pass
"""
for distance in range(0,rangedistance):
    for mass in range(0,rangemass):)
        for age in range(0,rangeage):
            matrix2=Magnitude(logL[mass,age],distance)"""
