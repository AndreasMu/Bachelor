import numpy as np
import math

Z=0.02
#The metalicity should be 0.02
numberfactor = 1.
#The numberfactor is what I use to regulate the total number of stars

rangedistance=10 #The amount of intervals I have for distance
maxdistance=40.#The maximum distance in kpc
stepdistance=float(maxdistance)/rangedistance

rangemass=10 #Amount of intervals I have for mass
maxmass=10.
minmass=1.

rangeage=10 #Amount of intervals I have for age

rangeall=rangedistance*rangemass*rangeage

array= np.array(range(rangeall))
matrix=array.reshape((rangedistance,rangemass,rangeage))
#First dimension: distance
#Second dimension: mass
#Third dimension: age

#This function is an excerpt from test.py. Since Z=0.02 in our case the formulas
#get rather simple in comparison.
def Mainsequenceage(mass):
    a1=1593.89
    a2=2706.708
    a3=146.6143
    a4=0.0414196
    a5=0.3426349
    a6=19.49814
    a7=4.90383
    a8=0.05212154
    a9=1.312179
    a10=0.8073972
    k=max(0.95,min(0.95-0.03*0.30103,0.99))
    mu=max(0.5,1.0-0.01*max(a6/mass**a7,a8+a9/mass**a10))
    tbgb=(a1+a2*mass**4+a3*mass**5.5+mass**7)/(a4*mass**2+a5*mass**7)
    thook=mu*tbgb
    tms=max(thook,k*tbgb)
    return tms

maxage=Mainsequenceage(1.)

def Distancerelation(distance):
    volume=(distance**3-(distance-1)**3)*(4./3)*math.pi*stepdistance**3
    #The Volume of a spherical shell is 4/3*pi*(a^3-b^3) a>b. So this is the
    #volume of the nth shell.
    return volume

def Massrelation(mass):
    #integral IMF from minmass to infty =e*l**-1.35/1.35. norm to 1 => 
    e=1.35*minmass**1.35
    massfactor = e*mass**(-2.35)
    return massfactor

def Agerelation(age,mass):
    tms=Mainsequenceage(mass)
    age=maxage*age
    if(age>tms):
        agefactor = 0
    else:
        agefactor = 1
    return agefactor

#This is just to check the Distance and Massrelation functions.
"""maxvolume=0
for distance in range(0,rangedistance):
    volume=Distancerelation(distance)
    print"volume= %g"%volume
    maxvolume+=volume
    
for mass in range(0,rangemass):
    massfactor=Massrelation(mass)
"""

for distance in range(1,rangedistance+1):
    for mass in range(0,rangemass):
        masss = (float(mass)/rangemass)*maxmass
        masss = (masss**8/maxmass**8)*maxmass +1
        print"%g"%masss
        #Irgendwie muss ich das logarithmisch kriegen...
        for age in range(1,rangeage+1):
            ages=float(age)/rangeage
            agefactor = Agerelation(ages,masss)
            volume = Distancerelation(distance)
            massfactor = Massrelation(masss)
            matrix[distance-1,mass,age-1]=volume*massfactor*agefactor*numberfactor







            
