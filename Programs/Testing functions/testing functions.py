#This will tell me how many entries there are in the agerelation matrix so
#essentially how many stars I can simulate of a given mass.
agerange = np.zeros(rangemass)
agerange = agerange.astype(int)
agerangee=0
for mass in range(0,rangemass):
    agerangee=0
    for age in range(0,rangeage):
        agerangee+=1
        if agerelation[mass,age]==0:
            agerange[mass]=agerangee-1
            break
        if age==rangeage-1:
            agerange[mass]=agerangee
