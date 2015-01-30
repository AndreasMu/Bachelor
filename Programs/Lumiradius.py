import math
import matplotlib.pyplot as plt
import numpy as np

#All those manual input routines are not needed anymore for this program.
"""M=raw_input("Greetings and welcome! My name is Richard and I am Andreas' first python program!\n\n\
Would you be so kind as to give me a number between 0.5 and 50? (This will be the mass of\
 the star we will look at. For the mass of the sun use 1)")
#Give me a Mass
try: #This checks whether the input is actually a number.
    M=float(M)
except:
    print "\n I don't have time for this... I said I need a number! What is that supposed to be?\
 Run me again when you have numbers for me!"
    exit()#End Program if its no number. Have to start over!
if (50>=M) and (M>=0.5):#This will check whether the mass is within the given parameters.
    print"\n Thank you very much!\n"
    X=0
else:
    print"\n Well aren't we clever? Okay, let's roll with it!\n"
    X=1
Z=raw_input("Now give me a number between 0.005 and 0.05\
 This will be our metalicity (The fraction of heavy elements in our star. The sun is at 0.02)")
try:#Same as with the mass I first check whether the input is actually a number.
    Z=float(Z)
except:
    print "Uhm... yeah... that is not a number. Lets try this again, shall we?"
    exit()
if not(0.05>=Z>=0.005):#Checking whether Z is in the given parameters again.
    if(X==0): #Just for fun: Give different answers for different outcomes. Only important thing here is "print"
        print"\n You know... %g is not between 0.005 and 0.05... But okay lets roll with it\n" %Z
        Y=0
        X=1
    elif(X==1):
        print"\n Noooo don't listen to me! Why would you? Okay let's move on...\n"
        Y=0
        X=2
else:
    print"\n We are moving along quite nicely, almost done!\n"
    Y=0
if(X==0):#Again different questions for different values of X. Nothing major.
    Age=raw_input("Awesome! You have been fantastic so far! Now we have to know how far along our star \
is on the main sequence (Basically his normal life). I need a number between 0 and 1.")
    try:#Try whether Age is actually a number.
        Age=float(Age)
    except:
        print"\n And you did so good too... That was not a number. \
 You will have to start over. I'm sorry!"
        exit()#End Program if it's not
    if not(1>=Age>=0):#If the input is not between 0 and 1 i ask for a different age.
        Age=raw_input("\n I'm sorry, but you seem to have typed an invalid age. I need a number between\
 0 and 1. I can't work with anything else unfortunately.")
        try:
            Age=float(Age)
        except:
            print"\n That is not a number! Sorry, but you will have to start over."
            exit()
        if not(1>=Age>=0):#Since an Age outside of zero to one does not make sense, I end the program, if it's not compatible.
            print("\n Wrong Again. I meant it: I can't work with anything else! Please try again. \n \nThis\
 program will now terminate itself\n \n... wait... WHAT?!?")
            exit()
        else:
            print"\n Thank you! Let's get crackin! \n"
    else:
        print"\n Hooray! We have our own little star! Let's get crackin!\n"
elif(X==1):#Again... different outputs for different outcomes. Essentially the same as above.
    Age=raw_input("Okay, we had some fun, but now we have no more room for that. Now we have to know \
how far along our star is on the main sequence (Basically his normal life). I need a number between 0 and 1.\
 and before you get any ideas: I can't work with anything else.")
    try:
        Age=float(Age)
    except:
        print"\n What in the heavens is that? I wanted a Number, not... whatever that is! \n"
        exit()
    if not(1>=Age>=0):
        print"\n Oh... How creative of you... Nope... I'm out"
        exit()
    else:
        print"\n Okay! We have our own... star I guess? Let's hope for the best!\n"
elif(X==2):
    Age=raw_input("... Well let's get this over with. Give me a number between 0 and 1. AND ONLY \
BETWEEN 0 AND 1! I WON'T BE ABLE TO WORK WITH ANYTHING ELSE.")
    try:
        Age=float(Age)
    except: #Someone is mad...
        print"OH COME ON!"
        exit()
    if not(1>=Age>=0):
        K=0
        while(K<20):
            print"...\n"
            K+=1
        print"No."
        exit()
    else:
        print"\n Congratulations it's an abomination! Don't look at me! You gave me these numbers"
"""


#Now that the fun stuff is out of the way, I start calculating.
#First I calculate the coefficients a to subsequently calculate everything else.
#This will be a whole lot of coefficients.

def Lumiradius(M,Z,Age):
    Squiggle= math.log10((Z/0.02))
    a45=0.23214+0.001828075*Squiggle-0.02232007*(Squiggle**2)-0.003378734*(Squiggle**3)
    a46=0.01163659+0.003427682*Squiggle+0.001421393*(Squiggle**2)-0.003710666*(Squiggle**3)
    a47=0.0104802-0.01231921*Squiggle-0.0168686*(Squiggle**2)-0.004234354*(Squiggle**3)
    a48=1.55559-0.3223927*Squiggle-0.5197429*(Squiggle**2)-0.1066441*(Squiggle**3)
    a49=max(0.145,0.0977-0.231*Squiggle-0.0753*(Squiggle**2))
    a50=min(0.306+0.053*Squiggle,0.24+0.18*Squiggle+0.595*(Squiggle**2))
    a51=min(0.3625+0.062*Squiggle,0.33+0.132*Squiggle+0.218*(Squiggle**2))
    if (Z>0.01):
        a52=min(1.0,1.1064+0.415*Squiggle+0.18*(Squiggle**2))
        a53=min(1.1,1.19+0.377*Squiggle+0.176*(Squiggle**2))
    else:
        a52=max(0.9,1.1064+0.415*Squiggle+0.18*(Squiggle**2))
        a53=max(1.0,1.19+0.377*Squiggle+0.176*(Squiggle**2))
    
    A=(a45+a46*(2**a48))/(2**0.4+a47*(2**1.9)) #A is used in one of the aL calculations
    if(M<0.5):
        aL=a49
    elif(0.5<=M<0.7):
        aL=a49+5*(0.3-a49)*(M-0.5)
    elif(0.7<=M<a52):
        aL=0.3+(a50-0.3)*(M-0.7)/(a52-0.7)
    elif(a52<=M<a53):
        aL=a50+(a51-a50)*(M-a52)/(a53-a52)
    elif(a53<=M<2.0):
        aL=a51+(A-a51)*(M-a53)/(2.0-a53)
    elif(M>=2.0):
        aL=(a45+a46*(M**a48))/((M**0.4)+a47*(M**1.9))
    #End of aL calculations. 

    #Begin of bL
    a54=0.3855707-0.6104166*Squiggle+5.676742*(Squiggle**2)+10.60894*(Squiggle**3)+5.284014*(Squiggle**4)
    a55=0.3579064-0.6442936*Squiggle+5.494644*(Squiggle**2)+10.54952*(Squiggle**3)+5.280991*(Squiggle**4)
    a56=0.9587587+0.8777464*Squiggle+0.2017321*(Squiggle**2)
    a57=max(0.6355-0.4192*Squiggle,max(1.25,min(1.4,1.5135+0.3769*Squiggle)))

    bL=max(0.0,a54-a55*(M**a56))
    if(M>a57) and(bL>0.0):
        B=max(0.0,a54-a55*(a57**a56))
        bL=max(0.0,B-10.0*(M-a57)*B)
        
    #End of bL
    
    #Begin of aR    
    a58=0.4907546-0.1683928*Squiggle-0.3108742*(Squiggle**2)-0.07202918*(Squiggle**3)
    a59=4.53707-4.465455*Squiggle-1.61269*(Squiggle**2)-1.623246*(Squiggle**3)
    a60=1.79622+0.281402*Squiggle+1.423325*(Squiggle**2)+0.3421036*(Squiggle**3)
    a61=2.256216+0.37734*Squiggle+1.537867*(Squiggle**2)+0.4396373*(Squiggle**3)
    a62=max(0.065,0.0843-0.0475*Squiggle-0.0352*(Squiggle**2))
    a63=0.0736+0.0749*Squiggle+0.04426*(Squiggle**2)
    if (Z<0.004):
        a63=min(0.055,a63)
    a64=max(0.091,min(0.121,0.136+0.0352*Squiggle))
    a65=0.001564231+0.001653042*Squiggle-0.004439786*(Squiggle**2)-0.004951011*(Squiggle**3)-0.00121653*(Squiggle**4)
    a66=max(1.477+0.296*Squiggle, min(1.6,-0.308-1.046*Squiggle))
    a66=max(0.8,min(0.8-2.0*Squiggle,a66))
    a67=5.210157-4.143695*Squiggle-2.12087*(Squiggle**2)
    a68=min(a66,1.116+0.166*Squiggle)
    
    A2=(a58*(a66**a60))/(a59+(a66**a61))
    A3=(a58*(a67**a60))/(a59+(a67**a61))
    if(M<0.5):
        aR=a62
    elif(0.5<=M<0.65):
        aR=a62+(a63-a62)*(M-0.5)/0.15
    elif(0.65<=M<a68):
        aR=a63+(a64-a63)*(M-0.65)/(a68-0.65)
    elif(a68<=M<a66):
        aR=a64+(A2-a64)*(M-a68)/(a66-a68)
    elif(a66<=M<=a67):
        aR=(a58*(M**a60))/(a59+(M**a61))
    elif(a67<M):
        aR=A3+a65*(M-a67)
    #End of aR
    
    #Begin of bR
    a69=1.071489-0.1164852*Squiggle-0.08623831*(Squiggle**2)-0.01582349*(Squiggle**3)
    a70=0.7108492+0.7935927*Squiggle+0.3926983*(Squiggle**2)+0.03622146*(Squiggle**3)
    a71=3.478514-0.02585474*Squiggle-0.01512955*(Squiggle**2)-0.00283369*(Squiggle**3)
    a72=0.9132108-0.1653695*Squiggle+0.0363784*(Squiggle**3)
    if(Z>0.01):
        a72=max(a72,0.95)
    a73=0.00396933+0.004539076*Squiggle+0.001720906*(Squiggle**2)+0.0001897857*(Squiggle**3)
    a74=max(1.4,min(1.6+0.764*Squiggle+0.3322*(Squiggle**2),1.6))
    
    B2=(a69*(2.0**3.5))/(a70+2.0**a71)
    B3=(a69*(16.0**3.5))/(a70+16.0**a71)
    if(M<=1.0):
        bR=1.06
    elif(1.0<M<a74):
        bR=1.06+(a72-1.06)*(M-1.0)/(a74-1.06)
    elif(a74<=M<2.0):
        bR=a72+(B2-a72)*(M-a74)/(2.0-a74)
    elif(2.0<=M<=16.0):
        bR=(a69*(M**3.5))/(a70+M**a71)
    elif(16<M):
        bR=B3+a73*(M-16.0)
    bR=bR-1    
    #End of bR
    
    #Begin of c
    a75=max(1.0,min(0.8109-0.6282*Squiggle,1.27))
    a75=max(a75,0.6355-0.4192*Squiggle)
    a76=0.01192334+0.01083057*Squiggle+1.230969*(Squiggle**2)+1.551656*(Squiggle**3)
    a76=max(-0.1015564-0.2161264*Squiggle-0.05182516*(Squiggle**2),a76)
    a77=-0.1668868+0.5818123*Squiggle-11.05027*(Squiggle**2)-16.68070*(Squiggle**3)
    a77=max(-0.3868776-0.5457078*Squiggle-0.1463472*(Squiggle**2),min(0.0,a77))
    a78=max(0.0,min(7.454+9.046*Squiggle,0.7615495+0.1068243*Squiggle-0.2011333*(Squiggle**2)-0.09371415*(Squiggle**3)))
    a79=min(max(2.0,-13.3-18.6*Squiggle),9.409838+1.522928*Squiggle)
    a80=max(0.0585542,-0.2711-0.5756*Squiggle-0.0838*(Squiggle**2))
    a81=min(1.5,max(0.4,2.493+1.1475*Squiggle))
    
    C1=a76+a77*(1.0-a78)**a79
    if(a75==1.0):
        C2=C1
    else:
        C2=a80
    if(M<=1.0):
        c=a76+a77*(M-a78)**a79
    elif(1.0<M<=a75):
        c=C1+(a80-C1)*((M-1.0)/(a75-1.0))**a81
    elif(a75<M<a75+0.1):
        c=C2-10.0*(M-a75)*C2
    elif(M>a75+0.1):
        c=0.0
    #End of c
    
    Mh=1.0185+0.16015*Squiggle+0.0892*(Squiggle**2)
    
    #First I just write down all the tier 2 cofficients I need for everything else. Then I calculate tier 1.
    a1=1593.89+2053.038*Squiggle+1231.226*(Squiggle**2)+232.7785*(Squiggle**3)
    a2=2706.708+1483.131*Squiggle+577.2723*(Squiggle**2)+74.1123*(Squiggle**3)
    a3=146.6143-104.8442*Squiggle-67.95374*(Squiggle**2)-13.91127*(Squiggle**3)
    a4=0.0414196+0.04564888*Squiggle+0.02958542*(Squiggle**2)+0.005571483*(Squiggle**3)
    a5=0.3426349
    a6=19.49814+1.758178*Squiggle-6.008212*(Squiggle**2)-4.470533*(Squiggle**3)
    a7=4.90383
    a8=0.05212154+0.03166411*Squiggle-0.002750074*(Squiggle**2)-0.002271549*(Squiggle**3)
    a9=1.312179-0.32944936*Squiggle+0.0923186*(Squiggle**2)+0.02610989*(Squiggle**3)
    a10=0.8073972
    a14=3858.911+2459.681*Squiggle-76.30093*(Squiggle**2)-348.6057*(Squiggle**3)-48.61703*(Squiggle**4)#Needed in a11 so earlier
    a11=a14*(1.031538-0.243448*Squiggle+7.732821*(Squiggle**2)+6.460705*(Squiggle**3)+1.374484*(Squiggle**4))
    a12=a14*(1.043715-1.577474*Squiggle-5.168234*(Squiggle**2)-5.596506*(Squiggle**3)-1.299394*(Squiggle**4))
    a13=785.9573-8.542048*Squiggle-26.42511*(Squiggle**2)-9.585707*(Squiggle**3)
    a15=288.872+295.2979*Squiggle+185.0341*(Squiggle**2)+37.97254*(Squiggle**3)
    a16=7.19658+0.5613746*Squiggle+0.3805871*(Squiggle**2)+0.08398728*(Squiggle**3)
    s=math.log10(Z)#tier 3 coefficient *yay*
    a20=26.52091+81.78458*Squiggle+115.6058*(Squiggle**2)+76.33811*(Squiggle**3)+19.50698*(Squiggle**4)#Needed in a18 so earlier
    a17=math.exp(max(0.097-0.1072*(s+3),max(0.097,min(0.1461,0.1461+0.1237*(s+2)))))
    a18=a20*(0.2187715-2.154437*Squiggle-3.768678*(Squiggle**2)-1.975518*(Squiggle**3)-0.3021475*(Squiggle**4))
    a19=a20*(1.46644+1.839725*Squiggle+6.442199*(Squiggle**2)+4.023635*(Squiggle**3)+0.6957529*(Squiggle**4))
    a21=1.472103-2.947609*Squiggle-3.312828*(Squiggle**2)-0.9945065*(Squiggle**3)
    a22=3.071048-5.679941*Squiggle-9.745523*(Squiggle**2)-3.594543*(Squiggle**3)
    a23=2.61789+1.019135*Squiggle-0.03292551*(Squiggle**2)-0.07445123*(Squiggle**3)
    a24=0.01075567+0.01773287*Squiggle+0.009610479*(Squiggle**2)+0.001732469*(Squiggle**3)
    a25=1.476246+1.899331*Squiggle+1.19501*(Squiggle**2)+0.3035051*(Squiggle**3)
    a26=5.502535-0.06601663*Squiggle+0.09968707*(Squiggle**2)+0.03599801*(Squiggle**3)
    a27=95.11033+68.19618*Squiggle-10.45625*(Squiggle**2)-14.74939*(Squiggle**3)
    a28=31.13458+10.12033*Squiggle-4.650511*(Squiggle**2)-2.463185*(Squiggle**3)
    a32=6.682518+0.2827718*Squiggle-0.07294429*(Squiggle**2) #Needed in a29
    a29=(1.413057+0.4578814*Squiggle-0.06850581*(Squiggle**2)-0.05588658*(Squiggle**3))**a32
    a30=39.10862+51.96646*Squiggle+22.6497*(Squiggle**2)+2.87368*(Squiggle**3)
    a31=4.597479-0.2855179*Squiggle+0.2709724*(Squiggle**2)
    a33=max(0.6355-0.4192*Squiggle,max(1.25,min(1.4,1.5135+0.3769*Squiggle)))
    a34=0.1910302+0.1158624*Squiggle+0.0334899*(Squiggle**2)+0.002599706*(Squiggle**3)
    a35=0.3931056+0.07277637*Squiggle-0.1366593*(Squiggle**2)-0.04508946*(Squiggle**3)
    a36=0.3267776+0.1204424*Squiggle+0.09988332*(Squiggle**2)+0.02455361*(Squiggle**3)
    a37=0.5990212+0.05570264*Squiggle+0.06207626*(Squiggle**2)+0.01777283*(Squiggle**3)
    a38=0.7330122+0.5192827*Squiggle+0.2316416*(Squiggle**2)+0.008346941*(Squiggle**3)
    a39=1.172768-0.1209262*Squiggle-0.1193023*(Squiggle**2)-0.02859837*(Squiggle**3)
    a40=0.3982622-0.2296279*Squiggle-0.2262539*(Squiggle**2)-0.05219837*(Squiggle**3)
    a41=3.571038-0.02223625*Squiggle-0.02611794*(Squiggle**2)-0.006359648*(Squiggle**3)
    a42=min(1.25,max(1.1,1.9848+1.1386*Squiggle+0.3564*(Squiggle**2)))
    a43=0.063+0.0481*Squiggle+0.00984*(Squiggle**2)
    a44=min(1.3,max(0.45,1.2+2.45*Squiggle))
    
    #Calculation of a few different tier 1 coefficients.#LZAMS M**3 TERM DIFFERENT FROM PAPER!
    LZAMS=(0.3970417*M**5.5+8.527626*M**11)/(0.00025546*M**3+5.432889*M**5+5.563579*M**7+0.7886606*M**8+0.00586685*M**9.5)
    RZAMS=(1.715359*M**2.5+6.597788*M**6.5+10.08855*M**11+1.012495*M**19+0.07490166*M**19.5)/(0.01077422+3.082234*M**2+17.84778*M**8.5+M**18.5+0.00022582*M**19.5)
        
    dA=min(a34/(a33**a35),a36/(a33**a37)) #Tier 2
    if(M<=Mh):
        dL=0.0
    elif(Mh<M<a33):
        dL=dA*(((M-Mh)/(a33-Mh))**0.4)
    elif(M>=a33):
        dL=min(a34/(M**a35),a36/(M**a37))
    
    dB=(a38+a39*(2.0**3.5))/(a40*(2.0**3)+(2.0**a41))-1.0#Tier 2
    if(M<Mh):
        dR=0.0
    elif(Mh<M<=a42):
        dR=a43*((M-Mh)/(a42-Mh))**0.5
    elif(a42<M<2.0):
        dR=a43+(dB-a43)*((M-a42)/(2.0-a42))**a44
    elif(M>=2.0):
        dR=(a38+a39*(M**3.5))/(a40*(M**3)+(M**a41))-1.0
    
    LTMS=(a11*M**3+a12*M**4+a13*M**(a16+1.8))/(a14+a15*M**5+M**a16)
    
    if(M<=a17):
        RTMS=(a18+a19*M**a21)/(a20+M**a22)
    elif(M>=a17+0.1):
        RTMS=(-0.08672073*M**3+a23*M**a26+a24*M**(a26+1.5))/(a25+M**5)
    else: #I interpolate linearly by multiplying every term with 10*(M-a17)(0 at M=a17, 1 at M=a17+.1) and (a17+.1-M)(0 at M=a17+0.1 and 1 at M=a17)
        RTMS=(10*(a17+0.1-M)*(a18+a19*M**a21)/(a20+M**a22)+(M-a17)*10*(-0.08672073*M**3+a23*M**a26+a24*M**(a26+1.5))/(a25+M**5))/(2.0)
    if(M<0.5):
        RTMS=max(RTMS, 1.5*RZAMS)

    k=max(0.95,min(0.95-0.03*(Squiggle+0.30103),0.99))
    mu=max(0.5,1.0-0.01*max(a6/M**a7,a8+a9/M**a10))
    tbgb=(a1+a2*M**4+a3*M**5.5+M**7)/(a4*M**2+a5*M**7)
    thook=mu*tbgb
    tms=max(thook,k*tbgb)
    
    t1=min(1.0,tms*Age/thook)
    t2=max(0.0,min(1.0,(tms*Age-(1.0-0.01)*thook)/(0.01*thook)))
    
    if(Z<=0.0009):
        if(M<=1.0):
            eta=10
        if(M>=1.1):
            eta=20
        else:
            eta=10+(M-1.0)*100 #This is a linear interpolation between both instances of eta.
    else:
        eta=10
    

    
    #This is now the actual calculation of L and R. 
    logL=aL*Age+bL*Age**eta+(math.log10(LTMS/LZAMS)-aL-bL)*Age**2-dL*(t1**2-t2**2)+math.log10(LZAMS)
    logR=aR*Age+bR*Age**10+c*Age**40+(math.log10(RTMS/RZAMS)-aR-bR-c)*Age**3-dR*(t1**3-t2**3)+math.log10(RZAMS)

    #Here to test the function
    #arrayR[Iteration]=logR
    #arrayL[Iteration]=logL
    #arrayT[Iteration]=Age
    #print"LTMS=%g"%LTMS
    #End of calculations!
    return (logL, logR)

#Testing whether this function actually works
"""arrayR=np.array(range(101))
arrayR=arrayR.astype(float)
arrayL=np.array(range(101))
arrayL=arrayL.astype(float)
arrayT=np.array(range(101))
arrayT=arrayT.astype(float)
Age=0
Iteration=0
M=1.
Z=0.02
while(Iteration<101):
    Lumiradius(M,Z,Age)
    Age+=0.01
    Iteration+=1
plt.plot(arrayT, arrayR)
plt.xlabel('t/tms')
plt.ylabel('log(R)')
plt.show()"""
