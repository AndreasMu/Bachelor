BolomCorr=0	;0=Balona(1994), 1=Flower (1996), 2=Nieva (2013) 15800<Teff<32000 K

;Type the magnitude in V band
V=4.235
;Type the temperature in K
T=12800.
;Type the distance in pc
D=120.5
;Type the reddening E(B-V)
R=0.00
;Type the bolometric correction
BC=0.0


;D=1/(D/1000.)
MV=V-5.*alog10(D)+5.-3.1*R
th=5040./T
if BC eq 0.0 then begin
    if BolomCorr eq 0 then BC=-5.5647+18.9446*th-19.8827*(th^2.)+6.1302*(th^3.)
    if BolomCorr eq 1 then begin
	print,'Teff: ',alog10(T)
	if alog10(T) le 3.7 then begin
		BC=-210.13793+0.19596489*T-7.4465325e-05*T^2+1.4337726e-08*T^3- $
			1.3955426e-12*T^4+5.4925758e-17*T^5
	endif
	if alog10(T) ge 4.0 then begin
		BC=4.1940953-0.00070441042*T+3.4516521e-08*T^2-9.5565244e-13*T^3+ $
			1.2790825e-17*T^4-6.4741275e-23*T^5
	endif
	if alog10(T) lt 4.0 and alog10(T) gt 3.7 then begin
		BC=-29.325541+0.018052720*T-4.4823439e-06*T^2.+5.5894085e-10*T^3.- $
			3.4753865e-14*T^4.+8.5372998e-19*T^5
	endif
    endif
    if BolomCorr eq 2 then BC=21.00-5.34*alog10(T)
endif
Mbol=MV+BC
L=0.4*(4.72-Mbol)
Lbol=(10^((4.72-Mbol)/2.5))*!L_sun

print,'Absolute magnitude:',MV
print,'Bolometric correction:',BC
print,'Bolometric magnitude:',Mbol
print,'Luminosity (log(L/L0)):',L
print,'Luminosity (Lbol):',Lbol

end
