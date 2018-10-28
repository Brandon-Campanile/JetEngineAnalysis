Rbar = 8314
Ta = 220
Pa = 10*10**3
Patm = 101325
M = 1.50
PRc = 30
PRf = 1.2
B = 2
b = 0.1
f = 0.018
fab = 0.010
nadd = 0.92
yd = 1.4
MW = 28.8
yf = 1.4
npf = 0.90
npc = 0.90
yc = 1.38
nb = 0.99
yb = 1.33
delH = 45*10**6
PRb = 0.98
Tmaxo = 1300
Cb1 = 700
bmax = 0.12
yt = 1.33
npt = 0.92
delp = 550*10**3
rhof = 780
nadp = 0.35
ytm = 1.34
npft = 0.92
yft = 1.33
nab = 0.96
yab = 1.32
PRab = 0.97
nadn = 0.95
yn = 1.35
nadfn = 0.97
yfn = 1.4
PRnm = 0.8
nadcn = 0.95
ycn = 1.37
u = M*(1.4*Rbar*Ta/MW)**0.5
CB1 = 245
deld = CB1*(M**2)*(Pa/Patm)*(B**1.5)

def diffuser(Ta, Pa, M, nad, y):
    if M<=1:
        rd = 1
    elif M>1 and M<5:
        rd = 1-0.075*(M-1)**1.35
    Po1 = Pa*rd*(1+nad*0.5*(y-1)*M**2)**(y/(y-1))
    To1 = Ta*(1+0.5*(y-1)*M**2)
    return Po1, To1

def fan(To1, Po1, npf, y, MW, PRf):
    cpf = y*(Rbar/MW)/(y-1)
    To2 = To1*PRf**((y-1)/(y*npf))
    Po2 = Po1*PRf
    Wf_ma = (1+B)*cpf*(To2-To1)
    return Po2, To2, Wf_ma

def compressor(To2, Po2, PRc, npc, y, MW):
    cpc = y*(Rbar/MW)/(y-1)
    To3 = To2*PRc**((y-1)/(y*npc))
    Po3 = Po2*PRc
    Wc_ma = cpc*(To3-To2)
    return Po3, To3, Wc_ma

def burner(To3, Po3, PRb, y, MW, nb, delH, Tmaxo, Cb1, b, bmax):
    cpb = y*(Rbar/MW)/(y-1)
    To4 = (1/(1-b+f))*((1-b)*To3+nb*delH*f/cpb)
    Po4 = Po3*PRb
    Tmax = Tmaxo +Cb1*(b/bmax)**0.5
    #fmax = (1-To3/To4)/(delH*nb/cpb/To4-1) #needs improvement
    return Po4, To4, Tmax

def fuel_pump(delp, rhof, nadp, f):
    Wp_ma = f*delp/rhof/nadp
    return Wp_ma

def turbine(To4, Po4, Wc, Wp, y, npt, MW):
    cpt = y*(Rbar/MW)/(y-1)
    To5_1 = To4 - Wc/(cpt*(1+f-b))-Wp/(cpt*(1+f-b))
    TR = To5_1/To4
    Po5_1 = Po4*(TR**(1/npt))**(y/(y-1))
    return Po5_1, To5_1

def tmixer(To3, To5_1, Po5_1, f, b, y):
    To5_m = (b*To3+(1+f-b)*To5_1)/(1+f)
    #Po5_m = Po5_1
    #Po5_m = Po5_1*((To5_m/To5_1)**(y/(y-1)))*((To5_1/To3)**(y*b/((y-1)*(1+f-b))))
    Po5_m = Po5_1*((To5_m/To5_1)**(y/(y-1)))*((To5_1/To3)**(y*b/((y-1)*(1+f))))
    return Po5_m, To5_m

def fan_turbine(To5_m, Po5_m, npft, y, Wf, f, MW):
    cpft = y*(Rbar/MW)/(y-1)
    To5_2 = To5_m - Wf/(cpft*(1+f))
    TR = To5_2/To5_m
    Po5_2 = Po5_m*(TR**(1/npft))**(y/(y-1))
    return Po5_2, To5_2

def afterburner(To5_2, Po5_2, PRab, y, MW, fab, f, delH, nab):
    cpab = y*(Rbar/MW)/(y-1)
    if fab>0:
        PR = PRab
    else:
        PR = 1
    Po6 = Po5_2*PR
    To6 = (1/(1+f+fab))*((1+f)*To5_2+nab*delH*fab/cpab)
    return Po6, To6

def nozzle(To6, Po6, Pa, nadn, y, MW):
    cpn = y*(Rbar/MW)/(y-1)
    Te = To6*(1-nadn*(1-(Pa/Po6)**((y-1)/y)))
    Pe = Pa
    ue = (2*cpn*(To6-Te))**0.5
    return Pe, Te, ue

def fan_nozzle(To2, Po2, Pa, nadfn, y, MW):
    cpfn = y*(Rbar/MW)/(y-1)
    Tef = To2*(1-nadfn*(1-(Pa/Po2)**((y-1)/y)))
    Pef = Pa
    uef = (2*cpfn*(To2-Tef))**0.5
    return Pef, Tef, uef

def nozzle_mixer(To6, Po6, To2, Po2, PRnm, B, f, fab):
    To7 = (B*To2+(1+f+fab)*To6)/(1+B+f+fab)
    y = 1.44-(1.39*10**-4)*To7+(3.57*10**-8)*To7**2
    Po7 = Po6*PRnm*((Po2/Po6)**(B/(1+B+f+fab)))*((To7/To6)**(y/(y-1)))*((To6/To2)**(y*B/((y-1)*(1+B+f+fab))))
    return Po7, To7, y

def combined_nozzle(To7, Po7, Pa, y, MW, nadcn):
    cpcn = y*(Rbar/MW)/(y-1)
    Tec = To7*(1-nadcn*(1-(Pa/Po7)**((y-1)/y)))
    Pec = Pa
    uec = (2*cpcn*(To7-Tec))**0.5
    return Pec, Tec, uec

def st(ue, uef, u, f, fab, B, deld):
    ST = (1+f+fab)*ue+B*uef-(1+B)*u-deld
    return ST

def stcn(uec, u, f, fab, B, deld):
    STCN = (1+B+f+fab)*uec-(1+B)*u-deld
    return STCN

def tsfc(ST, f, fab):
    TSFC = (f+fab)/ST
    return TSFC

def prop_eff(ST, ue, uef, u, f, fab, B):
    KE = (1+f+fab)*0.5*ue**2+B*0.5*uef**2-(1+B)*0.5*u**2
    np = ST*u/KE
    return np, KE

def prop_eff_cn(STCN, uec, u, f, fab, B):
    KECN = (1+B+f+fab)*0.5*uec**2-(1+B)*0.5*u**2
    npcn = STCN*u/KECN
    return npcn, KECN

def therm_eff(KE, f, fab, delH):
    nth = KE/((f+fab)*delH)
    return nth

Po1, To1 = diffuser(Ta, Pa, M, nadd, yd)
Po2, To2, Wf_ma = fan(To1, Po1, npf, yf, MW, PRf)
Po3, To3, Wc_ma = compressor(To2, Po2, PRc, npc, yc, MW)
Po4, To4, Tmax = burner(To3, Po3, PRb, yb, MW, nb, delH, Tmaxo, Cb1, b, bmax)
Wp_ma = fuel_pump(delp, rhof, nadp, f)
Po5_1, To5_1 = turbine(To4 , Po4, Wc_ma, Wp_ma, yt, npt, MW)
Po5_m, To5_m = tmixer(To3, To5_1, Po5_1, f, b, ytm)
Po5_2, To5_2 = fan_turbine(To5_m, Po5_m, npft, yft, Wf_ma, f, MW)
Po6, To6 = afterburner(To5_2, Po5_2, PRab, yab, MW, fab, f, delH, nab)
Pe, Te, ue = nozzle(To6, Po6, Pa, nadn, yn, MW)
Pef, Tef, uef = fan_nozzle(To2, Po2, Pa, nadfn, yfn, MW)
Po7, To7, ynm = nozzle_mixer(To6, Po6, To2, Po2, PRnm, B, f, fab)
Pec, Tec, uec = combined_nozzle(To7, Po7, Pa, ycn, MW, nadcn)
ST = st(ue, uef, u, f, fab, B, deld)
STCN = stcn(uec, u, f, fab, B, deld)
TSFC = tsfc(ST, f, fab)
TSFCCN = tsfc(STCN, f, fab)
np, KE = prop_eff(ST, ue, uef, u, f, fab, B)
npcn, KECN = prop_eff_cn(STCN, uec, u, f, fab, B)
nth = therm_eff(KE, f, fab, delH)
nthcn = therm_eff(KECN, f, fab, delH)
no = np*nth
nocn = npcn*nthcn