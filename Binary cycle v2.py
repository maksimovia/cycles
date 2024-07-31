import numpy as np
from CoolProp.CoolProp import PropsSI as prop
from data import nodes, blocks
from modules import comp, turb, heat_exch_2streams,comb_stoic,mix,heat_exch
from scipy.optimize import root_scalar
P3g = 3e6
dP_KU = 0.025
P10g = 1e5
P6g = P10g*(1+dP_KU)
T6g = 1060+273.15
KPDcomp = 0.8
KPDturb = 0.9

Pk = 0.004246971e6
P1 = 4.138998301e6
Pd = 0.12e6
dP_deair = 0.3
#dT_deair = 5
T8 = 105 + 273.15
dT_econ = 10
dP_econ = 0.25e6
T7 = 60 + 273.15
T1 = 515 + 273.15
dT_pinch = 15
dP_PE = 0.25e6
dP_GPK = 0.025


X1g = 'REFPROP::O2[0.20]&N2[0.77]&CO2[0.01]&H2O[0.01]&Ar[0.01]'
T1g = 15+273.15
P1g = 1e5
G1g = 120
H1g = prop("H", "P", P1g, "T", T1g, X1g)
S1g = prop("S", "P", P1g, "T", T1g, X1g)
Q1g = prop("Q", "P", P1g, "T", T1g, X1g)
nodes.loc['1g'] = [T1g, P1g, H1g, S1g, Q1g, G1g, X1g]
X2g = 'REFPROP::Methane[1]&H2[0]&CO[0]'
T2g = 50+273.15
P2g = 1.2e6
H2g = prop("H", "P", P2g, "T", T2g, X2g)
S2g = prop("S", "P", P2g, "T", T2g, X2g)
Q2g = prop("Q", "P", P2g, "T", T2g, X2g)
nodes.loc['2g'] = [T2g, P2g, H2g, S2g, Q2g, '', X2g]
def Troot(G2g):
    nodes.loc['2g','G'] = float(G2g)
    comp('AirCOMP', '1g', '3g', P3g, KPDcomp)
    comp('FuelCOMP', '2g', '4g', P3g, KPDcomp)
    comb_stoic('COMB', '3g', '4g','5g',dP=0)
    return nodes.loc['5g']['T'] - T6g
root_scalar(Troot,x0=1, xtol=10**-5)
turb('GTURB', '5g', '6g', P6g, KPDturb)

D0_ = input[0]
Gotb_ = input[1]
Grec_ = input[2]
T8g_ = input[3]
T7g_ = input[4]

def Root(input):
    D0_ = input[0]
    Gotb_ = input[1]
    Grec_ = input[2]
    T8g_ = input[3]
    T7g_ = input[4]

    X1 = 'H2O'#1
    T1 = T1
    P1 = P1
    G1 = D0_#11111111
    H1 = prop("H", "P", P1, "T", T1, X2g)
    S1 = prop("S", "P", P1, "T", T1, X2g)
    Q1 = prop("Q", "P", P1, "T", T1, X2g)
    nodes.loc['1'] = [T1, P1, H1, S1, Q1, G1, X1]
    turb('CVD','1','2',Pd*(1+dP_deair),KPDturb)#2
    nodes.loc['3'] = nodes.loc['2']#3
    nodes.loc['3','G'] = nodes.loc['2']['G'] - Gotb_ #2222222
    turb('CND','3','4',Pk,KPDturb)#4
    heat_exch('COND','4','5',x=0)#5
    comp('CONPUMP','5','6',Pd*(1+dP_GPK),KPDcomp)#6
    X8 = 'H2O'#8
    T8 = T8
    P8 = Pd
    G8 = nodes.loc['6']['G'] + Grec_#3333333
    H8 = prop("H", "P", P8, "T", T8, X8)
    S8 = prop("S", "P", P8, "T", T8, X8)
    Q8 = prop("Q", "P", P8, "T", T8, X8)
    nodes.loc['8'] = [T8, P8, H8, S8, Q8, G8, X8]
    nodes.loc['9'] = nodes.loc['8']#9
    nodes.loc['9','G'] = Grec_#3333333
    comp('RECPUMP', '9', '10', Pd * (1 + dP_GPK), KPDcomp)#10
    mix('REC', '6', '10', '7')#7
    nodes.loc['11'] = nodes.loc['8']#11
    nodes.loc['11','G'] = nodes.loc['6']['G']
    nodes.loc['12'] = nodes.loc['2']  # 12
    nodes.loc['12','G'] = Gotb_#22222222
    X13 = 'H2O'#13
    H13 = nodes.loc['2','H']
    P13 = Pd
    G13 = nodes.loc['12','G']#222222
    T13 = prop("T", "P", P13, "H", H13, X13)
    S13 = prop("S", "P", P13, "H", H13, X13)
    Q13 = prop("Q", "P", P13, "H", H13, X13)
    nodes.loc['13'] = [T13, P13, H13, S13, Q13, G13, X13]
    mix('DEAIR','13','11','14')#14
    comp('PPUMP','14','15',P1+dP_econ+dP_PE,KPDcomp)#15
    X8g = nodes.loc['6g','fluid']#8g
    T8g = T8g_ #44444
    P8g = P10g*(1+dP_KU*2/4)
    G8g = nodes.loc['6g','G']
    H8g = prop("H", "P", P8g, "T", T8g, X8g)
    S8g = prop("S", "P", P8g, "T", T8g, X8g)
    Q8g = prop("Q", "P", P8g, "T", T8g, X8g)
    nodes.loc['8g'] = [T8g, P8g, H8g, S8g, Q8g, G8g, X8g]
    X17 = 'H2O' # 17
    Q17 = 0
    P17 = nodes.loc['15','P'] - dP_econ
    G17 = D0_ #111111
    H17 = prop("H", "P", P17, "Q", Q17, X17)
    S17 = prop("S", "P", P17, "Q", Q17, X17)
    T17 = prop("T", "P", P17, "Q", Q17, X17)
    nodes.loc['17'] = [T17, P17, H17, S17, Q17, G17, X17]
    heat_exch_2streams('ECON','8g','9g','15','16',T22=T17-dT_econ)#16
    X7g = nodes.loc['6g','fluid']#7g
    T7g = T7g_ #55555
    P7g = P10g*(1+dP_KU*3/4)
    G7g = nodes.loc['6g','G']
    H7g = prop("H", "P", P7g, "T", T7g, X7g)
    S7g = prop("S", "P", P7g, "T", T7g, X7g)
    Q7g = prop("Q", "P", P7g, "T", T7g, X7g)
    nodes.loc['7g'] = [T7g, P7g, H7g, S7g, Q7g, G7g, X7g]
    heat_exch_2streams('EVAP','7g','8g','17','18',Q22=1)#18
    nodes.loc['19'] = nodes.loc['18']#19
    heat_exch_2streams('PP','6g','7g','19','1',T22=T1)

    eq1 =
    eq2 =
    eq3 =
    eq4 =
    eq5 = nodes.g_

#111
#222
#333
#444
#555

print(nodes)