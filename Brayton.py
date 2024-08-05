import numpy as np
from CoolProp.CoolProp import PropsSI as prop
from data import nodes, blocks
from modules import comp, heat, turb, cond,comb_stoic
from scipy.optimize import root_scalar

#def Sensetivity(P,T):
P2 = 3e6#P
P6 = 1e5
T5 = 1100+273.15#T+273.15
KPDcomp = 0.8
KPDturb = 0.9

X1 = 'REFPROP::O2[0.45]&N2[0.025]&CO2[0.50]&Ar[0.025]'
T1 = 15+273.15
P1 = 1e5
G1 = 120
H1 = prop("H", "P", P1, "T", T1, X1)
S1 = prop("S", "P", P1, "T", T1, X1)
Q1 = prop("Q", "P", P1, "T", T1, X1)
nodes.loc['1'] = [T1, P1, H1, S1, Q1, G1, X1]
X3 = 'REFPROP::Methane[1]&H2[0]&CO[0]'
T3 = 50+273.15
P3 = 1.2e6
G3 = 1
H3 = prop("H", "P", P3, "T", T3, X3)
S3 = prop("S", "P", P3, "T", T3, X3)
Q3 = prop("Q", "P", P3, "T", T3, X3)
nodes.loc['3'] = [T3, P3, H3, S3, Q3, G3, X3]

def T_CC(G3):
    nodes.loc['3','G'] = float(G3)
    comp('AirCOMP', '1', '2', P2, KPDcomp)
    comp('FuelCOMP', '3', '4', P2, KPDcomp)
    comb_stoic('COMB', '2', '4','5',dP=0)
    print(nodes.loc['5']['T'] - T5)
    return nodes.loc['5']['T'] - T5
root_scalar(T_CC,bracket=[0.2,5], xtol=10**-9,method='bisect')
turb('TURB', '5', '6', P6, KPDturb)

# print(nodes.iloc[:,0:6])
# print(blocks)
Nturb = blocks.loc['TURB','N']
Ncomp = (blocks.loc['AirCOMP','N']+blocks.loc['FuelCOMP','N'])
Q = blocks.loc['COMB','Q']
G = nodes.loc['5']['G']
KPD = (Nturb-Ncomp)/Q*100
print(P2,T5,KPD,(Nturb-Ncomp)/G)
print(nodes.iloc[:,0:5])
# for P in np.linspace(2e6,6e6,6):
#     for T in np.linspace(1060, 1600, 6):
#         Sensetivity(P, T)