import numpy as np
from CoolProp.CoolProp import PropsSI as prop
from data import nodes, blocks
from modules import comp, heat, turb, cond,comb_stoic
from scipy.optimize import root_scalar
M_N2 = prop('M', 'N2') * 1000
M_O2 = prop('M', 'O2') * 1000
M_CH4 = prop('M', 'CH4') * 1000
def Optimize(P,T):
    P3 = P
    P6 = 1e5
    T6 = T+273.15
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
    def T_CC(G2):
        X2 = 'REFPROP::Methane[0]&H2[1]&CO[0]'
        T2 = 50+273.15
        P2 = 1.2e6
        G2 = float(G2)
        H2 = prop("H", "P", P2, "T", T2, X2)
        S2 = prop("S", "P", P2, "T", T2, X2)
        Q2 = prop("Q", "P", P2, "T", T2, X2)
        nodes.loc['2'] = [T2, P2, H2, S2, Q2, G2, X2]
        comp('AirCOMP', '1', '3', P3, KPDcomp)
        comp('FuelCOMP', '2', '4', P3, KPDcomp)
        comb_stoic('COMB', '3', '4','5',dP=0)
        return nodes.loc['5']['T'] - T6
    root_scalar(T_CC,x0=1, xtol=10**-5)
    turb('TURB', '5', '6', P6, KPDturb)

    # print(nodes.iloc[:,0:6])
    # print(blocks)
    Nturb = blocks.loc['TURB','N']*0.98
    Ncomp = (blocks.loc['AirCOMP','N']+blocks.loc['FuelCOMP','N'])/0.98
    Q = blocks.loc['COMB','Q']
    G = nodes.loc['5']['G']
    KPD = (Nturb-Ncomp)/Q*100
    print(P3,T6,KPD,(Nturb-Ncomp)/G)

# for P in np.linspace(2e6,6e6,6):
#     for T in np.linspace(1060, 1600, 6):
#         Optimize(P, T)