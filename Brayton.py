import numpy as np
from CoolProp.CoolProp import PropsSI as prop
from data import nodes, blocks
from modules import comp, heat, turb, cond,comb_stoic_CH4_O2N2,comb_stoic_CH4_O2N22
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

    X1 = 'REFPROP::N2[0.79]&O2[0.21]'
    T1 = 15+273.15
    P1 = 1e5
    G1 = 980
    H1 = prop("H", "P", P1, "T", T1, X1)
    S1 = prop("S", "P", P1, "T", T1, X1)
    Q1 = prop("Q", "P", P1, "T", T1, X1)
    nodes.loc['1'] = [T1, P1, H1, S1, Q1, G1, X1]
    def T_CC(G2):
        X2 = 'REFPROP::Methane'
        T2 = 50+273.15
        P2 = 1.2e6
        G2 = float(G2)
        H2 = prop("H", "P", P2, "T", T2, X2)
        S2 = prop("S", "P", P2, "T", T2, X2)
        Q2 = prop("Q", "P", P2, "T", T2, X2)
        nodes.loc['2'] = [T2, P2, H2, S2, Q2, G2, X2]
        comp('AirCOMP', '1', '3', P3, KPDcomp)
        comp('FuelCOMP', '2', '4', P3, KPDcomp)
        comb_stoic_CH4_O2N22('COMB', '3', '4','5',dP=0)
        return nodes.loc['5']['T'] - T6
    root_scalar(T_CC,x0=G1*0.21*M_O2/(0.79*M_N2+0.21*M_O2)/(4*M_O2/M_CH4), xtol=10**-5)
    turb('TURB', '5', '6', P6, KPDturb)

    # print(nodes.iloc[:,0:6])
    # print(blocks)
    KPD = (blocks.loc['TURB','N']-blocks.loc['AirCOMP','N']-blocks.loc['FuelCOMP','N'])/(nodes.loc['2']['G']*50115000)*100
    KPD2 = (blocks.loc['TURB','N']*0.99-blocks.loc['AirCOMP','N']/0.99-blocks.loc['FuelCOMP','N']/0.99)/(nodes.loc['2']['G']*50115000)*100
    print(P3,T6,KPD2,blocks.loc['TURB','N']*0.99,blocks.loc['AirCOMP','N']/0.99,blocks.loc['FuelCOMP','N']/0.99,nodes.loc['1']['G'],nodes.loc['2']['G'])

for P in np.linspace(2e6,10e6,10):
    for T in np.linspace(1060, 1600, 10):
        Optimize(P, T)