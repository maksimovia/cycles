import numpy as np
from CoolProp.CoolProp import PropsSI as prop
from data import nodes, blocks
from modules import comp, heat, turb, cond,comb_stoic
from scipy.optimize import root_scalar

def Sensetivity(P2,T):
    P6 = 1e5
    T5 = T+273.15#1000+273.15#T+273.15
    KPDcomp = 0.8
    KPDturb = 0.9

    X1 = 'REFPROP::O2[0.20]&N2[0.77]&CO2[0.01]&Ar[0.01]&H2O[0.01]'
    T1 = 15+273.15
    P1 = 1e5
    G1 = 100
    H1 = prop("H", "P", P1, "T", T1, X1)
    S1 = prop("S", "P", P1, "T", T1, X1)
    Q1 = prop("Q", "P", P1, "T", T1, X1)
    nodes.loc['1'] = [T1, P1, H1, S1, Q1, G1, X1]
    X3 = 'REFPROP::Methane[1]&H2[0]&CO[0]'
    T3 = 50+273.15
    P3 = 1.2e5
    G3 = ''
    H3 = prop("H", "P", P3, "T", T3, X3)
    S3 = prop("S", "P", P3, "T", T3, X3)
    Q3 = prop("Q", "P", P3, "T", T3, X3)
    nodes.loc['3'] = [T3, P3, H3, S3, Q3, G3, X3]

    def T_CC(G3):
        nodes.loc['3','G'] = float(G3)
        comp('AirCOMP', '1', '2', P2, KPDcomp)
        comp('FuelCOMP', '3', '4', P2, KPDcomp)
        comb_stoic('COMB', '2', '4','5',dP=0)
        #print(nodes.loc['5']['T'] - T5)
        #print(G3)
        return nodes.loc['5']['T'] - T5
    root_scalar(T_CC,bracket=[0.2,3], xtol=10**-9,method='bisect')
    turb('TURB', '5', '6', P6, KPDturb)

    # print(nodes.iloc[:,0:6])
    # print(blocks)
    Nturb = blocks.loc['TURB','N']
    Ncomp = (blocks.loc['AirCOMP','N']+blocks.loc['FuelCOMP','N'])
    Q = blocks.loc['COMB','Q']
    G = nodes.loc['5']['G']
    KPD = (Nturb-Ncomp)/Q*100
    #print(nodes.iloc[:, 5])
    print(P2,T5-273.15,KPD,(Nturb-Ncomp)/G,Nturb/G)
    P = np.linspace(nodes.loc["1","P"],nodes.loc["2","P"],50)
    H = nodes.loc["1","H"] + (prop('H','P',P,'S',nodes.loc["1",'S'],nodes.loc["2", "fluid"]) - nodes.loc["1","H"])/KPDcomp
    X = prop("S","P",P,"H",H,nodes.loc["1","fluid"])/1000
    Y = prop("T","P",P,"H",H,nodes.loc["1","fluid"])-273.15

    X = [S/1000 for S in np.linspace(nodes.loc["2", "S"], nodes.loc["5", "S"], 50)]
    Y = [prop("T", "S", S, "P", nodes.loc["2", "P"], nodes.loc["5", "fluid"])-273.15 for S in np.linspace(nodes.loc["2", "S"], nodes.loc["5", "S"], 50)]
    P = np.linspace(nodes.loc["5", "P"], nodes.loc["6", "P"], 50)
    H = nodes.loc["5", "H"] - (
                nodes.loc["5", "H"] - prop('H', 'P', P, 'S', nodes.loc["5", 'S'], nodes.loc["6", "fluid"])) * KPDturb
    X = prop("S","P",P,"H",H,nodes.loc["5","fluid"])/1000
    Y = prop("T","P",P,"H",H,nodes.loc["5","fluid"])-273.15

    X = [S / 1000 for S in np.linspace(nodes.loc["6", "S"], nodes.loc["1", "S"], 50)]
    Y = [prop("T", "S", S, "P", nodes.loc["6", "P"], nodes.loc["6", "fluid"]) - 273.15 for S in
         np.linspace(nodes.loc["6", "S"], nodes.loc["1", "S"], 50)]


    print(*X)
    print(*Y)

Sensetivity(2.5e6, 1200)
# print(nodes.iloc[:,0:6])
# print(blocks)
# for P in np.linspace(0.4e6,6e6,30):
#     Sensetivity(P, 1200)
    #for T in np.linspace(1060, 1600, 6):
