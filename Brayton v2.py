import numpy as np
from CoolProp.CoolProp import PropsSI as prop
from data import nodes, blocks
from modules import Node, comp, turb,comb_stoic2
from scipy.optimize import root_scalar
def Sensetivity(P2,T):
    P5 = 1e5
    T4 = T+273.15#1000+273.15#T+273.15
    KPDcomp = 0.8
    KPDturb = 0.9

    T1 = 15 + 273.15
    P1 = 1e5
    G1 = 100
    Node('1', G=G1, P=P1, T=T1, fluid='REFPROP::O2[0.21]&N2[0.79]')
    comp('COMP', '1', '2', P2, KPDcomp)
    T3 = 100 + 273.15
    P3 = P2
    Node('3', G=None, P=P3, T=T3, fluid='REFPROP::Methane[1]&H2[0]&CO[0]')
    def Calc(G3):
        nodes.loc['3','G'] = float(G3)
        comb_stoic2('COMB', '2', '3','4',w_N2 = 0.79)
        #print(nodes.loc['4']['T'] - T4)
        return nodes.loc['4']['T'] - T4
    root_scalar(Calc,bracket=[0.2,3], xtol=10**-5,method='bisect')
    turb('TURB', '4', '5', P5, KPDturb)

    # print(nodes.iloc[:,0:6])
    # print(blocks)
    Nturb = blocks.loc['TURB','N']
    Ncomp = blocks.loc['COMP','N']
    Q = blocks.loc['COMB','Q']
    G = nodes.loc['1']['G']
    KPD = (Nturb-Ncomp)/Q*100
    #print(nodes.iloc[:, 5])
    print(P2,T4-273.15,KPD,(Nturb-Ncomp)/G,Nturb/G,nodes.loc['3','G'])
    import matplotlib.pyplot as plt
    P = np.linspace(nodes.loc["1","P"],nodes.loc["2","P"],50)
    H = nodes.loc["1","H"] + (prop('H','P',P,'S',nodes.loc["1",'S'],nodes.loc["2", "fluid"]) - nodes.loc["1","H"])/KPDcomp
    X = prop("S","P",P,"H",H,nodes.loc["1","fluid"])/1000
    Y = prop("T","P",P,"H",H,nodes.loc["1","fluid"])-273.15
    plt.plot(X,Y)
    # print(*X)
    # print(*Y)
    X = [S/1000 for S in np.linspace(nodes.loc["2", "S"], nodes.loc["4", "S"], 50)]
    Y = [prop("T", "S", S, "P", nodes.loc["2", "P"], nodes.loc["4", "fluid"])-273.15 for S in np.linspace(nodes.loc["2", "S"], nodes.loc["4", "S"], 50)]
    Y = [T - 273.15 for T in np.linspace(nodes.loc["2", "T"], nodes.loc["4", "T"], 50)]
    plt.plot(X,Y)
    # print(*X)
    # print(*Y)
    P = np.linspace(nodes.loc["4", "P"], nodes.loc["5", "P"], 50)
    H = nodes.loc["4", "H"] - (
                nodes.loc["4", "H"] - prop('H', 'P', P, 'S', nodes.loc["4", 'S'], nodes.loc["5", "fluid"])) * KPDturb
    X = prop("S","P",P,"H",H,nodes.loc["5","fluid"])/1000
    Y = prop("T","P",P,"H",H,nodes.loc["5","fluid"])-273.15
    plt.plot(X, Y)
    # print(*X)
    # print(*Y)
    X = [S / 1000 for S in np.linspace(nodes.loc["5", "S"], nodes.loc["1", "S"], 50)]
    Y = [prop("T", "S", S, "P", nodes.loc["5", "P"], nodes.loc["5", "fluid"]) - 273.15 for S in
         np.linspace(nodes.loc["5", "S"], nodes.loc["1", "S"], 50)]
    Y = [T - 273.15 for T in np.linspace(nodes.loc["5", "T"], nodes.loc["1", "T"], 50)]
    plt.plot(X, Y)

    X = [S / 1000 for S in np.linspace(nodes.loc["1", "S"], nodes.loc["5", "S"], 50)]
    Y = [prop("T", "S", S, "P", nodes.loc["1", "P"], nodes.loc["1", "fluid"]) - 273.15 for S in
         np.linspace(nodes.loc["1", "S"], nodes.loc["5", "S"], 50)]
    plt.show()
    # print(*X)
    # print(*Y)

# Sensetivity(2e6, 1200)
# print(nodes.iloc[:,0:6])
# print(blocks)
Sensetivity(2.5e6, 1200)
# for P in np.linspace(0.4e6,6e6,30):
#     Sensetivity(P, 1200)

