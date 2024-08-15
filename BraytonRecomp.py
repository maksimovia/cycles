import numpy as np
from CoolProp.CoolProp import PropsSI as prop
from data import nodes, blocks
from modules import comp, turb,heat_exch,heat_exch_2streams,split,mix
from scipy.optimize import root_scalar,root,minimize

def Sensitivity(inp):
    P1 = inp[0]
    x = inp[1]
    #x = 0.68
    G1 = 100
    P7 = 8e6
    T7 = 35+273.15
    #P1 = 22e6
    T2 = 600+273.15
    KPDcomp = 0.8
    KPDturb = 0.9
    dTh = 10
    dTl = 10
    def Calc(input):
        T1 = input[0]
        X1 = 'CO2'
        H1 = prop("H", "P", P1, "T", T1, X1)
        S1 = prop("S", "P", P1, "T", T1, X1)
        Q1 = prop("Q", "P", P1, "T", T1, X1)
        nodes.loc['1'] = [T1, P1, H1, S1, Q1, G1, X1]
        heat_exch('HEAT','1','2',T=T2)
        turb('TURB','2','3',P7,KPDturb)
        T12 = input[1]
        X12 = 'CO2'
        P12 = P1
        G12 = G1
        H12 = prop("H", "P", P12, "T", T12, X12)
        S12 = prop("S", "P", P12, "T", T12, X12)
        Q12 = prop("Q", "P", P12, "T", T12, X12)
        nodes.loc['12'] = [T12, P12, H12, S12, Q12, G12, X12]
        heat_exch_2streams('HTHE','3','4','12','1',T22=T1)
        T5 = input[2]
        X5 = 'CO2'
        P5 = P7
        G5 = G1
        H5 = prop("H", "P", P5, "T", T5, X5)
        S5 = prop("S", "P", P5, "T", T5, X5)
        Q5 = prop("Q", "P", P5, "T", T5, X5)
        nodes.loc['5'] = [T5, P5, H5, S5, Q5, G5, X5]
        T8 = input[3]
        X8 = 'CO2'
        T8 = T8
        P8 = P1
        G8 = G1*x
        H8= prop("H", "P", P8, "T", T8, X8)
        S8 = prop("S", "P", P8, "T", T8, X8)
        Q8 = prop("Q", "P", P8, "T", T8, X8)
        nodes.loc['8'] = [T8, P8, H8, S8, Q8, G8, X8]
        heat_exch_2streams('LTHE','4','5','8','9',T12=T5)
        split('SPLIT','5','6','10',x)
        heat_exch('COOL','6','7',T=T7)
        comp('MCOMP','7','8',P1,KPDcomp)
        comp('ACOMP','10','11',P1,KPDcomp)
        mix('MIX','9','11','12')

        eq1 = nodes.loc['8','T'] - T8
        eq2 = nodes.loc['12','T'] - T12
        eq3 = blocks.loc['HTHE', 'dT'] - dTh
        eq4 = blocks.loc['LTHE','dT'] - dTl
        #print(input)
        return eq3,eq2,eq1,eq4
    root(Calc,x0=[800, 500, 400, 350],method='hybr')

    # print(nodes)
    # print(blocks.iloc[:,0:2])
    N_TURB = blocks.loc['TURB']['N']
    N_MCOMP = blocks.loc['MCOMP']['N']
    N_RCOMP = blocks.loc['ACOMP']['N']
    Q_COND = blocks.loc['COOL','Q']
    Q_HEAT = blocks.loc['HEAT','Q']
    KPD1 = (N_TURB - N_MCOMP - N_RCOMP)/Q_HEAT*100
    KPD2 = 1 - Q_COND/Q_HEAT
    print(P1/1e6,x,KPD1)
    # X = [S/1000 for S in np.linspace(nodes.loc["6", "S"], nodes.loc["7", "S"], 50)]
    # Y = [prop("T", "S", S, "P", nodes.loc["6", "P"], nodes.loc["6", "fluid"])-273.15 for S in np.linspace(nodes.loc["6", "S"], nodes.loc["7", "S"], 50)]
    # P = np.linspace(nodes.loc["10","P"],nodes.loc["11","P"],50)
    # H = nodes.loc["10","H"] + (prop('H','P',P,'S',nodes.loc["10",'S'],nodes.loc["11", "fluid"]) - nodes.loc["10","H"])/KPDcomp
    # X = prop("S","P",P,"H",H,nodes.loc["7","fluid"])/1000
    # Y = prop("T","P",P,"H",H,nodes.loc["7","fluid"])-273.15
    # T = np.linspace(273.15,prop('Tcrit',nodes.loc["3", "fluid"]), 50)
    # T1 = T[-2]
    # T = np.append(T[0:-2],np.linspace(T1,prop('Tcrit',nodes.loc["3", "fluid"]),100))
    # Y = T -273.15
    # X = prop("S", "T", T, "Q", 0, nodes.loc["2", "fluid"])/1000
    # print(*X)
    # print(*Y)

    # print(*[blocks.loc['HTHE','Q']/20*i/1e6 for i in range(21)]) #+ blocks.loc['HTHE', 'Q']/1e6
    # print(*[x - 273.15 for x in blocks.loc['HTHE', 'T2']])
    return -KPD1
#print(Q_HEAT+N_MCOMP+N_RCOMP-N_TURB-Q_COND)
#Sensitivity(23.5e6,0.67)
minimize(Sensitivity,x0=[22e6,0.6], method='Nelder-Mead',tol=10**-2)

# print(nodes)
# print(blocks.iloc[:,0:2])
# for P in np.arange(18e6,28.1e6,0.5e6):
#     for x in np.arange(0.5,0.78,0.01):
#         Sensitivity(P,x)