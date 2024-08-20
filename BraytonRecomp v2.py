import numpy as np
from CoolProp.CoolProp import PropsSI as prop
from data import nodes, blocks
from modules import Node, comp, turb,heat_exch,heat_exch_2streams,split,mix
from scipy.optimize import root_scalar,root,minimize

def Sensitivity(inp):
    P1 = inp[0]
    x = inp[1]
    #x = 0.68
    G1 = 1
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
        Node('1',G=G1,P=P1,T=T1,fluid='CO2')
        heat_exch('HEAT','1','2',T=T2)
        turb('TURB','2','3',P7,KPDturb)
        T12 = input[1]
        Node('12',G=G1,P=P1,T=T12,fluid='CO2')
        heat_exch_2streams('HTHE','3','4','12','1',T22=T1)
        T5 = input[2]
        Node('5',G=G1,P=P7,T=T5,fluid='CO2')
        T8 = input[3]
        Node('8',G=G1*x,P=P1,T=T8,fluid='CO2')
        heat_exch_2streams('LTHE','4','5','8','9',T12=T5)
        split('SPLIT','5','6','10',x)
        heat_exch('COOL','6','7',T=T7)
        comp('MCOMP','7','8',P1,KPDcomp)
        comp('ACOMP','10','11',P1,KPDcomp)
        mix('MIX','9','11','12')
        Equation1 = nodes.loc['8','T'] - T8
        Equation2 = nodes.loc['12','T'] - T12
        Equation3 = blocks.loc['HTHE', 'dT'] - dTh
        Equation4 = blocks.loc['LTHE','dT'] - dTl
        return Equation3,Equation2,Equation1,Equation4
    root(Calc,x0=[800, 500, 400, 350],method='hybr')
    N_TURB = blocks.loc['TURB']['N']
    N_MCOMP = blocks.loc['MCOMP']['N']
    N_RCOMP = blocks.loc['ACOMP']['N']
    Q_COND = blocks.loc['COOL','Q']
    Q_HEAT = blocks.loc['HEAT','Q']
    KPD1 = (N_TURB - N_MCOMP - N_RCOMP)/Q_HEAT*100
    print(Q_HEAT + N_MCOMP + N_RCOMP - N_TURB - Q_COND)
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

Sensitivity([23.5e6,0.67])
print(blocks)
# minimize(Sensitivity,x0=[22e6,0.6], method='Nelder-Mead',tol=10**-2)

# print(nodes)
# print(blocks.iloc[:,0:2])
# for P in np.arange(18e6,28.1e6,0.5e6):
#     for x in np.arange(0.5,0.78,0.01):
#         Sensitivity(P,x)