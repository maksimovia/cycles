import numpy as np
from CoolProp.CoolProp import PropsSI as prop
from data import nodes, blocks
from modules import Node, comp, turb,heat_exch,heat_source
from scipy.optimize import root_scalar,root
import pandas as pd

T7 = 200+273.15
P7 = 1e5
Node('7',G=500,P=P7,T=T7,fluid='REFPROP::N2[0.77]&O2[0.14]&H2O[0.06]&CO2[0.03]')
T8 = 80+273.15
def Sensitivity(P6,P3):
    #P6 = 6e6
    KPDpump = 0.8
    KPDturb = 0.9
    dTh = 10
    dTr = 5

    def Calc(inp):
        Gorc = inp[0]
        T1 = inp[1]
        Node('1',G=Gorc,P=P6,T=T1,fluid='R236ea')
        heat_exch('HEAT', '7', '8', '1', '2', T12=T8)
        turb('TURB','2','3',P3,KPDturb)
        T6 = inp[2]
        Node('6', G=Gorc, P=P6, T=T6, fluid='R236ea')
        heat_exch('REGEN', '3', '4', '6', '1', T22=T1)
        heat_source('COND', '4', '5', x=0)
        comp('PUMP','5','6',P6,KPDpump)
        Equation1 = T6 - nodes.loc['6','T']
        Equation2 = blocks.loc['HEAT','dT'] - dTh
        Equation3 = blocks.loc['REGEN','dT'] - dTr
        return Equation1,Equation2,Equation3

    root(Calc,x0 =[300,300,350],method='hybr', tol = 10**-3)
    print(nodes.iloc[:,0:6])
    print(blocks.iloc[:,0:1])

    Qheat = blocks.loc['HEAT','Q']
    Qcond = blocks.loc['COND','Q']
    Nturb = blocks.loc['TURB','N']
    Npump = blocks.loc['PUMP','N']
    KPD = (Nturb-Npump)/Qheat*100
    print(P6,P3,KPD)
    # print(Qheat-Qcond-Nturb+Npump,'Balance')
    # print("KPD1",(Nturb-Npump)/Qheat)
    # print("KPD2",1-(Qcond)/Qheat)
    # print(*[blocks.loc['HEAT', 'Q'] / 20 / 1e6 * i for i in range(21)])
    # print(*[x-273.15 for x in blocks.loc['HEAT','T1']])
    # print(*[x-273.15 for x in blocks.loc['HEAT','T2']])


    # X = [S/1000 for S in np.linspace(nodes.loc["1", "S"], nodes.loc["2", "S"], 50)]
    # Y = [prop("T", "S", S, "P", nodes.loc["1", "P"], nodes.loc["1", "fluid"])-273.15 for S in np.linspace(nodes.loc["1", "S"], nodes.loc["2", "S"], 50)]

    # P = np.linspace(nodes.loc["2","P"],nodes.loc["3","P"],50)
    # P1 = P[-2]
    # P = np.append(P[0:-2],np.linspace(P1,nodes.loc["3","P"],100))
    # H = nodes.loc["2","H"] - (nodes.loc["2","H"] - prop('H','P',P,'S',nodes.loc["2",'S'],nodes.loc["3", "fluid"]))*KPDturb
    # X = prop("S","P",P,"H",H,nodes.loc["2","fluid"])/1000
    # Y = prop("T","P",P,"H",H,nodes.loc["2","fluid"])-273.15
    #
    # X = nodes.loc["5":"6", "S"]/1000
    # Y = nodes.loc["5":"6", "T"]-273.15
    #
    # T = np.linspace(273.15,prop('Tcrit',nodes.loc["3", "fluid"]), 50)
    # T1 = T[-2]
    # T = np.append(T[0:-2],np.linspace(T1,prop('Tcrit',nodes.loc["3", "fluid"]),100))
    # Y = T -273.15
    # X = prop("S", "T", T, "Q", 1, nodes.loc["2", "fluid"])/1000
    # print(*X)
    # print(*Y)
Sensitivity(3e6, 0.245e6)
# for P0 in np.linspace(2e6,7e6,10):
#     Sensitivity(P0, 0.35e6)

    # for Pk in np.arange(0.2e6,0.41e6,0.05e6):
