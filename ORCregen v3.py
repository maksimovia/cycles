import numpy as np
from CoolProp.CoolProp import PropsSI as prop
from data import nodes, blocks
from modules import Node, comp, heat, turb, cond,heat_exch_2streams,heat_exch
from scipy.optimize import root_scalar,root
import pandas as pd
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
        heat_exch_2streams('HEAT','7','8','1','2',T12=T8)
        turb('TURB','2','3',P3,KPDturb)
        T6 = inp[2]
        Node('6', G=Gorc, P=P6, T=T6, fluid='R236ea')
        heat_exch_2streams('REGEN','3','4','6','1',T22=T1)
        heat_exch('COND','4','5',x=0)
        comp('PUMP','5','6',P6,KPDpump)
        Equation1 = T6 - nodes.loc['6','T']
        Equation2 = blocks.loc['HEAT','dT'] - dTh
        Equation3 = blocks.loc['REGEN','dT'] - dTr
        return Equation1,Equation2,Equation3

    T7 = 200+273.15
    P7 = 1e5
    Node('7',G=500,P=P7,T=T7,fluid='REFPROP::N2[0.77]&O2[0.14]&H2O[0.06]&CO2[0.03]')
    T8 = 80+273.15
    root(Calc,x0 =[300,300,350],method='hybr')
    # print(nodes.iloc[:,0:6])
    # print(blocks)

    Qheat = blocks.loc['HEAT','Q']
    Qcond = blocks.loc['COND','Q']
    Nturb = blocks.loc['TURB','N']
    Npump = blocks.loc['PUMP','N']
    KPD = (Nturb-Npump)/Qheat*100
    print(P6,P3,KPD)
    # print(Qheat-Qcond-Nturb+Npump,'Balance')
    # print("KPD1",(Nturb-Npump)/Qheat)
    # print("KPD2",1-(Qcond)/Qheat)
    # # print(blocks.loc['REGEN','T1'])
    # print(blocks.loc['REGEN','T2'])
    # print([blocks.loc['REGEN','Q']/20*i for i in range(21)])

    # X = [S/1000 for S in np.linspace(nodes.loc["4", "S"], nodes.loc["5", "S"], 50)]
    # Y = [prop("T", "S", S, "P", nodes.loc["4", "P"], nodes.loc["4", "fluid"])-273.15 for S in np.linspace(nodes.loc["4", "S"], nodes.loc["5", "S"], 50)]
    #
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
    #
    #
    #
    # print(*X)
    # print(*Y)

for P0 in np.linspace(3e6,8e6,6):
    for Pk in np.linspace(0.1e6,0.6e6,6):
        Sensitivity(6e6, 0.245e6)
