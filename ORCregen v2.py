import numpy as np
from CoolProp.CoolProp import PropsSI as prop
from data import nodes, blocks
from modules import comp, heat, turb, cond,heat_exch_2streams,heat_exch
from scipy.optimize import root_scalar,root

def Optimize(P6,Pk):
    # Исходные данные блоков
    #P6 = 6e6
    KPDpump = 0.8
    KPDturb = 0.9
    dTh = 10
    dTr = 5

    def Calc(inp):
        Gorc = inp[0]
        T1 = inp[1]

        X1 = 'R236ea'
        P1 = P6
        H1 = prop("H", "P", P1, "T", T1, X1)
        S1 = prop("S", "P", P1, "T", T1, X1)
        Q1 = prop("Q", "P", P1, "T", T1, X1)
        nodes.loc['1'] = [T1, P1, H1, S1, Q1, Gorc, X1]
        heat_exch_2streams('HEAT','7','8','1','2',T12=T8)
        turb('TURB','2','3',Pk,KPDturb)
        T6 = inp[2]
        X6 = 'R236ea'
        H6 = prop("H", "P", P6, "T", T6, X6)
        S6 = prop("S", "P", P6, "T", T6, X6)
        Q6 = prop("Q", "P", P6, "T", T6, X6)
        nodes.loc['6'] = [T6, P6, H6, S6, Q6, Gorc, X6]
        heat_exch_2streams('REGEN','3','4','6','1',T22=T1)
        heat_exch('COND','4','5',x=0)
        comp('PUMP','5','6',P6,KPDpump)

        eq1 = T6 - nodes.loc['6','T']
        eq2 = blocks.loc['HEAT','dT'] - dTh
        eq3 = blocks.loc['REGEN','dT'] - dTr
        #print(eq1,eq2,eq3)
        return eq1,eq2,eq3

    X7 = 'REFPROP::N2[0.77]&O2[0.14]&H2O[0.06]&CO2[0.03]'
    G7 = 500
    T7 = 200+273.15
    P7 = 1e5
    H7 = prop("H", "P", P7, "T", T7, X7)
    S7 = prop("S", "P", P7, "T", T7, X7)
    Q7 = prop("Q", "P", P7, "T", T7, X7)
    nodes.loc['7'] = [T7, P7, H7, S7, Q7, G7, X7]
    T8 = 80+273.15
    #Pk = 245000

    root(Calc,x0 =[300,300,350],method='hybr')
    # print(nodes.iloc[:,0:6])
    # print(blocks)

    Qheat = blocks.loc['HEAT','Q']
    Qcond = blocks.loc['COND','Q']
    Nturb = blocks.loc['TURB','N']
    Npump = blocks.loc['PUMP','N']
    KPD1 = (Nturb-Npump)/Qheat
    print(P6,Pk,KPD1)
    # print(Qheat-Qcond-Nturb+Npump,'Balance')
    # print("KPD1",(Nturb-Npump)/Qheat)
    # print("KPD2",1-(Qcond)/Qheat)
    # # print(blocks.loc['REGEN','T1'])
    # print(blocks.loc['REGEN','T2'])
    # print([blocks.loc['REGEN','Q']/20*i for i in range(21)])

for P0 in np.linspace(3e6,8e6,6):
    for Pk in np.linspace(0.1e6,0.6e6,6):
        Optimize(P0,Pk)