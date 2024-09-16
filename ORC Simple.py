import numpy as np
from CoolProp.CoolProp import PropsSI as prop
from data import nodes, blocks
from modules import Node, comp, turb,heat_source,heat_exch
from scipy.optimize import root_scalar
T5 = 200+273.15
P5 = 1e5
G5 = 500
Node('5',G=G5,P=P5,T=T5,fluid='Air')#REFPROP::N2[0.77]&O2[0.14]&H2O[0.06]&CO2[0.03]
T6 = 80+273.15
def Sensitivity(P1,T1):
    def Calc(G1):
        P2 = 0.25e6
        KPDturb = 0.85
        KPDpump = 0.8
        Node('1',G=G1,P=P1,T=T1,fluid='R236ea')
        turb('TURB', '1', '2', P2, eff=KPDturb)
        heat_source('COND', '2', '3', x=0)
        comp('PUMP', '3', '4', P1, eff=KPDpump)
        heat_exch('HEAT','5','6','4','1',T12=T6)
        return T1-nodes.loc['1','T']
    root_scalar(Calc,x0=400)

    N_TURB = blocks.loc['TURB','N']
    N_PUMP = blocks.loc['PUMP','N']
    Q_COND = blocks.loc['COND','Q']
    Q_HEAT = blocks.loc['HEAT','Q']
    print('Тепловой баланс:',Q_HEAT+N_PUMP-N_TURB-Q_COND)
    print(nodes)
    # print(blocks)
    KPD = (N_TURB-N_PUMP)/Q_HEAT*100
    # print('КПД',KPD)
    print(P1/1e6,T1-273.15,KPD,nodes.loc['3']['Q'],nodes.loc['1','G'])

Sensitivity(5e6, 150+273.15)

