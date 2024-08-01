import numpy as np
from CoolProp.CoolProp import PropsSI as prop
from data import nodes, blocks
from modules import comp, heat, turb, cond,heat_exch_2streams,heat_exch
from scipy.optimize import root_scalar,root

# Исходные данные блоков
P4 = 6e6
KPDpump = 0.8
KPDturb = 0.9
dTh = 10
dTr = 5


X1 = 'REFPROP::N2[0.77]&O2[0.14]&H2O[0.06]&CO2[0.03]'
G1 = 500
T1 = 200+273.15
P1 = 1e5
H1 = prop("H", "P", P1, "T", T1, X1)
S1 = prop("S", "P", P1, "T", T1, X1)
Q1 = prop("Q", "P", P1, "T", T1, X1)
nodes.loc['1'] = [T1, P1, H1, S1, Q1, G1, X1]
T2 = 100+273.15
Pk = 245000

def Root1(inp):
    Gorc_ = inp[0]
    T7_ = inp[1]
    T8_ = inp[2]
    X8 = 'R236ea'#8
    P8 = Pk
    G8 = Gorc_
    T8 = T8_
    H8 = prop("H", "P", P8, "T", T8, X8)
    S8 = prop("S", "P", P8, "T", T8, X8)
    Q8 = prop("Q", "P", P8, "T", T8, X8)
    nodes.loc['8'] = [T8, P8, H8, S8, Q8, G8, X8]
    heat_exch('COND','8','3',x=0)#3
    print(nodes.loc['3'])
    comp('PUMP','3','4',P4,KPDpump)
    X7 = 'PR::R236ea'#7
    P7 = Pk
    G7 = Gorc_
    T7 = T7_
    H7 = prop("H", "P", P7, "T", T7, X7)
    S7 = prop("S", "P", P7, "T", T7, X7)
    Q7 = prop("Q", "P", P7, "T", T7, X7)
    nodes.loc['7'] = [T7, P7, H7, S7, Q7, G7, X7]
    heat_exch_2streams('RECUP','7','8','4','5',T12=T8)#5
    heat_exch_2streams('EVAP','1','2','5','6',T12=T2)#6
    turb('TURB','6','7',Pk,KPDturb)#7

    eq1 = blocks.loc['EVAP']/(nodes.loc['6']['H']-nodes.loc['5']['H']) #?
    eq2 = T8_ - nodes.loc['8']['T']
    eq3 = T7_ - nodes.loc['7']['T']
    return eq1,eq2,eq3

root(Root1,x0 =([300,450,350]),method='hybr')