import numpy as np
from CoolProp.CoolProp import PropsSI as prop
from data import nodes, blocks
from modules import comp, turb,heat_exch,heat_exch_2streams,split,mix
from scipy.optimize import root_scalar,root
x = 0.68
Pk = 8e6
Tk = 35+273.15
Pcomp = 32e6
T2 = 600+273.15
KPDcomp = 0.8
KPDturb = 0.9
dTh = 10
dTl = 10
def Calc(inp):
    T1 = inp[0]
    T12 = inp[1]
    T5 = inp[2]
    T8 = inp[3]
    G = inp[4]
    X1 = 'CO2'
    P1 = Pcomp
    G1 = 100
    H1 = prop("H", "P", P1, "T", T1, X1)
    S1 = prop("S", "P", P1, "T", T1, X1)
    Q1 = prop("Q", "P", P1, "T", T1, X1)
    nodes.loc['1'] = [T1, P1, H1, S1, Q1, G1, X1]
    heat_exch('HEAT','1','2',T=T2)
    turb('TURB','2','3',Pk,KPDturb)
    X12 = 'CO2'
    P12 = Pcomp
    G12 = 100
    H12 = prop("H", "P", P12, "T", T12, X12)
    S12 = prop("S", "P", P12, "T", T12, X12)
    Q12 = prop("Q", "P", P12, "T", T12, X12)
    nodes.loc['12'] = [T12, P12, H12, S12, Q12, G12, X12]
    heat_exch_2streams('HTHE','3','4','12','1',T22=T1)
    nodes.loc['5'] = nodes.loc['4']
    nodes.loc['5','T'] = T5
    X8 = 'CO2'
    T8 = T8
    P8 = Pcomp
    G8 = G
    H8= prop("H", "P", P8, "T", T8, X8)
    S8 = prop("S", "P", P8, "T", T8, X8)
    Q8 = prop("Q", "P", P8, "T", T8, X8)
    nodes.loc['8'] = [T8, P8, H8, S8, Q8, G8, X8]
    heat_exch_2streams('LTHE','4','5','8','9',T12=T5)
    split('SPLIT','5','6','10',x)
    heat_exch('COND','6','7',T=Tk)
    comp('MCOMP','7','8',Pcomp,KPDcomp)
    comp('ACOMP','10','11',Pcomp,KPDcomp)
    mix('MIX','9','11','12')
    eq3 = blocks.loc['HTHE', 'dT'] - dTh
    eq1 = nodes.loc['8','T'] - T8
    eq2 = nodes.loc['12','T'] - T12
    eq4 = blocks.loc['LTHE','dT'] - dTl
    eq5 = G - G1*x
    return eq3,eq2,eq1,eq4,eq5
root(Calc,x0=[600,500,550,550,50],method='hybr')

print(nodes)
print(blocks.iloc[:,0:2])
N_TURB= blocks.loc['TURB']['N']
N_MCOMP = blocks.loc['MCOMP']['N']
N_RCOMP = blocks.loc['ACOMP']['N']
Q_COND = blocks.loc['COND','Q']
Q_HEAT = blocks.loc['HEAT','Q']
print(Q_HEAT+N_MCOMP+N_RCOMP-N_TURB-Q_COND)