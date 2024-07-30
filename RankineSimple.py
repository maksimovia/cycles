import numpy as np
from CoolProp.CoolProp import PropsSI as prop
from data import nodes, blocks
from modules import comp, turb, cond,heater
from scipy.optimize import root_scalar
def Optimize(P,T):
    KPDturb = 0.9
    KPDpump = 0.8
    P2 = 5e3
    dP = 1e6

    X1 = 'Water'
    T1 = T+273.15
    P1 = P
    G1 = 100
    H1 = prop("H", "P", P1, "T", T1, X1)
    S1 = prop("S", "P", P1, "T", T1, X1)
    Q1 = prop("Q", "P", P1, "T", T1, X1)
    nodes.loc['1'] = [T1, P1, H1, S1, Q1, G1, X1]
    turb('TURB', '1', '2', P2, eff=KPDturb)
    heater('COND', '2', '3',x=0)
    comp('PUMP','3','4',P1+dP,eff=KPDpump)
    def Qroot(Qin):
        heater('BOILER','4','5',dP=dP,Q=Qin)
        dh=nodes.loc['1']['H']-nodes.loc['5']['H']
        return dh
    root_scalar(Qroot,x0=100000, xtol=10**-5)
    Nturb = blocks.loc['TURB','N']
    Npump = blocks.loc['PUMP','N']
    Q = blocks.loc['BOILER','Q']
    KPD = (Nturb-Npump)/Q*100
    print(P,T,KPD,nodes.loc['2']['Q'])

for P in np.linspace(5e6,22.5e6,50):
    for T in np.linspace(300,700,10):
        Optimize(P, T)