import numpy as np
from CoolProp.CoolProp import PropsSI as prop
from data import nodes, blocks
from modules import comp, heat, turb, cond
from scipy.optimize import root_scalar



def cycle_calc(Ppump):
    # Исходные данные блоков
    KPDpump = 0.8
    KPDturb = 0.9
    dTh = 10
    dTr = 5

    # Ввод исходных данных в (3)
    X3 = 'REFPROP::R236ea'
    G3 = 300
    Pcond = 245000
    T3 = prop("T", "P", Pcond, "Q", 0, X3)
    H3 = prop("H", "P", Pcond, "Q", 0, X3)
    S3 = prop("S", "P", Pcond, "Q", 0, X3)
    Q3 = prop("Q", "P", Pcond, "Q", 0, X3)
    nodes.loc['3'] = [T3, Pcond, H3, S3, Q3, G3, X3]

    # Ввод исходных данных в (1)
    X1 = 'REFPROP::N2[0.77]&O2[0.14]&H2O[0.06]&CO2[0.03]'
    G1 = 500
    T1 = 200+273.15
    P1 = 1e5
    H1 = prop("H", "P", P1, "T", T1, X1)
    S1 = prop("S", "P", P1, "T", T1, X1)
    Q1 = prop("Q", "P", P1, "T", T1, X1)
    nodes.loc['1'] = [T1, P1, H1, S1, Q1, G1, X1]
    Tout = 80+273.15


    def regen_root(T8):
        heat('REGEN', '7', '8', '4', '5', T8)
        return dTr - min(blocks.loc['REGEN', 'T1']-blocks.loc['REGEN', 'T2'])

    def Groot(G5):
        nodes.loc['5', 'G'] = G5
        heat('HEAT', '1', '2', '5', '6', Tout)
        return dTh - min(blocks.loc['HEAT', 'T1']-blocks.loc['HEAT', 'T2'])

    for i in range(9999):
        comp('PUMP', '3', '4', Ppump, KPDpump)
        if i == 0:
            nodes.loc['5'] = nodes.loc['4']
        else:
            root_scalar(regen_root, bracket=[nodes.loc['7', 'T'], nodes.loc['4', 'T']], xtol=10 ** -5)

        root_scalar(Groot, x0=nodes.loc['5', 'G'], xtol=10**-5)
        turb('TURB', '6', '7', Pcond, KPDturb)
        root_scalar(regen_root, bracket=[nodes.loc['7', 'T'], nodes.loc['4', 'T']], xtol=10 ** -5)
        cond('COND', '8','3')
        balance = (blocks.loc['HEAT','Q'][-1]+blocks.loc['PUMP','N']-blocks.loc['COND','Q']-blocks.loc['TURB','N'])/blocks.loc['HEAT','Q'][-1]
        print(balance)
        if abs(balance) < 10**-5:
            break

    KPD = (blocks.loc['TURB','N']-blocks.loc['PUMP','N'])/blocks.loc['HEAT','Q'][-1]
    return KPD

KPDlist = []
for Ppump in np.linspace(2e6,6e6,10):
    KPDlist.append(cycle_calc(Ppump))

print(KPDlist)
