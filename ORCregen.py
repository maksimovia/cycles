from CoolProp.CoolProp import PropsSI as prop
from data import nodes, blocks
from modules import comp, heat, turb, cond
from scipy.optimize import root_scalar

# Исходные данные блоков
Ppump = 6e6
KPDpump = 0.8
KPDturb = 0.9
dTh = 20
dTr = 5

# Ввод исходных данных в (3)
X3 = 'R236ea'
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


for i in range(9999):
    comp('PUMP', '3', '4', Ppump, KPDpump)
    if i == 0:
        nodes.loc['5'] = nodes.loc['4']
    else:
        def Groot1(T8):
            heat('REGEN', '7', '8', '4', '5', T8)
            return dTr - blocks.loc['REGEN', 'dT']
        root_scalar(Groot1, bracket=[nodes.loc['4', 'T'], nodes.loc['7', 'T']], xtol=10 ** -5)

    def Groot2(G5):
        nodes.loc['5', 'G'] = G5
        heat('HEAT', '1', '2', '5', '6', Tout)
        return dTh - blocks.loc['HEAT', 'dT']
    root_scalar(Groot2, x0=nodes.loc['5', 'G'], xtol=10**-5)

    turb('TURB', '6', '7', Pcond, KPDturb)

    def Groot3(T8):
        heat('REGEN', '7', '8', '4', '5', T8)
        return dTr - blocks.loc['REGEN', 'dT']
    root_scalar(Groot3, bracket=[nodes.loc['4', 'T'], nodes.loc['7', 'T']], xtol=10**-5)

    cond('COND', '8','3')

    balance = abs(blocks.loc['HEAT','Q']+blocks.loc['PUMP','N']-blocks.loc['COND','Q']-blocks.loc['TURB','N'])/blocks.loc['HEAT','Q']
    print(balance)
    if balance < 10**-5:
        break

KPD = (blocks.loc['TURB','N']-blocks.loc['PUMP','N'])/blocks.loc['HEAT','Q']
print(KPD)
print(nodes.loc[:, "T":"G"])


