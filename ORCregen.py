import numpy as np
from CoolProp.CoolProp import PropsSI as prop
from data import nodes, blocks
from modules import comp, heat, turb, cond
from scipy.optimize import root_scalar

# Исходные данные блоков
Ppump = 1e6
KPDpump = 0.8
KPDturb = 0.9
dTh = 10
dTr = 5

# Ввод исходных данных в (1)
X1 = 'REFPROP::N2[0.77]&O2[0.14]&H2O[0.06]&CO2[0.03]'
G1 = 500
T1 = 200+273.15
P1 = 1e5
H1 = prop("H", "P", P1, "T", T1, X1)
S1 = prop("S", "P", P1, "T", T1, X1)
Q1 = prop("Q", "P", P1, "T", T1, X1)
nodes.loc['1'] = [T1, P1, H1, S1, Q1, G1, X1]
T2 = 100+273.15


heat_







# Ввод исходных данных в (3)
X3 = 'PR::R236ea'
Pcond = 245000

G3 = G1*(prop("H","T",T1,"P",P1,X1)-prop("H","T",T2,"P",P1,X1))/(prop("H","Q",1,"P",Ppump,X3)-prop("H","Q",0,"P",Pcond,X3))


T3 = prop("T", "P", Pcond, "Q", 0, X3)
H3 = prop("H", "P", Pcond, "Q", 0, X3)
S3 = prop("S", "P", Pcond, "Q", 0, X3)
Q3 = prop("Q", "P", Pcond, "Q", 0, X3)
nodes.loc['3'] = [T3, Pcond, H3, S3, Q3, G3, X3]


def regen_root(T8):
    heat('REGEN', '7', '8', '4', '5', T8)
    return dTr - min(blocks.loc['REGEN', 'T1']-blocks.loc['REGEN', 'T2'])

def Groot(G5):
    nodes.loc['5', 'G'] = G5
    heat('HEAT', '1', '2', '5', '6', T2)
    return dTh - min(blocks.loc['HEAT', 'T1']-blocks.loc['HEAT', 'T2'])

for i in range(9999):
    comp('PUMP', '3', '4', Ppump, KPDpump)
    if i == 0:
        nodes.loc['5'] = nodes.loc['4']
    else:
        root_scalar(regen_root, bracket=[nodes.loc['4', 'T'],nodes.loc['7', 'T']], xtol=10 ** -5)
    print(nodes.loc['5', 'G'])
    root_scalar(Groot, x0=nodes.loc['5', 'G'], xtol=10**-5)

    # fi = 0.5
    # nodes.loc['6', 'G'] = nodes.loc['5', 'G']*fi + nodes.loc['4', 'G']*(1-fi)


    turb('TURB', '6', '7', Pcond, KPDturb)
    root_scalar(regen_root, bracket=[nodes.loc['4', 'T'],nodes.loc['7', 'T']], xtol=10 ** -5)
    cond('COND', '8','3')
    balance = (blocks.loc['HEAT','Q'][-1]+blocks.loc['PUMP','N']-blocks.loc['COND','Q']-blocks.loc['TURB','N'])/blocks.loc['HEAT','Q'][-1]
    print(balance)
    if abs(balance) < 10**-5:
        break

KPD = (blocks.loc['TURB','N']-blocks.loc['PUMP','N'])/blocks.loc['HEAT','Q'][-1]
# print(KPD)
print(nodes.loc[:, "T":"G"])
# print(blocks)

import matplotlib.pyplot as plt


X = nodes.loc["3":"4","S"]
Y = nodes.loc["3":"4","T"]
plt.plot(X,Y)

X = [S for S in np.linspace(nodes.loc["4","S"],nodes.loc["5","S"],50)]
Y = [prop("T","S",S,"P",nodes.loc["5","P"],nodes.loc["5","fluid"]) for S in np.linspace(nodes.loc["4","S"],nodes.loc["5","S"],50)]
plt.plot(X,Y)

X = [S for S in np.linspace(nodes.loc["5","S"],nodes.loc["6","S"],50)]
Y = [prop("T","S",S,"P",nodes.loc["5","P"],nodes.loc["6","fluid"]) for S in np.linspace(nodes.loc["5","S"],nodes.loc["6","S"],50)]
plt.plot(X,Y)

X = nodes.loc["6":"7","S"]
Y = nodes.loc["6":"7","T"]
plt.plot(X,Y)

X = [S for S in np.linspace(nodes.loc["7","S"],nodes.loc["8","S"],50)]
Y = [prop("T","S",S,"P",nodes.loc["7","P"],nodes.loc["8","fluid"]) for S in np.linspace(nodes.loc["7","S"],nodes.loc["8","S"],50)]
plt.plot(X,Y)

X = [S for S in np.linspace(nodes.loc["8","S"],nodes.loc["3","S"],50)]
Y = [prop("T","S",S,"P",nodes.loc["8","P"],nodes.loc["3","fluid"]) for S in np.linspace(nodes.loc["8","S"],nodes.loc["3","S"],50)]
plt.plot(X,Y)
plt.show()


X = blocks.loc["REGEN", "Q"]
Y = blocks.loc["REGEN", "T1"]
plt.plot(X,Y)
X = blocks.loc["REGEN", "Q"]
Y = blocks.loc["REGEN", "T2"]
plt.plot(X,Y)
plt.show()

X = blocks.loc["HEAT", "Q"]
Y = blocks.loc["HEAT", "T1"]
plt.plot(X,Y)
X = blocks.loc["HEAT", "Q"]
Y = blocks.loc["HEAT", "T2"]
plt.plot(X,Y)
plt.show()