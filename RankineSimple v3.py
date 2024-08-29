import numpy as np
from CoolProp.CoolProp import PropsSI as prop
from data import nodes, blocks
from modules import Node, comp, turb,heat_source
from scipy.optimize import root_scalar
def Sensitivity(P1,T1):

    P2 = 5e3
    G1 = 1
    #P1 = 14e6
    #T1 = 560+273.15
    KPDturb = 0.85
    KPDpump = 0.8
    Node('1',G=G1,P=P1,T=T1,fluid='Water')
    turb('TURB', '1', '2', P2, eff=KPDturb)
    heat_source('COND', '2', '3', x=0)
    comp('PUMP', '3', '4', P1, eff=KPDpump)
    heat_source('HEAT', '4', '1', T=T1)

    N_TURB = blocks.loc['TURB','N']
    N_PUMP = blocks.loc['PUMP','N']
    Q_COND = blocks.loc['COND','Q']
    Q_HEAT = blocks.loc['HEAT','Q']
    # print('Тепловой баланс:',Q_HEAT+N_PUMP-N_TURB-Q_COND)
    # print(nodes)
    # print(blocks)
    KPD = (N_TURB-N_PUMP)/Q_HEAT*100
    # print('КПД',KPD)
    print(P1/1e6,T1-273.15,KPD,nodes.loc['3']['Q'],nodes.loc['1']['S'])

Sensitivity(5e6, 550+273.15)
# for P in np.arange(5e6,25.1e6,1e6):
#     for T in np.arange(300,651,50):
#         Sensitivity(P, T+273.15)

KPDturb = 0.85

import matplotlib.pyplot as plt
X = [S for S in np.linspace(nodes.loc["4","S"],nodes.loc["1","S"],50)]
Y = [prop("T","S",S,"P",nodes.loc["4","P"],nodes.loc["4","fluid"])-273.15 for S in np.linspace(nodes.loc["4","S"],nodes.loc["1","S"],50)]
plt.plot(X,Y)
print(*[x/1000 for x in X])
print(*Y)
P = np.linspace(nodes.loc["1","P"],nodes.loc["2","P"],50)
P1 = P[-2]
P = np.append(P[0:-2],np.linspace(P1,nodes.loc["2","P"],100))
H = nodes.loc["1","H"] - (nodes.loc["1","H"] - prop('H','P',P,'S',nodes.loc["1",'S'],'Water'))*KPDturb
X = prop("S","P",P,"H",H,nodes.loc["2","fluid"])
Y = prop("T","P",P,"H",H,nodes.loc["2","fluid"])-273.15
plt.plot(X,Y)
print(*[x/1000 for x in X])
print(*Y)
X = nodes.loc["2":"3","S"]
Y = nodes.loc["2":"3","T"]-273.15
plt.plot(X,Y)
print(*[x/1000 for x in X])
print(*Y)
X = nodes.loc["3":"4","S"]
Y = nodes.loc["3":"4","T"]-273.15
plt.plot(X,Y)
print(*[x/1000 for x in X])
print(*Y)
T = np.linspace(273.15,prop('Tcrit',nodes.loc["3", "fluid"]), 50)
T1 = T[-2]
T = np.append(T[0:-2],np.linspace(T1,prop('Tcrit',nodes.loc["3", "fluid"]),100))
Y = T-273.15
X = prop("S", "T", T, "Q", 0, nodes.loc["2", "fluid"])
plt.plot(X,Y)
print(*[x/1000 for x in X])
print(*Y)
Y = T-273.15
X = prop("S", "T", T, "Q", 1, nodes.loc["2", "fluid"])
plt.plot(X,Y)
print(*[x/1000 for x in X])
print(*Y)
plt.show()