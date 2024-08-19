import numpy as np
from CoolProp.CoolProp import PropsSI as prop
from data import nodes, blocks
from modules import Node, comp, turb,heat_exch
from scipy.optimize import root_scalar
def Sensitivity(P1,T2):
    #T2 = T+273.15 #560+273.15#14e6
    KPDturb = 0.85
    KPDpump = 0.8
    P3 = 5e3
    G = 1
    P2 = 10e6
    T2 = 500+273.15
    Node('2',G=G,P=P2,T=T2,fluid='Water')
    turb('TURB', '2', '3', P3, eff=KPDturb)
    heat_exch('COND', '3', '4', x=0)
    comp('PUMP', '4', '1', P1, eff=KPDpump)
    heat_exch('HEAT', '1', '2',T=T2)

    N_TURB = blocks.loc['TURB','N']
    N_PUMP = blocks.loc['PUMP','N']
    Q_COND = blocks.loc['COND','Q']
    Q_HEAT = blocks.loc['HEAT','Q']
    print('Тепловой баланс:',Q_HEAT+N_PUMP-N_TURB-Q_COND)
    print(nodes)
    print(blocks)
    KPD = (N_TURB-N_PUMP)/Q_HEAT*100
    print('КПД',KPD)
    print(P1,T2,KPD,nodes.loc['3']['Q'])

Sensitivity(10e6, 500)
# for P in np.arange(10e6,30.1e6,2.5e6):
#     for T in np.arange(400,651,50):
#         Sensitivity(P, T)




KPDturb = 0.85


# import matplotlib.pyplot as plt
#
#
# X = [S for S in np.linspace(nodes.loc["1","S"],nodes.loc["2","S"],50)]
# Y = [prop("T","S",S,"P",nodes.loc["1","P"],nodes.loc["2","fluid"])-273.15 for S in np.linspace(nodes.loc["1","S"],nodes.loc["2","S"],50)]
# plt.plot(X,Y)
# print(*[x/1000 for x in X])
# print(*Y)
# P = np.linspace(nodes.loc["2","P"],nodes.loc["3","P"],50)
# P1 = P[-2]
# P = np.append(P[0:-2],np.linspace(P1,nodes.loc["3","P"],100))
# H = nodes.loc["2","H"] - (nodes.loc["2","H"] - prop('H','P',P,'S',nodes.loc["2",'S'],'Water'))*KPDturb
# X = prop("S","P",P,"H",H,nodes.loc["2","fluid"])
# Y = prop("T","P",P,"H",H,nodes.loc["2","fluid"])-273.15
# plt.plot(X,Y)
# print(*[x/1000 for x in X])
# print(*Y)
# X = nodes.loc["3":"4","S"]
# Y = nodes.loc["3":"4","T"]-273.15
# plt.plot(X,Y)
# print(*[x/1000 for x in X])
# print(*Y)
# X = nodes.loc["4":"1","S"]
# Y = nodes.loc["4":"1","T"]-273.15
# plt.plot(X,Y)
# print(*[x/1000 for x in X])
# print(*Y)
# T = np.linspace(273.15,prop('Tcrit',nodes.loc["3", "fluid"]), 50)
# T1 = T[-2]
# T = np.append(T[0:-2],np.linspace(T1,prop('Tcrit',nodes.loc["3", "fluid"]),100))
# Y = T -273.15
# X = prop("S", "T", T, "Q", 0, nodes.loc["2", "fluid"])/1000
#
#
#
# print(*X)
# print(*Y)



# plt.show()