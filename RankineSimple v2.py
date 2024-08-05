import numpy as np
from CoolProp.CoolProp import PropsSI as prop
from data import nodes, blocks
from modules import comp, turb,heat_exch
from scipy.optimize import root_scalar
def Sensitivity(P,T):
    T2 = T+273.15 #560+273.15
    P = P #14e6
    KPDturb = 0.85
    KPDpump = 0.8
    Pcond = 5e3
    G = 100
    def Calc(T1in):
        X1 = 'Water'
        P1 = P
        T1 = T1in
        H1 = prop("H", "P", P1, "T", T1, X1)
        S1 = prop("S", "P", P1, "T", T1, X1)
        Q1 = prop("Q", "P", P1, "T", T1, X1)
        nodes.loc['1'] = [T1, P1, H1, S1, Q1, G, X1]
        heat_exch('BOILER', '1', '2',T=T2)
        turb('TURB', '2', '3', Pcond, eff=KPDturb)
        heat_exch('COND', '3', '4',x=0)
        comp('PUMP','4','1',P1,eff=KPDpump)
        #print(T1in - nodes.loc['1']['T'])
        return T1in - nodes.loc['1']['T']
    root_scalar(Calc,bracket=[273.15,T2],xtol=10**-5,method='bisect')
    N_turb = blocks.loc['TURB','N']
    N_pump = blocks.loc['PUMP','N']
    Q_cond = blocks.loc['COND','Q']
    Q_boiler = blocks.loc['BOILER','Q']
    #print('Тепловой баланс:',Q_boiler+N_pump-N_turb-Q_cond)
    #print(nodes)
    #print(blocks)
    KPD1 = (N_turb-N_pump)/Q_boiler*100
    KPD2 = (1 - Q_cond/Q_boiler)*100
    #print('КПД',KPD1,KPD2)
    print(P,T2,KPD1,nodes.loc['3']['Q'])

for P in np.arange(10e6,30.1e6,2.5e6):
    for T in np.arange(400,651,50):
        Sensitivity(P, T)







#
# import matplotlib.pyplot as plt
#
#
# X = [S for S in np.linspace(nodes.loc["1","S"],nodes.loc["2","S"],50)]
# Y = [prop("T","S",S,"P",nodes.loc["1","P"],nodes.loc["2","fluid"])-273.15 for S in np.linspace(nodes.loc["1","S"],nodes.loc["2","S"],50)]
# plt.plot(X,Y)
#
#
# P = np.linspace(nodes.loc["2","P"],nodes.loc["3","P"],50)
# P1 = P[-2]
# P = np.append(P[0:-2],np.linspace(P1,nodes.loc["3","P"],100))
# H = nodes.loc["2","H"] - (nodes.loc["2","H"] - prop('H','P',P,'S',nodes.loc["2",'S'],'Water'))*KPDturb
# X = prop("S","P",P,"H",H,nodes.loc["2","fluid"])
# Y = prop("T","P",P,"H",H,nodes.loc["2","fluid"])-273.15
# plt.plot(X,Y)
#
# X = nodes.loc["3":"4","S"]
# Y = nodes.loc["3":"4","T"]-273.15
# plt.plot(X,Y)
#
# X = nodes.loc["4":"1","S"]
# Y = nodes.loc["4":"1","T"]-273.15
# plt.plot(X,Y)
#
# #plt.show()